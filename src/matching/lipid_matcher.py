import time
from typing import Set, Union
import logging
from src.parsing.lipid_parser import LipidParser
from src.data.config_loader import ConfigLoader
from src.matching.matching_strategies import MatchingStrategy
from src.models.lipid_components import *
from src.models.models import MatchResult, ChainPosition

logger = logging.getLogger(__name__)


class LipidMatcher:
    def __init__(self, lipid_parser: LipidParser, config_loader: ConfigLoader):
        self.lipid_parser = lipid_parser
        self.matching_strategies = self._load_matching_strategies(config_loader.get_matching_strategies())
        self.super_matching_strategies = self._load_matching_strategies(config_loader.get_super_matching_strategies())
        self.strategy_functions = self._initialize_strategy_functions()
        self.super_strategy_functions = self._initialize_super_strategy_functions()
        self.headgroup_map = None

    def set_headgroup_map(self, headgroup_map: Dict[str, Set[str]]) -> None:
        """Set the headgroup map for looking up lipids by headgroup."""
        self.headgroup_map = headgroup_map

    def _load_matching_strategies(self, strategy_config: Dict[str, List[str]]) -> Dict[MatchingStrategy, List[str]]:
        return {MatchingStrategy[key.upper()]: value for key, value in strategy_config.items()}

    def _initialize_super_strategy_functions(self) -> Dict[MatchingStrategy, Callable]:
        return {
            MatchingStrategy.CER_PC_TO_ACER_LPC: self.match_cer_pc_to_acer_lpc,
            MatchingStrategy.LYSOGPL_GPL_TO_GPL: self.match_lysogpl_gpl_to_gpl,
            MatchingStrategy.PC_PE_TO_NAPE: self.match_pc_pe_to_nape,
            MatchingStrategy.PC_PE_TO_PC: self.match_pc_pe_to_pc,
            MatchingStrategy.PC_PE_TO_PE: self.match_pc_pe_to_pe,
            MatchingStrategy.PCO_SPH_TO_LPCO_CER: self.match_pco_sph_to_lpco_cer,
            # Add more super matching strategies as needed
        }

    def _initialize_strategy_functions(self) -> Dict[MatchingStrategy, Union[Callable, MatchResult]]:
        return {
            MatchingStrategy.COMPLETE: self.match_lipids_completely,
            MatchingStrategy.PC_TO_LPC: self.match_pc_to_lpc_lipids,
            MatchingStrategy.UNMARKED_LPC_TO_PC: self.match_unmarked_lpc_to_pc_using_input_lipids,
            MatchingStrategy.LPC_TO_PC: self.match_lpc_to_pc_using_input_lipids,
            MatchingStrategy.N_ACYL_MATCH: self.match_n_acyl_chain_using_input_lipids,
            MatchingStrategy.BACKBONE_MATCH: self.match_backbone,
            MatchingStrategy.N_ACYL_MATCH_DESATURASE: self.match_n_acyl_chain_for_desaturase,
            MatchingStrategy.NAPE_PA_MATCH: self.match_nape_to_pa,
            MatchingStrategy.NAPE_GPNAE_MATCH: self.match_nape_gpnae,
            MatchingStrategy.NAPE_NAE_MATCH: self.match_nape_nae,
            MatchingStrategy.NAPE_NALPE_MATCH: self.match_nape_nalpe,
            MatchingStrategy.DG_TO_TG: self.match_dg_to_tg_using_input,
            MatchingStrategy.NO_MATCH: MatchResult.NO_MATCH,
            MatchingStrategy.SUPER_MATCH: MatchResult.SUPER_MATCH,
        }

    def match_sort(self, reactant: str, product: str, input_lipids: Set[str], check_super_matched=False) -> Union[
        str, List[str], MatchResult]:
        try:
            reactant_family = self.lipid_parser.parse_lipid(reactant).headgroup
            product_family = self.lipid_parser.parse_lipid(product).headgroup
        except (ValueError, TypeError) as e:
            # --- ROBUSTNESS FIX: Catch parsing errors for malformed DB entries ---
            logger.warning(
                f"Could not parse reactant ('{reactant}') or product ('{product}') for matching. "
                f"Skipping this pair. Error: {e}"
            )
            # Return NO_MATCH to signal that this pair cannot be processed.
            return MatchResult.NO_MATCH
            # --------------------------------------------------------------------

        combined_families = f"{reactant_family}_{product_family}"
        logger.debug(f"Searching for strategies for the following families: {combined_families}")

        if check_super_matched:
            if combined_families in self.matching_strategies.get(MatchingStrategy.SUPER_MATCH, []):
                logger.debug(f"SUPER matching strategy found for {combined_families}")
                return MatchResult.SUPER_MATCH

        for strategy, pairs in self.matching_strategies.items():
            if combined_families in pairs:
                if strategy == MatchingStrategy.NO_MATCH:
                    logger.debug(f"NO_MATCH returned for {combined_families}, will grab other reactant to match with")
                    return MatchResult.NO_MATCH
                elif strategy == MatchingStrategy.SUPER_MATCH:
                    logger.debug(f"SUPER matching strategy found for {combined_families}")
                    return MatchResult.SUPER_MATCH
                else:
                    logger.debug(f"{strategy} found for {combined_families}")
                    return self.strategy_functions[strategy](reactant, product, input_lipids)

        logger.warning(f"No strategy was found for the following families: {combined_families}, product returned as is")
        return product

    # def has_direct_match(self, reactant_class: str, product_class: str) -> bool:
    #     """
    #     NEW HELPER: Checks if a direct, non-super, non-no-match strategy exists
    #     for a given reactant-product class pair.
    #     """
    #     combined_families = f"{reactant_class}_{product_class}"
    #     for strategy, pairs in self.matching_strategies.items():
    #         if strategy not in [MatchingStrategy.NO_MATCH, MatchingStrategy.SUPER_MATCH]:
    #             if combined_families in pairs:
    #                 return True
    #     return False

    def is_valid_pair(self, reactant_class: str, product_class: str) -> bool:
        """
        NEW, UNIFIED HELPER: Checks if a reactant-product class pair corresponds to ANY
        valid matching strategy (direct OR super), excluding NO_MATCH.
        """
        if not reactant_class or not product_class:
            return False

        combined_families = f"{reactant_class}_{product_class}"

        # Check all normal matching strategies
        for strategy, pairs in self.matching_strategies.items():
            if strategy != MatchingStrategy.NO_MATCH:
                if combined_families in pairs:
                    return True  # Found a valid direct or super match rule

        # This check might be redundant if SUPER_MATCH is in the main strategies, but it's safe.
        # Check dedicated super-matching strategies by constructing the key
        # This part is complex and might not be needed if your config is structured well.
        # Let's assume for now the main matching_strategies list is comprehensive.

        return False

    def super_match_sort(self, reactant1: str, reactant2: str, all_product_classes: Set[str],
                         input_lipids: Set[str]) -> List[Dict[str, str]]:
        """
        CORRECTED: Returns a LIST of product dictionaries to handle cases
        where multiple valid product combinations can be formed from one set of reactants.
        """
        reactant1_family = self.lipid_parser.parse_lipid(reactant1).headgroup
        reactant2_family = self.lipid_parser.parse_lipid(reactant2).headgroup
        sorted_reactants_family = sorted([reactant1_family, reactant2_family])

        # We only need to find the correct strategy once.
        for product_class in all_product_classes:
            combined_families_key = f"{sorted_reactants_family[0]}${sorted_reactants_family[1]}_{product_class}"

            for strategy, pairs in self.super_matching_strategies.items():
                if combined_families_key in pairs:
                    # Found the right strategy. Execute it and return its result.
                    sorted_reactants = sorted([reactant1, reactant2])
                    logger.debug(f"Executing SUPER strategy {strategy.name} for {combined_families_key}")

                    # The dedicated handler will check all possibilities.
                    return self.super_strategy_functions[strategy](
                        sorted_reactants[0], sorted_reactants[1], input_lipids
                    )

        logger.warning(
            f"No super-matching strategy found for {reactant1_family}, {reactant2_family} -> {all_product_classes}")
        return []

    def match_lipids_completely(self, reactant: str, product: str, input_lipids: Set[str]) -> str:
        # TODO: Initialize product lipid using its component class
        reactant_components = self.lipid_parser.parse_lipid(reactant)
        product_components = self.lipid_parser.parse_lipid(product)
        product_headgroup = product_components.headgroup

        # Handle species-level matching first
        if reactant_components.level in (LipidLevel.SPECIES, LipidLevel.MOLECULAR_SPECIES):
            return f"{product_headgroup}({reactant_components.species_composition})"

        if isinstance(reactant_components, SphingoLipidComponents):
            # product_components = self.lipid_parser.parse_lipid(product)
            # synthesized_product = SphingoLipidComponents(product_components.headgroup,
            #                                              reactant_components.backbone,
            #                                              reactant_components.n_acyl_chain,
            #                                              reactant_components.hydroxy,
            #                                              )
            hydroxy_str = f"({reactant_components.hydroxy})" if reactant_components.hydroxy else ""
            return f"{product_headgroup}({reactant_components.backbone}/{reactant_components.n_acyl_chain}{hydroxy_str})"
        elif isinstance(reactant_components, (GlycerophosphoLipidComponents, TriAcylComponents)):
            hydroxy_str = f"({reactant_components.hydroxy})" if reactant_components.hydroxy else ""
            return f"{product_headgroup}({reactant_components.sn1_chain}/{reactant_components.sn2_chain}{hydroxy_str})"
        elif isinstance(reactant_components, LysoGlycerophospholipidComponents) and isinstance(product_components, MonoAcylComponents):
            return f"{product_headgroup}({reactant_components.active_chain})"
        elif isinstance(reactant_components, LysoGlycerophospholipidComponents):
            active_chain = reactant_components.active_chain
            active_chain_position = reactant_components.active_chain_position
            if '_' in product_headgroup:  # if lyso lipid has _sn1 or _sn2 info, get rid of them
                product_headgroup = product_headgroup.split('_')[0]
            elif '(' in product_headgroup: # lyso lipid has (O-) or (P-)
                product_headgroup = product_headgroup.split('(')[0]
            if active_chain_position == ChainPosition.SN1:
                return f"{product_headgroup}({active_chain}/0:0)"
            else:
                return f"{product_headgroup}(0:0/{active_chain})"
            # return f"{product_headgroup}({reactant_components.sn1_chain}/{reactant_components.sn2_chain})"
        else:
            raise ValueError(f"Unsupported lipid type for complete matching: {type(reactant_components)}")

    def match_pc_to_lpc_lipids(self, reactant: str, product: str, input_lipids: Set[str]) -> str:
        reactant_components = self.lipid_parser.parse_lipid(reactant)
        product_components = self.lipid_parser.parse_lipid(product)

        if not isinstance(reactant_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(reactant_components)}")
        # if not isinstance(product_components, LysoGlycerophospholipidComponents):
        #     raise ValueError(f"Expected LysoGlycerophospholipidComponents, got {type(product_components)})")

        if product_components.headgroup.endswith('_sn1'):
            headgroup_without_sn_location = product_components.headgroup.split('_')[0]
            reactant_sn1_chain = reactant_components.sn1_chain
            matched_product = f"{headgroup_without_sn_location}({reactant_sn1_chain}/0:0)"
        elif product_components.headgroup.endswith('_sn2'):
            headgroup_without_sn_location = product_components.headgroup.split('_')[0]
            reactant_sn2_chain = reactant_components.sn2_chain
            matched_product = f"{headgroup_without_sn_location}(0:0/{reactant_sn2_chain})"
        elif product_components.headgroup.endswith('(O-)') or product_components.headgroup.endswith('(P-)'):
            headgroup_without_O_P = product_components.headgroup.split('(')[0]
            reactant_sn1_chain = reactant_components.sn1_chain
            matched_product = f"{headgroup_without_O_P}({reactant_sn1_chain}/0:0)"
        else:
            raise ValueError(f"Unexpected lyso type as a product: {product}")

        return matched_product

    def match_unmarked_lpc_to_pc_using_input_lipids(self, reactant: str, product: str, input_lipids: Set[str]) -> str:
        raise NotImplementedError("match_unmarked_lpc_to_pc_using_input_lipids not implemented yet")

    def match_lpc_to_pc_using_input_lipids(self, reactant: str, product: str, input_lipids: Set[str]) -> Set[str]:
        start_time = time.time()
        reactant_components = self.lipid_parser.parse_lipid(reactant)
        product_components = self.lipid_parser.parse_lipid(product)

        if not isinstance(reactant_components, LysoGlycerophospholipidComponents):
            raise ValueError(f"Expected LysoGlycerophospholipidComponents, got {type(reactant_components)}")

        # --- NEW LOGIC for single search mode ---
        if not input_lipids:
            product_headgroup = product_components.headgroup
            active_chain = reactant_components.active_chain

            # Synthesize a template with a placeholder for the unknown chain.
            if reactant_components.active_chain_position == ChainPosition.SN1:
                synthesized_template = f"{product_headgroup}({active_chain}/X:Y)"
            else:  # SN2
                synthesized_template = f"{product_headgroup}(X:Y/{active_chain})"

            logger.debug(f"Single search: Synthesized template '{synthesized_template}' for {reactant}")
            return {synthesized_template}

        reactant_active_chain = reactant_components.active_chain
        reactant_active_chain_position = reactant_components.active_chain_position
        product_headgroup = product_components.headgroup
        matched_products = set()
        input_pcs = self.headgroup_map.get(product_headgroup, set())
        for possible_gp_lipid in input_pcs:
            gp_lipid_components = self.lipid_parser.parse_lipid(possible_gp_lipid)

            if not isinstance(gp_lipid_components, GlycerophosphoLipidComponents):
                raise ValueError(f"Expected GlycerophospholipidComponents, got {type(gp_lipid_components)}")

            # If reactant has active chain at sn1, then it needs to match sn1 chain in the product
            # If reactant has active chain at sn2, then it needs to match sn2 chain in the product
            if reactant_active_chain_position == ChainPosition.SN1:
                if gp_lipid_components.sn1_chain == reactant_active_chain:
                    matched_products.add(possible_gp_lipid)
            else:  # sm2 position
                if gp_lipid_components.sn2_chain == reactant_active_chain:
                    matched_products.add(possible_gp_lipid)

        # if reactant_components.headgroup.endswith('_sn2'): # Look for any PC with matching sn2 chain
        #     reactant_active_chain = reactant_components.active_chain
        #     reactant_active_chain_position = reactant_components.active_chain_position
        #     product_headgroup = product_components.headgroup
        #     matched_products = set()
        #     input_pcs = self.headgroup_map.get(product_headgroup, set())
        #     for possible_gp_lipid in input_pcs:
        #         gp_lipid_components = self.lipid_parser.parse_lipid(possible_gp_lipid)
        #         if not isinstance(gp_lipid_components, GlycerophosphoLipidComponents):
        #             raise(f"Expected GlycerophospholipidComponents, got {type(gp_lipid_components)}")
        #         if gp_lipid_components.sn2_chain == reactant_active_chain:
        #             matched_products.add(possible_gp_lipid)
        #         # if (lipid.startswith(starting_string_of_product) and lipid.endswith(ending_string_of_product) and not
        #         #         (lipid.startswith(f"{product_headgroup}(O-") or lipid.startswith(f"{product_headgroup}(P-"))):
        #         #     matched_products.add(lipid)
        #
        # else:
        #     reactant_sn1_chain = reactant_components.sn1_chain
        #     product_headgroup = product_components.headgroup
        #     if '(' in product_headgroup:
        #         product_headgroup = product_headgroup.split('(')[0]
        #     matched_products = set()
        #     input_pcs = self.headgroup_map.get(product_headgroup, set())
        #     for possible_gp_lipid in input_pcs:
        #         gp_lipid_components = self.lipid_parser.parse_lipid(possible_gp_lipid)
        #         if not isinstance(gp_lipid_components, GlycerophosphoLipidComponents):
        #             raise(f"Expected GlycerophospholipidComponents, got {type(gp_lipid_components)}")
        #         if gp_lipid_components.sn1_chain == reactant_sn1_chain:
        #             matched_products.add(possible_gp_lipid)
        #     starting_string_of_product = f"{product_headgroup}({reactant_sn1_chain}/"
        #     for lipid in input_lipids:
        #         if lipid.startswith(starting_string_of_product):
        #             if '/O-' in lipid:
        #                 logging.warning(f"Skipping di-ether")
        #                 continue
        #             matched_products.add(lipid)
        end_time = time.time()
        logger.info(f"match_lpc_to_pc_using_input_lipids execution time: {end_time - start_time} seconds")
        return matched_products

    def match_n_acyl_chain_using_input_lipids(self, reactant: str, product: str, input_lipids: Set[str]) -> Set[str]:
        start_time = time.time()
        # For Sph -> Cer reactions
        reactant_components = self.lipid_parser.parse_lipid(reactant)
        product_components = self.lipid_parser.parse_lipid(product)

        if not isinstance(reactant_components, SphingoidBaseComponents):
            raise ValueError(f"Expected SphingoidBaseComponents, got {type(reactant_components)}")

        # --- NEW LOGIC for single search mode ---
        if not input_lipids:
            # We don't have an N-acyl chain to add, so we create a backbone-level template.
            # Example: reactant Sph(d18:1) -> product template "Cer(d18:1)"
            product_headgroup = product_components.headgroup
            reactant_backbone = reactant_components.backbone
            synthesized_template = f"{product_headgroup}({reactant_backbone})"
            logger.debug(f"Single search: Synthesized template '{synthesized_template}' for {reactant}")
            return {synthesized_template}

        ceramide_headgroup = product_components.headgroup
        sphingoid_base_backbone = reactant_components.backbone
        matched_products = set()
        input_pcs = self.headgroup_map.get(ceramide_headgroup, set())
        for possible_cer_lipid in input_pcs:
            possible_cer_components = self.lipid_parser.parse_lipid(possible_cer_lipid)
            if not isinstance(possible_cer_components, SphingoLipidComponents):
                raise ValueError(f"Expected SphingoLipidComponents, got {type(possible_cer_components)}")
            if possible_cer_components.backbone == sphingoid_base_backbone:
                matched_products.add(possible_cer_lipid)

        end_time = time.time()
        logger.info(f"match_n_acyl_chain_using_input_lipids execution time: {end_time - start_time} seconds")
        return matched_products

    def match_dg_to_tg_using_input(self, reactant: str, product: str, input_lipids: Set[str]) -> Set[str]:
        reactant_components = self.lipid_parser.parse_lipid(reactant)
        product_components = self.lipid_parser.parse_lipid(product)

        if not isinstance(reactant_components, TriAcylComponents):
            raise ValueError(f"Expected TriAcylComponents, got {type(reactant_components)}")

        # --- NEW LOGIC for single search mode ---
        if not input_lipids:
            product_headgroup = product_components.headgroup
            dg_sn1_chain = reactant_components.sn1_chain
            dg_sn2_chain = reactant_components.sn2_chain

            # Synthesize a template with a placeholder for the unknown sn3 chain.
            synthesized_template = f"{product_headgroup}({dg_sn1_chain}/{dg_sn2_chain}/X:Y)"
            logger.debug(f"Single search: Synthesized template '{synthesized_template}' for {reactant}")
            return {synthesized_template}

        tg_headgroup = product_components.headgroup
        dg_sn1_chain = reactant_components.sn1_chain
        dg_sn2_chain = reactant_components.sn2_chain
        matched_products = set()
        input_tgs = self.headgroup_map.get(tg_headgroup, set())
        for possible_tg_lipid in input_tgs:
            tg_components = self.lipid_parser.parse_lipid(possible_tg_lipid)
            if not isinstance(tg_components, TriAcylComponents):
                raise ValueError(f"Expected TriAcylComponents, got {type(tg_components)}")
            if dg_sn1_chain != "" and tg_components.sn1_chain == dg_sn1_chain and tg_components.sn2_chain == dg_sn2_chain:
                matched_products.add(possible_tg_lipid)

        return matched_products


    def match_backbone(self, reactant: str, product: str, input_lipids: Set[str]) -> str:
        # For Cer -> Sph
        # or HexSph -> Sph
        reactant_components = self.lipid_parser.parse_lipid(reactant)
        product_components = self.lipid_parser.parse_lipid(product)

        if not (isinstance(reactant_components, SphingoLipidComponents) or
                isinstance(reactant_components, SphingoidBaseComponents)):
            raise ValueError(f"Expected SphingoLipidComponents, got {type(reactant_components)}")


        reactant_backbone = reactant_components.backbone
        matched_product = f"{product_components.headgroup}({reactant_backbone})"
        return matched_product

    def match_n_acyl_chain_for_desaturase(self, reactant: str, product: str, input_lipids: Set[str]) -> Set[str]:
        """
        Synthesizes a desaturated product from a saturated backbone.

        MODIFIED: Now returns a Set. If the reactant backbone ends with :0, it becomes :1.
        If it ends with :1, it becomes :0. Otherwise returns MatchResult.NO_MATCH.
        """
        reactant_components = self.lipid_parser.parse_lipid(reactant)
        product_components = self.lipid_parser.parse_lipid(product)

        if not isinstance(reactant_components, SphingoLipidComponents):
            raise ValueError(f"Expected SphingoLipidComponents, got {type(reactant_components)}")

        reactant_backbone = reactant_components.backbone

        # --- BACKBONE TRANSFORMATION LOGIC ---
        if reactant_backbone.endswith(":0"):
            # Change :0 to :1
            backbone_base = reactant_backbone.split(':')[0]
            matched_product_backbone = f"{backbone_base}:1"
        elif reactant_backbone.endswith(":1"):
            # Change :1 to :0
            backbone_base = reactant_backbone.split(':')[0]
            matched_product_backbone = f"{backbone_base}:0"
        else:
            # Neither :0 nor :1 - return no match
            logger.debug(
                f"Desaturase match skipped: Reactant backbone '{reactant_backbone}' does not end with :0 or :1.")
            return MatchResult.NO_MATCH
        # ----------------------------------------

        reactant_n_acyl = reactant_components.n_acyl_chain
        matched_product = f"{product_components.headgroup}({matched_product_backbone}/{reactant_n_acyl})"

        logger.debug(f"Desaturase match successful: Synthesized '{matched_product}' from '{reactant}'.")
        # Return the valid new product inside a set.
        return {matched_product}

    def match_cer_gpl_to_acer(self, cer_reactant: str, gpl_reactant: str, acer_product: str,
                              input_lipids: Set[str]) -> Union[str, List[str]]:
        # Implementation for Cer + GPL -> AcylCer matching
        cer_reactant_components = self.lipid_parser.parse_lipid(cer_reactant)
        gpl_reactant_components = self.lipid_parser.parse_lipid(gpl_reactant)

        if not isinstance(cer_reactant_components, SphingoLipidComponents):
            raise ValueError(f"Expected SphingoLipidComponents, got {type(cer_reactant_components)}")
        if not isinstance(gpl_reactant_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(gpl_reactant_components)}")

        cer_backbone = cer_reactant_components.backbone
        cer_n_acyl = cer_reactant_components.n_acyl_chain
        gpl_sn1_chain = gpl_reactant_components.sn1_chain
        gpl_sn2_chain = gpl_reactant_components.sn2_chain

        sn1_acyl_cer_product = f"1-O-{gpl_sn1_chain}-Cer({cer_backbone}/{cer_n_acyl})"
        sn2_acyl_cer_product = f"1-O-{gpl_sn2_chain}-Cer({cer_backbone}/{cer_n_acyl})"

        return [sn1_acyl_cer_product, sn2_acyl_cer_product]

    def match_lysogpl_gpl_to_gpl(self, lysogpl_reactant: str, gpl_reactant: str, gpl_product: str,
                                 input_lipids: Set[str]) -> str:
        # Implementation for LPC_sn1 + PE -> PC
        lysogpl_reactant_components = self.lipid_parser.parse_lipid(lysogpl_reactant)
        gpl_reactant_components = self.lipid_parser.parse_lipid(gpl_reactant)
        gpl_product_components = self.lipid_parser.parse_lipid(gpl_product)

        if not isinstance(lysogpl_reactant_components, LysoGlycerophospholipidComponents):
            raise ValueError(f"Expected LysoGlycerophospholipidComponents, got {type(lysogpl_reactant_components)}")
        if not isinstance(gpl_reactant_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(gpl_reactant_components)}")

        lysogpl_sn1_chain = lysogpl_reactant_components.active_chain
        gpl_sn2_chain = gpl_reactant_components.sn2_chain

        matched_product = f"{gpl_product_components.headgroup}({lysogpl_sn1_chain}/{gpl_sn2_chain})"
        return matched_product

    def match_nape_to_pa(self, nape_reactant: str, pa_product: str, input_lipids: Set[str]) -> str:
        # implementation for NAPE -> PA, NAPE -> LPA_sn1
        nape_reactant_components = self.lipid_parser.parse_lipid(nape_reactant)
        pa_product_components = self.lipid_parser.parse_lipid(pa_product)

        if not isinstance(nape_reactant_components, TriAcylComponents):
            raise ValueError(f"Expected TriAcylComponents, got {type(nape_reactant_components)}")

        nape_s1_chain = nape_reactant_components.sn1_chain
        nape_s2_chain = nape_reactant_components.sn2_chain

        matched_product = f"{pa_product_components.headgroup}({nape_s1_chain}/{nape_s2_chain})"
        return matched_product

    def match_nape_gpnae(self, nape_reactant: str, nae_product: str, input_lipids: Set[str]) -> str:
        # implementation for NAPE -> GP-NAE
        nape_reactant_components = self.lipid_parser.parse_lipid(nape_reactant)
        gpnae_product_components = self.lipid_parser.parse_lipid(nae_product)

        if not isinstance(nape_reactant_components, TriAcylComponents):
            raise ValueError(f"Expected TriAcylComponents, got {type(nape_reactant_components)}")

        nape_acyl_chain = nape_reactant_components.sn3_chain

        matched_product = f"{gpnae_product_components.headgroup}(0:0/0:0/{nape_acyl_chain})"
        return matched_product

    def match_nape_nalpe(self, nape_reactant: str, nalpe_product: str, input_lipids: Set[str]) -> str:
        # implementation for NAPE -> NALPE (N-acyl-1-acyl-sn-glycero-3-phosphoethanolamine)
        nape_reactant_components = self.lipid_parser.parse_lipid(nape_reactant)
        nalpe_product_components = self.lipid_parser.parse_lipid(nalpe_product)

        if not isinstance(nape_reactant_components, TriAcylComponents):
            raise ValueError(f"Expected TriAcylComponents, got {type(nape_reactant_components)}")

        nape_sn1_chain = nape_reactant_components.sn1_chain
        nape_n_acyl_chain = nape_reactant_components.sn3_chain

        matched_product = f"{nalpe_product_components.headgroup}({nape_sn1_chain}/0:0/{nape_n_acyl_chain})"
        return matched_product

    def match_nape_nae(self, nape_reactant, nae_product: str, input_lipids: Set[str]) -> str:
        # implementation for NAPE -> NAE
        nape_reactant_components = self.lipid_parser.parse_lipid(nape_reactant)
        nae_product_components = self.lipid_parser.parse_lipid(nae_product)

        if not isinstance(nape_reactant_components, TriAcylComponents):
            raise ValueError(f"Expected TriAcylComponents, got {type(nape_reactant_components)}")

        nape_acyl_chain = nape_reactant_components.sn3_chain

        matched_product = f"{nae_product_components.headgroup}({nape_acyl_chain})"
        return matched_product

    def match_pc_pe_to_nape(self, pc_reactant: str, pe_reactant: str, nape_product: str, input_lipids: Set[str]) -> str:
        # implementation for PC + PE -> NAPE
        pc_reactant_components = self.lipid_parser.parse_lipid(pc_reactant)
        pe_reactant_components = self.lipid_parser.parse_lipid(pe_reactant)
        nape_product_components = self.lipid_parser.parse_lipid(nape_product)

        if not isinstance(pc_reactant_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(pc_reactant_components)}")
        if not isinstance(pe_reactant_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(pe_reactant_components)}")

        pe_sn1_chain = pe_reactant_components.sn1_chain
        pe_sn2_chain = pe_reactant_components.sn2_chain
        pc_sn2_chain = pc_reactant_components.sn2_chain
        nape_headgroup = nape_product_components.headgroup

        matched_product = f"{nape_headgroup}({pe_sn1_chain}/{pe_sn2_chain}/{pc_sn2_chain})"
        return matched_product

    def match_cl_to_mlcl(self, cl_reactant: str, mlcl_product: str, input_lipids: Set[str]) -> str:
        # implement for CL -> MLCL
        cl_reactant_components = self.lipid_parser.parse_lipid(cl_reactant)
        mlcl_product_components = self.lipid_parser.parse_lipid(mlcl_product)

        if not isinstance(cl_reactant_components, CardiolipinComponents):
            raise ValueError(f"Expected CardiopinComponents, got {type(cl_reactant_components)}")

        cl_sn1_chain = cl_reactant_components.sn1_chain
        cl_sn2_chain = cl_reactant_components.sn2_chain
        cl_sn1_prime_chain = cl_reactant_components.sn1_prime_chain
        mlcl_headgroup = mlcl_product_components.headgroup

        matched_product = f"{mlcl_headgroup}({cl_sn1_chain}/{cl_sn2_chain}/{cl_sn1_prime_chain}/0:0)"
        return matched_product


    def match_pc_pe_to_pc(self, pc_reactant: str, pe_reactant: str, pc_product: str, input_lipids: Set[str]) -> str:
        pc_reactant_components = self.lipid_parser.parse_lipid(pc_reactant)
        pe_reactant_components = self.lipid_parser.parse_lipid(pe_reactant)

        if not isinstance(pc_reactant_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(pc_reactant_components)}")
        if not isinstance(pe_reactant_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(pe_reactant_components)}")


        pe_sn1_chain = pe_reactant_components.sn1_chain
        pe_sn2_chain = pe_reactant_components.sn2_chain
        pc_sn1_chain = pc_reactant_components.sn1_chain
        pc_sn2_chain = pc_reactant_components.sn2_chain

        if any(var == '0:0' for var in [pe_sn1_chain, pe_sn2_chain, pc_sn1_chain, pc_sn2_chain]):
            return 'PC'

        matched_pc_product = f"{pc_reactant_components.headgroup}({pe_sn1_chain}/{pc_sn2_chain})"
        return matched_pc_product

    def match_pc_pe_to_pe(self, pc_reactant: str, pe_reactant: str, pe_product: str, input_lipids: Set[str]) -> str:
        pc_reactant_components = self.lipid_parser.parse_lipid(pc_reactant)
        pe_reactant_components = self.lipid_parser.parse_lipid(pe_reactant)

        if not isinstance(pc_reactant_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(pc_reactant_components)}")
        if not isinstance(pe_reactant_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(pe_reactant_components)}")

        pe_sn1_chain = pe_reactant_components.sn1_chain
        pe_sn2_chain = pe_reactant_components.sn2_chain
        pc_sn1_chain = pc_reactant_components.sn1_chain
        pc_sn2_chain = pc_reactant_components.sn2_chain

        if any(var == '0:0' for var in [pe_sn1_chain, pe_sn2_chain, pc_sn1_chain, pc_sn2_chain]):
            return 'PE'

        matched_pe_product = f"{pe_reactant_components.headgroup}({pc_sn1_chain}/{pe_sn2_chain})"
        return matched_pe_product

    # def match_pco_spb_to_cer(self, pco_reactant: str, spb_reactant: str, cer_product: str, input_lipids: Set[str]) -> str:
    #     pco_reactant_components = self.lipid_parser.parse_lipid(pco_reactant)
    #     spb_reactant_components = self.lipid_parser.parse_lipid(spb_reactant)
    #     cer_product_components = self.lipid_parser.parse_lipid(cer_product)
    #     cer_product_headgroup = cer_product_components.headgroup
    #
    #     if not isinstance(pco_reactant_components, GlycerophosphoLipidComponents):
    #         raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(pco_reactant_components)}")
    #     if not isinstance(spb_reactant_components, SphingoidBaseComponents):
    #         raise ValueError(f"Expected SphingoidBaseComponents, got {type(spb_reactant_components)}")
    #
    #     pco_reactant_sn2 = pco_reactant_components.sn2_chain
    #     spb_reactant_backbone = spb_reactant_components.backbone
    #
    #     matched_cer_product = f"{cer_product_headgroup}({spb_reactant_backbone}/{pco_reactant_sn2})"
    #     return matched_cer_product

    def match_cer_pc_to_acer_lpc(self, cer_reactant: str, pc_reactant: str,
                                         input_lipids: Set[str]) -> List[Dict[str, str]]:
        """
        CORRECTED SUPER-STRATEGY: Checks ALL possible outcomes for Cer + PC => ACer + LPC
        and returns a list of all valid product set dictionaries.
        """
        potential_acers = self.match_cer_gpl_to_acer(cer_reactant, pc_reactant, "placeholder", input_lipids)
        potential_lpc_sn1 = self.match_pc_to_lpc_lipids(pc_reactant, "LPC_sn1", input_lipids)
        potential_lpc_sn2 = self.match_pc_to_lpc_lipids(pc_reactant, "LPC_sn2", input_lipids)

        valid_product_combinations = []

        # Check Scenario 1 (sn1 transfer from PC)
        acer_from_sn1 = potential_acers[0]
        lpc_from_sn2 = potential_lpc_sn2
        if acer_from_sn1 in input_lipids and lpc_from_sn2 in input_lipids:
            valid_product_combinations.append({"1-O-acyl-Cer": acer_from_sn1, "LPC_sn2": lpc_from_sn2})

        # Check Scenario 2 (sn2 transfer from PC)
        acer_from_sn2 = potential_acers[1]
        lpc_from_sn1 = potential_lpc_sn1
        if acer_from_sn2 in input_lipids and lpc_from_sn1 in input_lipids:
            valid_product_combinations.append({"1-O-acyl-Cer": acer_from_sn2, "LPC_sn1": lpc_from_sn1})

        return valid_product_combinations

    def match_pco_sph_to_lpco_cer(self, pco_reactant: str, sph_reactant: str,
                                  input_lipids: Set[str]) -> List[Dict[str, str]]:
        """
        NEW DEDICATED SUPER-STRATEGY: Handles PC(O-) + Sph => LPC(O-) + Cer.
        Synthesizes both products and returns them in a dictionary if valid.
        """
        try:
            pco_components = self.lipid_parser.parse_lipid(pco_reactant)
            sph_components = self.lipid_parser.parse_lipid(sph_reactant)
        except ValueError as e:
            logger.error(f"Could not parse reactants for PCO_SPH_TO_LPCO_CER strategy: {e}")
            return []

        if not isinstance(pco_components, GlycerophosphoLipidComponents):
            raise ValueError(f"Expected GlycerophosphoLipidComponents, got {type(pco_components)}")


        # Logic: The sn2 chain from PC(O-) is transferred to the Sphingoid base to make Ceramide.
        # The remaining PC(O-) becomes an LPC(O-).

        # 1. Synthesize the Ceramide product
        # The headgroup is 'Cer'. The backbone is from Sph. The N-acyl chain is the sn2 from PC(O-).
        synthesized_cer = f"Cer({sph_components.backbone}/{pco_components.sn2_chain})"

        # 2. Synthesize the LPC(O-) product
        # It retains the sn1 chain of the original PC(O-).
        synthesized_lpc = f"LPC({pco_components.sn1_chain}/0:0)"

        # 3. Validate that BOTH products exist in the user's input list.
        if synthesized_cer in input_lipids and synthesized_lpc in input_lipids:
            logger.debug(f"Found valid complex match: {synthesized_cer} and {synthesized_lpc}")
            # Return the result as a list containing one dictionary, matching the expected return type.
            return [{
                "Cer": synthesized_cer,
                "LPC(O-)": synthesized_lpc
            }]

        return []
