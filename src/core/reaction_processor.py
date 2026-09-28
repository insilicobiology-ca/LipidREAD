import time
import logging
from typing import List, Set, Union, Dict, Optional
from src.models.models import Reaction, ReactionType, MatchResult
from src.models.lipid_level import LipidLevel
from src.matching.lipid_matcher import LipidMatcher
from src.parsing.lipid_parser import LipidParser
from src.data.data_access import LipidDataAccess
from src.models.lipid_components import LipidComponents, SphingoLipidComponents

logger = logging.getLogger(__name__)


class ReactionProcessor:
    """
    A class responsible for processing individual reactions.

    This class handles the logic for matching reaction products, processing
    super reactions, and managing various reaction-specific operations.
    """

    def __init__(self, lipid_matcher: LipidMatcher, lipid_parser: LipidParser, data_access: LipidDataAccess):
        """
        Initialize the ReactionProcessor.

        Args:
            lipid_matcher (LipidMatcher): An instance of LipidMatcher for lipid matching operations.
            lipid_parser (LipidParser): An instance of LipidParser for parsing lipid structures.
            data_access (LipidDataAccess): An instance of LipidDataAccess for database operations.
        """
        self.lipid_matcher = lipid_matcher
        self.lipid_parser = lipid_parser
        self.data_access = data_access
        self.headgroup_map = None

    def set_headgroup_map(self, headgroup_map: Dict[str, Set[str]]) -> None:
        """Set the headgroup map for looking up lipids by headgroup."""
        self.headgroup_map = headgroup_map

    # def set_preprocessed_lipids(self, preprocessed_lipids: Dict[str, Dict[str, Dict[str, Set[str]]]]):
    #     self.preprocessed_lipids = preprocessed_lipids

    # def process_reaction(self, reaction: Reaction, original_lipid: str, molecule_id: int,
    #                      input_lipids: Set[str]) -> List[Reaction]:
    #     """
    #     Process a single reaction and return a list of matched reactions.
    #
    #     This method handles the matching of reaction products and creates new reactions
    #     based on the matching results. It also handles super reactions differently.
    #
    #     Args:
    #         reaction (Reaction): The reaction to process.
    #         original_lipid (str): The original lipid involved in the reaction.
    #         molecule_id (int): The ID of the molecule associated with the reaction.
    #         input_lipids (Set[str]): A set of all input lipids.
    #
    #     Returns:
    #         List[Reaction]: A list of processed reactions with matched products.
    #     """
    #     logger.debug(f"Processing reaction_enzyme_id: {reaction.reaction_enzyme_id}")
    #     try:
    #         product_components = self.lipid_parser.parse_lipid(reaction.product)
    #
    #         if reaction.reaction_type == ReactionType.SUPER:
    #             logger.debug(f"Processing reaction as SUPER reaction")
    #             return self._process_super_reaction(reaction, original_lipid, molecule_id, input_lipids,
    #                                                 product_components)
    #         else:
    #             logger.debug(f"Processing reaction as REGULAR reaction")
    #             return self._process_regular_reaction(reaction, original_lipid, input_lipids, product_components)
    #     except ValueError as e:
    #         logger.warning(f"Error parsing lipid in reaction {reaction.reaction_id}: {str(e)}")
    #         return []
    #     except Exception as e:
    #         logger.error(f"Unexpected error in process_reaction for reaction {reaction.reaction_id}: {str(e)}")
    #         return []

    def process_reaction_group(self, reaction_pairs: List[Reaction], primary_reactant: str,
                               molecule_id: int, input_lipids: Optional[Set[str]],
                               validate_products: bool = True) -> List[Reaction]:
        """
        Orchestrates the processing of a single biochemical event.
        MODIFIED: Now logs a warning for any None participants before filtering them.
        """
        if not reaction_pairs:
            return []

        base_info = reaction_pairs[0]

        # --- KEY CHANGE: Check for and log None values before filtering ---
        all_db_participants_raw = [p for pair in reaction_pairs for p in pair.products] + \
                                  [r for pair in reaction_pairs for r in pair.reactants]

        if None in all_db_participants_raw:
            # Provide a detailed warning with context to help with debugging.
            logger.warning(
                f"Data Quality Issue: Found a NULL participant in reaction event. "
                f"Context - Reaction ID: {base_info.reaction_id}, "
                f"Enzyme: '{base_info.enzyme_name}', "
                f"Primary Lipid: '{primary_reactant}' (Molecule ID: {molecule_id}). "
                f"This likely indicates a molecule in the 'reaction_pairs' table is missing a 'cleaned_abbreviation'. "
                f"The NULL participant will be ignored for this reaction."
            )
            # Filter out the None values to prevent a crash.
            all_db_participants = [p for p in all_db_participants_raw if p is not None]
        else:
            all_db_participants = all_db_participants_raw
        # --------------------------------------------------------------------

        # If after filtering, there are no participants left, we can't proceed.
        if not all_db_participants:
            return []

        is_pre_specified = all(self._is_fully_specified(participant) for participant in all_db_participants)

        if is_pre_specified:
            return self._handle_pre_specified_molecular_reaction(base_info, input_lipids, validate_products)

        if base_info.reaction_type == ReactionType.REGULAR:
            return self._process_regular_reaction(base_info, primary_reactant, input_lipids, validate_products)

        elif base_info.reaction_type == ReactionType.SUPER:
            return self._process_super_reaction(base_info, primary_reactant, molecule_id, reaction_pairs, input_lipids,
                                                validate_products)

        return []

    def _process_regular_reaction(self, base_info: Reaction, primary_reactant: str,
                                  input_lipids: Optional[Set[str]], validate_products: bool = True) -> List[Reaction]:
        """Handles simple A -> B reactions, with conditional product validation."""

        # --- ROBUSTNESS FIX: Check for empty or None product ---
        if not base_info.products or base_info.products[0] is None:
            logger.warning(
                f"Skipping regular reaction for reactant '{primary_reactant}' due to missing or NULL product. "
                f"Context - Reaction ID: {base_info.reaction_id}, Enzyme: '{base_info.enzyme_name}'"
            )
            return []
        # --------------------------------------------------------


        product_template = base_info.products[0]

        synthesized_products = []

        # --- KEY CHANGE: Use the new, more precise check ---
        if self._is_fully_specified(product_template):
            # If the product from the DB is already a complete molecule, treat it as is.
            # No matching is needed.
            synthesized_products = [product_template]
        else:
            # If the product is a template (class, backbone, etc.),
            # we must perform matching to synthesize a specific product.
            synthesized_products = self._handle_headgroup_only_product(
                primary_reactant,
                product_template,
                input_lipids or set()
            )

        valid_reactions = []
        for p in synthesized_products:
            # The conditional validation logic remains the same and is correct.
            if not validate_products or (input_lipids and p in input_lipids):
                valid_reactions.append(base_info._replace(reactants=[primary_reactant], products=[p]))

        return valid_reactions

    def _process_super_reaction(self, base_info: Reaction, primary_reactant: str, molecule_id: int,
                                reaction_pairs: List[Reaction], input_lipids: Optional[Set[str]],
                                validate_products: bool = True) -> List[Reaction]:
        """
        Handles SUPER reactions with conditional product validation.
        """

        # --- ROBUSTNESS FIX: Safely parse participants and log errors --
        all_reactant_classes = set()
        all_product_classes = set()
        try:
            for pair in reaction_pairs:
                for r in pair.reactants:
                    if r: all_reactant_classes.add(self.lipid_parser.parse_lipid(r).headgroup)
                for p in pair.products:
                    if p: all_product_classes.add(self.lipid_parser.parse_lipid(p).headgroup)
        except ValueError as e:
            logger.warning(
                f"Skipping super reaction due to unparsable participant. "
                f"Context - Reaction ID: {base_info.reaction_id}, Enzyme: '{base_info.enzyme_name}'. "
                f"Error: {e}"
            )
            return []  # Skip this entire reaction group if a participant is unknown
        # ---------------------------------------------------------------

        final_reactions_for_output = []

        for product_class in all_product_classes:
            # Pass the correct input_lipids (or an empty set for unconstrained mode) to the matcher
            match_result = self.lipid_matcher.match_sort(primary_reactant, product_class, input_lipids or set(),
                                                         check_super_matched=True)

            processed_pairs = []

            # --- ROUTING LOGIC ---
            if isinstance(match_result, (str, list)):
                # Direct match found. Pass the flag to the helper.
                processed_pairs = self._handle_super_reaction_matched_product(
                    base_info, primary_reactant, match_result, molecule_id, input_lipids, validate_products
                )

            elif match_result == MatchResult.SUPER_MATCH:
                # SUPER_MATCH logic. Pass the flag to the helper.
                complex_reactions = self._handle_super_reaction_super_match(
                    base_info, primary_reactant, all_product_classes, molecule_id, input_lipids, validate_products
                )
                final_reactions_for_output.extend(complex_reactions)
                # For SUPER_MATCH, we assume all products are handled at once, so we can break.
                break
            else:
                # NO_MATCH, continue to the next product class.
                continue

            # --- RE-ASSEMBLY FOR DIRECT MATCHES ---
            # This block only runs if a direct match was found.
            for pair in processed_pairs:
                # `pair.reactants` is e.g., {'PC(16:0/18:0)', 'Cer'}
                # `pair.products` is e.g., {'DG(16:0/18:0)'}

                specific_reactants = set(pair.reactants)
                specific_products = set(pair.products)

                used_reactant_classes = {self.lipid_parser.parse_lipid(r).headgroup for r in specific_reactants}
                used_product_classes = {self.lipid_parser.parse_lipid(p).headgroup for p in specific_products}

                # Find the remaining blackboxed participants
                secondary_reactants = all_reactant_classes - used_reactant_classes
                secondary_products = all_product_classes - used_product_classes

                # Combine to form the final multi-participant reaction
                final_reactants = sorted(list(specific_reactants.union(secondary_reactants)))
                final_products = sorted(list(specific_products.union(secondary_products)))

                final_full_reaction = base_info._replace(reactants=final_reactants, products=final_products)
                final_reactions_for_output.append(final_full_reaction)

        return final_reactions_for_output

    # def _process_super_reaction_pair(self, base_info: Reaction, primary_reactant: str, product_class: str,
    #                                  molecule_id: int, input_lipids: Set[str]) -> List[Reaction]:
    #     """
    #     Helper that contains the logic from your old `_process_super_reaction`.
    #     It determines how to handle a single reactant-product pair within a super reaction.
    #     """
    #     match_result = self.lipid_matcher.match_sort(primary_reactant, product_class, input_lipids,
    #                                                  check_super_matched=True)
    #
    #     # We need the reaction_enzyme_id for the handlers.
    #     reaction_enzyme_id = base_info.reaction_enzyme_id
    #
    #     if isinstance(match_result, (str, list)):
    #         return self._handle_super_reaction_matched_product(base_info, primary_reactant, match_result, molecule_id,
    #                                                            input_lipids)
    #     elif match_result == MatchResult.SUPER_MATCH:
    #         return self._handle_super_reaction_super_match(base_info, primary_reactant, product_class, molecule_id,
    #                                                        input_lipids)
    #     elif match_result == MatchResult.NO_MATCH:
    #         return self._handle_super_reaction_no_match(base_info, primary_reactant, product_class, molecule_id,
    #                                                     input_lipids)
    #
    #     return []

    # def _handle_super_reaction_full_product(self, base_info: Reaction, input_lipids: Set[str]) -> List[Reaction]:
    #     """
    #     NEW HELPER: Handles reactions that are already defined with specific molecular species.
    #     It validates that all participants are present in the input list before creating the reaction.
    #     """
    #     all_participants = self.data_access.get_all_participants(base_info.reaction_enzyme_id)
    #     if not all_participants:
    #         return []
    #
    #     all_reactants, all_products = all_participants
    #
    #     # Validate that all required specific molecules are in the input list.
    #     if all(r in input_lipids for r in all_reactants) and all(p in input_lipids for p in all_products):
    #         logger.debug(
    #             f"Validated pre-specified reaction (ID: {base_info.reaction_id}) with all participants present in input.")
    #         # Return a list containing the single, valid, fully specified reaction.
    #         return [base_info._replace(reactants=sorted(list(all_reactants)), products=sorted(list(all_products)))]
    #
    #     # If not all participants are present, this specific reaction is not possible with the given input.
    #     return []

    def _handle_pre_specified_molecular_reaction(self, base_info: Reaction,
                                               input_lipids: Optional[Set[str]],
                                               validate_products: bool = True) -> List[Reaction]:
        """Handles pre-specified reactions, with conditional product validation."""
        all_participants = self.data_access.get_all_participants(base_info.reaction_enzyme_id)
        if not all_participants:
            return []

        all_reactants, all_products = all_participants

        # --- KEY CHANGE: Filter out None values before sorting ---
        valid_reactants = sorted([r for r in all_reactants if r is not None])
        valid_products = sorted([p for p in all_products if p is not None])

        # If validation is turned off, we can just return the reaction as defined.
        if not validate_products:
            logger.debug(f"Returning pre-specified reaction (ID: {base_info.reaction_id}) without validation.")
            return [base_info._replace(reactants=valid_reactants, products=valid_products)]

        # Original validation logic for the standard workflow.
        if input_lipids and all(r in input_lipids for r in valid_reactants) and all(p in input_lipids for p in valid_products):
            logger.debug(f"Validated pre-specified reaction (ID: {base_info.reaction_id}) with all participants present in input.")
            return [base_info._replace(reactants=valid_reactants, products=valid_products)]

        logger.debug(f"Skipping pre-specified reaction (ID: {base_info.reaction_id}) because not all participants were found in the input list.")
        return []

    def _handle_super_reaction_matched_product(self, base_info: Reaction, primary_reactant: str,
                                               matched_product: Union[str, List[str]],
                                               molecule_id: int, input_lipids: Optional[Set[str]],
                                               validate_products: bool = True) -> List[Reaction]:
        """
        Handles direct product matches, with conditional validation.
        The secondary reactant is required and remains blackboxed.
        """
        second_reactant_info = self.data_access.get_second_reactant(molecule_id, base_info.reaction_enzyme_id)

        # --- ROBUSTNESS FIX: Check for valid co-reactant info ---
        if not second_reactant_info or second_reactant_info['cleaned_abbreviation'] is None:
            logger.warning(
                f"Skipping super reaction pair for '{primary_reactant}' because its co-reactant could not be "
                f"identified or is missing an abbreviation. "
                f"Context - Reaction ID: {base_info.reaction_id}"
            )
            return []
        # -----------------------------------------------------------

        if not second_reactant_info:
            return []

        second_reactant_class = self.lipid_parser.parse_lipid(second_reactant_info['cleaned_abbreviation']).headgroup

        reactions = []
        products_to_check = [matched_product] if isinstance(matched_product, str) else matched_product

        for p in products_to_check:
            # KEY CHANGE: Conditionally validate the synthesized product.
            if not validate_products or (input_lipids and p in input_lipids):
                # The "pair" now correctly includes the primary reactant and the blackboxed secondary reactant.
                reactions.append(base_info._replace(reactants=[primary_reactant, second_reactant_class], products=[p]))

        return reactions

    # def _handle_super_reaction_no_match(self, base_info: Reaction, primary_reactant: str, product_class: str,
    #                                     molecule_id: int,
    #                                     input_lipids: Set[str]) -> List[Reaction]:
    #     reaction_enzyme_id = base_info.reaction_enzyme_id
    #     second_reactant_info = self.data_access.get_second_reactant(molecule_id, reaction_enzyme_id)
    #     if not second_reactant_info:
    #         return []
    #
    #     # CORRECTED: Access by key name 'cleaned_abbreviation'
    #     second_reactant_abbr = second_reactant_info['cleaned_abbreviation']
    #     matching_lipids = self._find_matching_lipids_in_input(second_reactant_abbr, input_lipids)
    #
    #     reactions = []
    #     for partner in matching_lipids:
    #         if partner == self.lipid_parser.parse_lipid(second_reactant_abbr).headgroup:
    #             continue
    #
    #         matched_product = self.lipid_matcher.match_sort(partner, product_class, input_lipids)
    #         if isinstance(matched_product, (str, list)):
    #             products = [matched_product] if isinstance(matched_product, str) else matched_product
    #             for p in products:
    #                 reactions.append(base_info._replace(reactants=[primary_reactant, partner], products=[p]))
    #     return reactions

    def _handle_super_reaction_super_match(self, base_info: Reaction, primary_reactant: str,
                                           all_product_classes: Set[str],
                                           molecule_id: int, input_lipids: Optional[Set[str]],
                                           validate_products: bool = True) -> List[Reaction]:
        """
        Handles SUPER_MATCH cases.
        If validating, it finds a specific partner in the input list.
        If not validating, it returns a blackboxed, class-level reaction.
        """
        second_reactant_info = self.data_access.get_second_reactant(molecule_id, base_info.reaction_enzyme_id)
        if not second_reactant_info:
            return []

        second_reactant_abbr = second_reactant_info['cleaned_abbreviation']

        # --- KEY CHANGE: Divert logic for single search mode ---
        if not validate_products:
            # In single search, we cannot find a partner. Return a blackboxed reaction.
            second_reactant_class = self.lipid_parser.parse_lipid(second_reactant_abbr).headgroup

            final_reactants = sorted([primary_reactant, second_reactant_class])
            final_products = sorted(list(all_product_classes))

            logger.debug(f"Single search SUPER_MATCH: returning class-level reaction for {primary_reactant}")
            return [base_info._replace(reactants=final_reactants, products=final_products)]

        # --- Original logic for normal, list-based processing ---
        # This part requires a valid input_lipids set.
        if not input_lipids:
            return []  # Should not happen if validate_products is True, but a safe guard.

        partner_lipids = self._find_matching_lipids_in_input(second_reactant_abbr, input_lipids)

        final_reactions = []
        for partner in partner_lipids:
            if self._is_class(partner):
                continue

            list_of_matched_products_dicts = self.lipid_matcher.super_match_sort(
                primary_reactant, partner, all_product_classes, input_lipids
            )

            for matched_products_dict in list_of_matched_products_dicts:
                # The validation here remains, as it's part of the normal workflow.
                if not any(p in input_lipids for p in matched_products_dict.values()):
                    continue

                final_products = [p_specific if p_specific in input_lipids else p_class
                                  for p_class, p_specific in matched_products_dict.items()]
                final_reactions.append(base_info._replace(
                    reactants=sorted([primary_reactant, partner]),
                    products=sorted(final_products)
                ))
        return final_reactions

    def _handle_headgroup_only_product(self, primary_reactant: str, product_class: str, input_lipids: Set[str]) -> List[
        str]:
        """
        Helper to call match_sort for a simple pair and return a list of matched product strings.
        Now correctly handles str, list, and set return types from the matcher.
        """
        match_result = self.lipid_matcher.match_sort(primary_reactant, product_class, input_lipids)

        # Case 1: Matcher returned a single specific product string.
        if isinstance(match_result, str):
            return [match_result]

        # Case 2: Matcher returned a list or set of specific products.
        # This handles functions like `match_lpc_to_pc_using_input_lipids`.
        if isinstance(match_result, (list, set)):
            return list(match_result)

        # Case 3: Matcher returned an enum (e.g., NO_MATCH) or something else.
        # In this context, it means no valid product was synthesized.
        return []

    def _find_matching_lipids_in_input(self, second_reactant_class: str, input_lipids: Set[str]) -> List[str]:
        """
        Finds specific lipids from the input list that match a given reactant class.
        MODIFIED: Now robustly handles cases where second_reactant_class is None.
        """
        # --- ROBUSTNESS FIX: Handle None input gracefully ---
        if second_reactant_class is None:
            # This can happen if the co-reactant in a super reaction has a NULL abbreviation.
            # We cannot find any matches for a None class, so we return an empty list.
            logger.warning(
                "Attempted to find matching lipids for a None reactant class. This co-reactant will be ignored.")
            return []
        # ---------------------------------------------------

        try:
            headgroup = self.lipid_parser.parse_lipid(second_reactant_class).headgroup
            # The headgroup_map is derived from input_lipids, so if it's empty/None, this will correctly return []
            matching_lipids = list(self.headgroup_map.get(headgroup, []))

            # If no specific lipids match, the "partner" is the class itself (for blackboxing).
            if not matching_lipids:
                return [headgroup]
            return matching_lipids
        except ValueError:
            # If the class name itself is unparsable, return it as-is for blackboxing.
            return [second_reactant_class]

    def _is_fully_specified(self, lipid_name: str) -> bool:
        """
        Determines if a lipid name represents a fully specified molecular species.
        Returns False if it's a class (e.g., "PC"), a species (e.g., "PC(38:4)"),
        or a backbone-only representation (e.g., "Cer(d18:1)").
        """
        try:
            components = self.lipid_parser.parse_lipid(lipid_name)
            # A lipid is not fully specified if it's just a headgroup or if any of its
            # expected chains are missing. Your component models can help here.
            if components.is_headgroup_only:
                return False
            if isinstance(components, SphingoLipidComponents) and not components.n_acyl_chain:
                return False  # This is a backbone-only lipid like Cer(d18:1)
            if components.level == LipidLevel.SPECIES:
                return False

            # Add other similar checks if you have other component types that can be partial
            # e.g., if a Glycerophospholipid can be represented by total carbons like PC(34:1)

        except ValueError:
            # If it can't be parsed, it's definitely not a specific molecule we can handle.
            return False

        return True

    def _is_class(self, lipid_name: str) -> bool:
        """Helper to check if a lipid name is a class name."""
        # Assuming headgroup-only lipids are class representations
        return self.lipid_parser.parse_lipid(lipid_name).is_headgroup_only