import logging
import re
from typing import List, Dict, Tuple, Set, Callable, Optional
from src.models.models import Reaction, ReactionType, ReactionValue, EnzymeInfo
from src.data.data_access import LipidDataAccess
from src.models.lipid_components import SphingoLipidComponents, SphingoidBaseComponents
from src.core.reaction_processor import ReactionProcessor
from src.utils.translation_mapper import TranslationMapper

logger = logging.getLogger(__name__)

class LipidProcessor:
    """
    A class responsible for processing lipids and their associated reactions.

    This class orchestrates the overall lipid processing workflow, including
    retrieving molecule information, processing reactions, and managing the
    reaction dictionaries.
    """
    HEX_PREFIX = "Hex"
    BETA_GLC_PREFIX = "beta-Glc"
    BETA_GAL_PREFIX = "beta-Gal"
    HEX2_PREFIX = "Hex2"
    LAC_PREFIX = "Lac"
    GALA_PREFIX = "Gala"


    def __init__(self, data_access: LipidDataAccess, reaction_processor: ReactionProcessor,
                 translation_mapper: Optional[TranslationMapper] = None):
        """
        Initialize the LipidProcessor.

        Args:
            data_access (LipidDataAccess): An instance of LipidDataAccess for database operations.
            reaction_processor (ReactionProcessor): An instance of ReactionProcessor for reaction processing.
        """
        self.data_access = data_access
        self.reaction_processor = reaction_processor
        self.reaction_dict: Dict[Tuple, ReactionValue] = {}
        self.super_reaction_dict: Dict[Tuple, ReactionValue] = {}
        self.enzyme_dict: Dict[str, EnzymeInfo] = {}
        self.hex_lipid_present = False
        self.translation_mapper: Optional[TranslationMapper] = None

    def set_translation_mapper(self, translation_mapper: TranslationMapper) -> None:
        """Set or update the translation mapper."""
        self.translation_mapper = translation_mapper


    def process_lipids(self, lipids: List[str], organism_id: int) -> Dict[str, Dict[Tuple, ReactionValue]]:
        """Main method for batch processing a list of lipids with mutual constraints."""
        logger.info(f"Starting to process {len(lipids)} lipids for organism {organism_id}")
        input_lipids = set(lipids)

        self.reaction_dict.clear()
        self.super_reaction_dict.clear()
        self.enzyme_dict.clear()

        for lipid in lipids:
            try:
                unique_reactions = self._perform_hierarchical_search(
                    lipid=lipid,
                    organism_id=organism_id,
                    search_action=self._forward_search_action,
                    search_context={'input_lipids': input_lipids, 'validate_products': True}
                )
                for reaction in unique_reactions.values():
                    self._add_reaction_to_dict(reaction, input_lipids)
            except Exception as e:
                logger.error(f"Error processing lipid {lipid} in batch mode: {e}", exc_info=True)
                continue

        return {
            'reaction_dict': self.reaction_dict,
            'super_reaction_dict': self.super_reaction_dict,
            'enzyme_dict': self.enzyme_dict
        }

    def analyze_single_lipid(self, lipid: str, organism_id: int) -> Dict[str, List[Reaction]]:
        """Performs a comprehensive forward and reverse analysis for a single lipid."""
        logger.info(f"Analyzing single lipid connections for: '{lipid}'")

        # Step 1: Perform the raw forward and reverse searches
        forward_reactions = self._perform_hierarchical_search(
            lipid=lipid,
            organism_id=organism_id,
            search_action=self._forward_search_action,
            search_context={'input_lipids': None, 'validate_products': False}
        )

        reverse_reactions = self._perform_hierarchical_search(
            lipid=lipid,
            organism_id=organism_id,
            search_action=self._reverse_search_action,
            search_context={'input_lipids': None, 'validate_products':False}  # Context not needed for reverse search
        )

        # Step 2: Apply the blackboxing rules to the results
        final_forward = self._apply_blackboxing_to_reactions(list(forward_reactions.values()), lipid)
        final_reverse = self._apply_blackboxing_to_reactions(list(reverse_reactions.values()), lipid)

        return {
            "forward_reactions": final_forward,
            "reverse_reactions": final_reverse
        }

    def _blackbox_participant(self, lipid_name: str, searched_lipid: str) -> str:
        """
        Applies backbone-centric blackboxing rules to a lipid participant.
        - The searched lipid itself is never changed.
        - Sphingolipids are reduced to their backbone representation.
        - Other lipids are reduced to their headgroup class.
        """
        # Rule 1: Never alter the lipid that was the original search term.
        # We also check if the name contains the searched lipid's structure (for cases like LPC -> PC)
        if lipid_name == searched_lipid:
            return lipid_name

        try:
            components = self.reaction_processor.lipid_parser.parse_lipid(lipid_name)

            # Rule 2: If it's a sphingolipid, reduce to its backbone form.
            if isinstance(components, (SphingoLipidComponents, SphingoidBaseComponents)):
                if hasattr(components, 'backbone') and components.backbone:
                    base_headgroup = components.headgroup
                    return f"{base_headgroup}({components.backbone})"
                else:  # It's a class-level name like "Cer" TODO: maybe unnecessary??
                    return components.headgroup

            # Rule 3: If not a sphingolipid, reduce to its base headgroup class.
            else:
                base_headgroup = components.headgroup
                return base_headgroup

        except (ValueError, TypeError):
            # If it's not a parsable lipid (e.g., ATP, H2O), return it as is.
            return lipid_name

    def _apply_blackboxing_to_reactions(self, reactions: List[Reaction], searched_lipid: str) -> List[Reaction]:
        """
        Takes a list of Reaction objects and applies blackboxing rules to all
        participants, returning a new list of modified Reaction objects.
        """
        blackboxed_reactions = []
        for reaction in reactions:
            new_reactants = sorted([self._blackbox_participant(r, searched_lipid) for r in reaction.reactants])
            new_products = sorted([self._blackbox_participant(p, searched_lipid) for p in reaction.products])

            # Create a new Reaction object with the modified participants
            new_reaction = reaction._replace(
                reactants=list(dict.fromkeys(new_reactants)),  # Use dict.fromkeys to get unique
                products=list(dict.fromkeys(new_products))
            )
            blackboxed_reactions.append(new_reaction)

        return blackboxed_reactions

    def _perform_hierarchical_search(self,
                                     lipid: str,
                                     organism_id: int,
                                     search_action: Callable,
                                     search_context: Dict) -> Dict[Tuple, Reaction]:
        """
        Performs a generic hierarchical search with robust, granular error handling.
        """
        unique_reactions = {}

        def _add_reactions(reactions: List[Reaction]):
            for r in reactions:
                unique_reactions[r.get_key()] = r

        # --- Step 1 & 2: Molecule and Family Level Search ---
        try:
            molecule_infos = self.data_access.get_molecule_info(lipid)
            if molecule_infos:
                for molecule_id, _, sl_lipid_class in molecule_infos:
                    # Search by direct molecule ID
                    _add_reactions(search_action(lipid, molecule_id, organism_id, **search_context))

                    # Search by family class
                    family_molecule = self.data_access.get_family_molecule(sl_lipid_class)
                    if family_molecule:
                        _add_reactions(search_action(lipid, family_molecule[0], organism_id, **search_context))
            else:
                # --- Step 3: Fallback Search ---
                # This only runs if the lipid itself was not found in the DB.
                components = self.reaction_processor.lipid_parser.parse_lipid(lipid)
                if isinstance(components, SphingoLipidComponents):
                    backbone_family = components.get_backbone_family()
                    backbone_infos = self.data_access.get_molecule_info(backbone_family)
                    for molecule_id, _, _ in backbone_infos:
                        _add_reactions(search_action(lipid, molecule_id, organism_id, **search_context))
        except Exception as e:
            logger.error(f"Error during molecule, family, or fallback search for '{lipid}': {e}", exc_info=True)

        # --- Step 4: Headgroup Level Search (always runs) ---
        try:
            headgroup = self.reaction_processor.lipid_parser.parse_lipid(lipid).headgroup
            if headgroup:
                headgroup_molecules = self.data_access.get_molecule_info(headgroup)
                for hg_molecule_id, _, _ in headgroup_molecules:
                    _add_reactions(search_action(lipid, hg_molecule_id, organism_id, **search_context))
        except Exception as e:
            logger.error(f"Error during headgroup search for '{lipid}': {e}", exc_info=True)

        return unique_reactions

    def _forward_search_action(self, lipid: str, molecule_id: int, organism_id: int,
                               input_lipids: Optional[Set[str]], validate_products: bool) -> List[Reaction]:
        """Action for forward search: gets reactions and synthesizes products."""
        all_reactions = []
        grouped_reactions = self.data_access.get_reactions(molecule_id, organism_id)
        for _, reaction_pairs in grouped_reactions.items():
            processed = self.reaction_processor.process_reaction_group(
                reaction_pairs, lipid, molecule_id, input_lipids, validate_products)
            all_reactions.extend(processed)
        return all_reactions

    def _reverse_search_action(self, target_lipid: str, molecule_id: int, organism_id: int, **kwargs) -> List[Reaction]:
        """Action for reverse search: gets synthesis reactions and synthesizes precursors."""
        synthesis_reactions = []
        grouped_reactions = self.data_access.get_synthesis_reactions_for_molecule(molecule_id, organism_id)
        for _, reaction_pairs in grouped_reactions.items():
            reaction = self._synthesize_precursors(target_lipid, reaction_pairs)
            if reaction:
                synthesis_reactions.append(reaction)
        return synthesis_reactions

    def _synthesize_precursors(self, target_lipid: str, reaction_pairs: List[Reaction]) -> Optional[Reaction]:
        """Synthesizes precursors using the 'reversed-role matching' strategy."""
        if not reaction_pairs: return None

        base_info = reaction_pairs[0]

        # --- KEY CHANGE: Filter out None values before creating the sorted list ---
        all_reactant_templates = sorted(list({
            p.reactants[0] for p in reaction_pairs if p.reactants and p.reactants[0] is not None
        }))
        all_product_templates = sorted(list({
            p.products[0] for p in reaction_pairs if p.products and p.products[0] is not None
        }))
        # --------------------------------------------------------------------------

        synthesized_reactants = []
        target_product_class = self.reaction_processor.lipid_parser.parse_lipid(target_lipid).headgroup

        primary_reactant_template = None
        for rt in all_reactant_templates:
            try:
                reactant_class = self.reaction_processor.lipid_parser.parse_lipid(rt).headgroup
                if self.reaction_processor.lipid_matcher.is_valid_pair(target_product_class, reactant_class):
                    primary_reactant_template = rt
                    break
            except (ValueError, TypeError):  # Catch errors if template is not a valid lipid
                continue

        if primary_reactant_template:
            try:
                matched_result = self.reaction_processor.lipid_matcher.match_sort(
                    reactant=target_lipid, product=primary_reactant_template, input_lipids=set()
                )
                if isinstance(matched_result, (str, list, set)):
                    specific_precursor = list(matched_result) if isinstance(matched_result, (list, set)) else [
                        str(matched_result)]
                    if specific_precursor:
                        synthesized_reactants.append(specific_precursor[0])
                else:
                    synthesized_reactants.append(primary_reactant_template)
            except Exception as e:
                logger.error(
                    f"Error during reversed-role synthesis of '{primary_reactant_template}' from '{target_lipid}': {e}")
                synthesized_reactants.append(primary_reactant_template)

        for rt in all_reactant_templates:
            if rt != primary_reactant_template:
                synthesized_reactants.append(rt)

        final_products = []
        for pt in all_product_templates:
            try:
                product_class = self.reaction_processor.lipid_parser.parse_lipid(pt).headgroup
                if product_class == target_product_class:
                    final_products.append(target_lipid)
                else:
                    final_products.append(pt)
            except (ValueError, TypeError):
                final_products.append(pt)

        return base_info._replace(
            reactants=sorted([r for r in set(synthesized_reactants) if r is not None]),
            products=sorted([p for p in set(final_products) if p is not None])
        )

    def _process_single_lipid(self, lipid: str, organism_id: int, input_lipids: Set[str]) -> None:
        """Process a single lipid."""
        logger.debug(f"Processing lipid: {lipid}")

        # if lipid.startswith(self.HEX2_PREFIX):
        #     self._process_hex2_lipid(lipid, organism_id, input_lipids)
        # elif lipid.startswith(self.HEX_PREFIX):
        #     self._process_hex_lipid(lipid, organism_id, input_lipids)
        # else:
        self._process_non_hex_lipid(lipid, organism_id, input_lipids)

    def _process_hex_lipid(self, lipid: str, organism_id: int, input_lipids: Set[str]) -> None:
        """Process a Hex lipid by creating and processing both Glc and Gal versions."""
        self.hex_lipid_present = True
        try:
            glc_version = lipid.replace(self.HEX_PREFIX, self.BETA_GLC_PREFIX, 1)
            gal_version = lipid.replace(self.HEX_PREFIX, self.BETA_GAL_PREFIX, 1)

            self._process_non_hex_lipid(glc_version, organism_id, input_lipids)
            self._process_non_hex_lipid(gal_version, organism_id, input_lipids)
        finally:
            self.hex_lipid_present = False

    def _process_hex2_lipid(self, lipid: str, organism_id: int, input_lipids: Set[str]) -> None:
        """Process a Hex2 lipid by creating and processing both Lac and Gala versions."""
        self.hex_lipid_present = True
        try:
            lac_version = lipid.replace(self.HEX2_PREFIX, self.LAC_PREFIX, 1)
            gala_version = lipid.replace(self.HEX2_PREFIX, self.GALA_PREFIX, 1)

            self._process_non_hex_lipid(lac_version, organism_id, input_lipids)
            self._process_non_hex_lipid(gala_version, organism_id, input_lipids)
        finally:
            self.hex_lipid_present = False


    def _process_non_hex_lipid(self, lipid: str, organism_id: int, input_lipids: Set[str]) -> None:
        molecule_infos = self.data_access.get_molecule_info(lipid)

        if molecule_infos:
            self._process_molecule_infos(molecule_infos, lipid, organism_id, input_lipids)
        else:
            self._process_lipid_fallback(lipid, organism_id, input_lipids)

        self._process_headgroup_reactions(lipid, organism_id, input_lipids)

    def _process_molecule_infos(self, molecule_infos: List[Tuple[int, str, str]], lipid: str, organism_id: int,
                                input_lipids: Set[str]) -> None:
        logger.debug(f"Processing {len(molecule_infos)} molecule(s) for lipid {lipid}")
        for molecule_id, _, sl_lipid_class in molecule_infos:
            self._safe_process(self._process_molecule_reactions, molecule_id, lipid, organism_id, input_lipids)
            self._safe_process(self._process_family_reactions, sl_lipid_class, lipid, organism_id, input_lipids)

    def _process_lipid_fallback(self, lipid: str, organism_id: int, input_lipids: Set[str]) -> None:
        logger.debug(f"No molecules found for lipid {lipid}. Attempting fallback processing.")
        try:
            lipid_components = self.reaction_processor.lipid_parser.parse_lipid(lipid)

            if isinstance(lipid_components, SphingoLipidComponents):
                self._process_sphingolipid_fallback(lipid_components, lipid, organism_id, input_lipids)
            else:
                logger.debug(f"Lipid {lipid} is not a SphingoLipid. Skipping backbone family processing.")

        except ValueError as e:
            logger.error(f"Error parsing lipid {lipid}: {str(e)}")

    def _process_sphingolipid_fallback(self, lipid_components: SphingoLipidComponents, lipid: str, organism_id: int,
                                       input_lipids: Set[str]) -> None:
        logger.debug(f"Processing SphingoLipid fallback for {lipid}")
        backbone_family = lipid_components.get_backbone_family()
        molecule_infos = self.data_access.get_molecule_info(backbone_family)

        if molecule_infos:
            self._process_molecule_infos(molecule_infos, lipid, organism_id, input_lipids)
        else:
            logger.debug(f"No molecules found for backbone family of {lipid}")

    def _safe_process(self, process_func: Callable, *args) -> None:
        try:
            logger.debug(f"Executing {process_func.__name__}")
            process_func(*args)
        except Exception as e:
            logger.error(f"Error in {process_func.__name__}: {str(e)}")

    def _process_molecule_reactions(self, molecule_id: int, original_lipid: str, organism_id: int,
                                    input_lipids: Set[str]) -> None:
        """
        MODIFIED: Process reactions for a specific molecule by handling them as grouped events.
        """
        # This call now gets the grouped dictionary from the data access layer.
        grouped_reactions = self.data_access.get_reactions(molecule_id, organism_id)

        # Loop over each complete reaction event, identified by its reaction_enzyme_id.
        for reaction_enzyme_id, reaction_pairs in grouped_reactions.items():
            try:  # <-- NEW: Granular error handling for each reaction event.
                # Pass the list of pairs to the new processor method.
                processed_full_reactions = self.reaction_processor.process_reaction_group(
                    reaction_pairs=reaction_pairs,
                    primary_reactant=original_lipid,
                    molecule_id=molecule_id,
                    input_lipids=input_lipids
                )

                # The processor returns a list of fully formed, multi-participant reactions.
                for full_reaction in processed_full_reactions:
                    self._add_reaction_to_dict(full_reaction, input_lipids)

            except Exception as e:
                # If one reaction event fails (e.g., due to an unparseable product), log it and continue.
                logger.error(
                    f"Failed to process reaction event (REID: {reaction_enzyme_id}) for lipid '{original_lipid}'. Error: {e}")
                # The loop will now continue to the next reaction_enzyme_id.

    def _process_family_reactions(self, sl_lipid_class: str, original_lipid: str, organism_id: int,
                                  input_lipids: Set[str]) -> None:
        """
        Process reactions associated with the lipid family.

        This method retrieves and processes reactions for the lipid family of the original lipid.

        Args:
            sl_lipid_class (str): The SwissLipids class of the lipid family.
            original_lipid (str): The original lipid associated with the family.
            organism_id (int): The ID of the organism for which to process reactions.
            input_lipids (Set[str]): A set of all input lipids.
        """
        family_molecule = self.data_access.get_family_molecule(sl_lipid_class)
        if family_molecule:
            # CORRECTED: Access the molecule_id by its key name from the dict_row.
            molecule_id = family_molecule[0]
            self._process_molecule_reactions(molecule_id, original_lipid, organism_id, input_lipids)

    def _process_headgroup_reactions(self, original_lipid: str, organism_id: int, input_lipids: Set[str]) -> None:
        """
        Process reactions associated with the lipid headgroup.

        This method retrieves and processes reactions for the headgroup of the original lipid.

        Args:
            original_lipid (str): The original lipid from which to extract the headgroup.
            organism_id (int): The ID of the organism for which to process reactions.
            input_lipids (Set[str]): A set of all input lipids.
        """
        logger.debug(f"Processing headgroup reactions for lipid: {original_lipid}")
        try:
            headgroup = self.reaction_processor.lipid_parser.parse_lipid(original_lipid).headgroup
            if headgroup:
                headgroup_molecules = self.data_access.get_molecule_info(headgroup)
                for hg_molecule_id, _, _ in headgroup_molecules:
                    self._process_molecule_reactions(hg_molecule_id, original_lipid, organism_id, input_lipids)
            else:
                logger.debug(f"No headgroup found for lipid: {original_lipid}")
        except Exception as e:
            logger.error(f"Error processing headgroup reactions for {original_lipid}: {str(e)}")

    def find_products_for_single_lipid(self, lipid: str, organism_id: int) -> List[Reaction]:
        """
        Finds all possible reactions and synthesizes products for a single lipid,
        without the constraint that products must be in a predefined list.
        """
        logger.info(f"Initiating unconstrained search for lipid: '{lipid}' in organism {organism_id}")

        # Use a set to store unique reactions by their key to avoid duplicates
        unique_reactions = {}

        try:
            molecule_infos = self.data_access.get_molecule_info(lipid)

            # Primary processing based on the specific lipid's molecule ID
            if molecule_infos:
                for molecule_id, _, sl_lipid_class in molecule_infos:
                    reactions = self._process_single_lipid_unconstrained(lipid, molecule_id, organism_id)
                    for r in reactions:
                        unique_reactions[r.get_key()] = r

                    # Also process reactions defined at the family level
                    family_reactions = self._process_family_reactions_unconstrained(sl_lipid_class, lipid, organism_id)
                    for r in family_reactions:
                        unique_reactions[r.get_key()] = r
            else:
                # If the lipid is not in the DB, try to find reactions via its structural components
                fallback_reactions = self._process_lipid_fallback_unconstrained(lipid, organism_id)
                for r in fallback_reactions:
                    unique_reactions[r.get_key()] = r

            # Finally, process reactions defined at the headgroup level
            headgroup_reactions = self._process_headgroup_reactions_unconstrained(lipid, organism_id)
            for r in headgroup_reactions:
                unique_reactions[r.get_key()] = r

        except Exception as e:
            logger.error(f"An error occurred during single lipid search for '{lipid}': {e}", exc_info=True)
            raise

        return list(unique_reactions.values())

    def _process_single_lipid_unconstrained(self, lipid: str, molecule_id: int, organism_id: int) -> List[Reaction]:
        """Helper to process a molecule's reactions without input constraints."""
        all_reactions = []
        grouped_reactions = self.data_access.get_reactions(molecule_id, organism_id)

        for reaction_enzyme_id, reaction_pairs in grouped_reactions.items():
            try:
                # Call the core processor with validation turned OFF
                processed_reactions = self.reaction_processor.process_reaction_group(
                    reaction_pairs=reaction_pairs,
                    primary_reactant=lipid,
                    molecule_id=molecule_id,
                    input_lipids=None,
                    validate_products=False
                )
                all_reactions.extend(processed_reactions)
            except Exception as e:
                logger.error(
                    f"Failed to process unconstrained reaction group {reaction_enzyme_id} for lipid {lipid}: {e}")

        return all_reactions

    def _process_family_reactions_unconstrained(self, sl_lipid_class: str, lipid: str, organism_id: int) -> List[
        Reaction]:
        """Helper to process family reactions without constraints."""
        family_molecule = self.data_access.get_family_molecule(sl_lipid_class)
        if family_molecule:
            molecule_id = family_molecule[0]
            return self._process_single_lipid_unconstrained(lipid, molecule_id, organism_id)
        return []

    def _process_lipid_fallback_unconstrained(self, lipid: str, organism_id: int) -> List[Reaction]:
        """Helper for fallback processing without constraints."""
        try:
            lipid_components = self.reaction_processor.lipid_parser.parse_lipid(lipid)
            if isinstance(lipid_components, SphingoLipidComponents):
                backbone_family = lipid_components.get_backbone_family()
                molecule_infos = self.data_access.get_molecule_info(backbone_family)

                all_reactions = []
                for molecule_id, _, _ in molecule_infos:
                    reactions = self._process_single_lipid_unconstrained(lipid, molecule_id, organism_id)
                    all_reactions.extend(reactions)
                return all_reactions
        except ValueError as e:
            logger.error(f"Error parsing lipid for fallback search {lipid}: {e}")
        return []

    def _process_headgroup_reactions_unconstrained(self, lipid: str, organism_id: int) -> List[Reaction]:
        """Helper for headgroup processing without constraints."""
        try:
            headgroup = self.reaction_processor.lipid_parser.parse_lipid(lipid).headgroup
            if headgroup:
                headgroup_molecules = self.data_access.get_molecule_info(headgroup)
                all_reactions = []
                for hg_molecule_id, _, _ in headgroup_molecules:
                    reactions = self._process_single_lipid_unconstrained(lipid, hg_molecule_id, organism_id)
                    all_reactions.extend(reactions)
                return all_reactions
        except Exception as e:
            logger.error(f"Error processing headgroup reactions for {lipid}: {e}")
        return []

    def _add_reaction_to_dict(self, reaction: Reaction, input_lipids: Set[str]) -> None:
        """
        MODIFIED: Validation logic is updated to check a list of products.
        """
        # This check is now more complex because a "self-reaction" might be valid
        # if other reactants/products are involved. We will rely on the new `get_key`
        # and validation in `_add_single_reaction_to_dict` to handle this correctly.
        # For now, we remove the simple check. A more robust check could be added later if needed.
        # if reaction.reaction_type == ReactionType.REGULAR and reaction.reactants[0] in reaction.products:
        #    ...

        target_dict = self.super_reaction_dict if reaction.reaction_type == ReactionType.SUPER else self.reaction_dict

        # The core logic of calling the single reaction handler remains the same.
        self._add_single_reaction_to_dict(reaction, target_dict, input_lipids, self.translation_mapper)

    def _add_single_reaction_to_dict(self, reaction: Reaction, target_dict: Dict, input_lipids: Set[str],
                                     translation_mapper: Optional[TranslationMapper] = None) -> None:
        """
        MODIFIED: Validates that at least one of the specific, non-class-level products
        of a reaction exists in the input lipid set before adding it.
        """
        # New validation logic for multi-product reactions
        has_valid_product_in_input = False
        for p in reaction.products:
            # A product is considered for validation if it's a specific molecular species
            # (i.e., not just a class name like 'DG') AND it's present in the input lipids.
            is_class_level = self.reaction_processor.lipid_parser.parse_lipid(p).is_headgroup_only
            if not is_class_level and p in input_lipids:
                has_valid_product_in_input = True
                break # Found at least one valid product, no need to check further.

        if has_valid_product_in_input:
            if translation_mapper:
                # The reverse_translate_reaction method will also need a small update
                # to handle the list of products.
                translated_reactions = translation_mapper.reverse_translate_reaction(reaction, input_lipids)
                for trans_reaction in translated_reactions:
                    self._add_reaction_to_dict_internal(trans_reaction, target_dict)
            else:
                logger.warning(f"Failed to use reverse translation mapper. Reaction ID: {reaction.reaction_id}")
                self._add_reaction_to_dict_internal(reaction, target_dict)
        else:
            logger.debug(f"Skipping reaction (ID: {reaction.reaction_id}) as none of its specific products "
                         f"({[p for p in reaction.products if not self.reaction_processor.lipid_parser.parse_lipid(p).is_headgroup_only]}) "
                         f"were found in the input lipid list.")



# #TODO: fix to use translation mapper correctly and fix translation mapper
#             self._add_enzyme_info(reaction.enzyme_name, reaction.uniprot_id)
#             key = reaction.get_key()
#             reaction_value = ReactionValue(reaction.rhea_id, reaction)
#
#             if key in target_dict:
#                 existing_value = target_dict[key]
#                 updated_value = existing_value.update_rhea_id(reaction_value.rhea_id)
#                 target_dict[key] = updated_value
#             else:
#                 target_dict[key] = reaction_value
#         else:
#             logger.warning(f"Attempted to add reaction with product {reaction.product} not in input list. "
#                            f"Reaction ID: {reaction.reaction_id}")

    def _add_reaction_to_dict_internal(self, reaction: Reaction, target_dict: Dict) -> None:
        """
        Internal helper method to add a single reaction to the dictionary.

        Args:
            reaction: The reaction to add
            target_dict: Dictionary to add the reaction to
        """
        self._add_enzyme_info(reaction.enzyme_name, reaction.uniprot_id, reaction.subcellular_locations)
        key = reaction.get_key()
        new_reaction_value = ReactionValue(reaction.rhea_id, reaction)

        if key in target_dict:
            existing_value = target_dict[key]
            merged_value = existing_value.merge_with(new_reaction_value)
            target_dict[key] = merged_value
        else:
            target_dict[key] = new_reaction_value

    def _get_hex_version_of_lipid(self, lipid: str) -> str:
        if lipid.startswith("beta-Glc") or (lipid.startswith("beta-Gal") and not lipid.startswith("beta-Gala")):
            return "Hex" + lipid[8:]
        elif lipid.startswith("Lac"):
            return "Hex2" + lipid[3:]
        elif lipid.startswith("Gala"):
            return "Hex2" + lipid[4:]
        return lipid

    def _convert_reaction_to_hex(self, reaction: Reaction) -> Reaction:
        """
        Convert the reactants and product of a reaction back to Hex form if they were originally Hex lipids.

        Args:
            reaction (Reaction): The reaction to convert.

        Returns:
            Reaction: A new Reaction object with converted lipid names where necessary.
        """
        return reaction._replace(
            reactants=[self._convert_to_hex(r) for r in reaction.reactants],
            product=self._convert_to_hex(reaction.product)
        )

    def _convert_to_hex(self, lipid: str) -> str:
        """
        Convert a lipid back to Hex form if it was originally a Hex lipid.

        Args:
            lipid (str): A lipid name.

        Returns:
            str: The converted lipid name or the original lipid if no conversion was needed.
        """
        if lipid.startswith("beta-Glc") or (lipid.startswith("beta-Gal") and not lipid.startswith("beta-Gala")):
            return "Hex" + lipid[8:]
        return lipid

    def _add_enzyme_info(self, enzyme_name: str, uniprot_id: str, locations_str: Optional[str] = None) -> None:
        """
        Add or update enzyme information in the enzyme dictionary.

        This method adds a new enzyme to the enzyme dictionary if it doesn't exist,
        or updates the existing entry if the enzyme is already in the dictionary.

        Args:
            enzyme_name (str): The name of the enzyme.
            uniprot_id (str): The UniProt ID associated with the enzyme.
        """
        if enzyme_name not in self.enzyme_dict:
            # First time seeing this enzyme, create a new entry
            self.enzyme_dict[enzyme_name] = EnzymeInfo(
                enzyme_name=enzyme_name,
                uniprot_id=uniprot_id,
                subcellular_locations=self._parse_and_clean_locations(locations_str)
            )

    def _parse_and_clean_locations(self, locations_str: Optional[str]) -> Optional[str]:
        """
        Parses the raw subcellular location string from the database into a clean,
        pipe-separated list that preserves form information when present.

        FINAL VERSION: Handles multiple patterns:
        - [Isoform X]: locations...
        - [Named form]: locations...
        - Direct locations without form names

        Returns format like:
        - "Isoform 1: Membrane, Cytoplasm | Neutral ceramidase: Cell membrane, Golgi apparatus"
        - "Cell membrane, Cytoplasm | Isoform 2: Nucleus"
        - "Cell membrane | Cytoplasm" (simple case)
        """
        if not locations_str:
            return None

        # STEP 1: Remove evidence codes (everything in curly braces) globally
        cleaned_str = re.sub(r'\s*\{[^}]*\}\s*', '', locations_str)

        # STEP 2: Split by "SUBCELLULAR LOCATION:" boundaries FIRST
        # This regex finds "SUBCELLULAR LOCATION:" followed by content until the next one or end
        location_sections = re.findall(
            r'SUBCELLULAR LOCATION:\s*(.*?)(?=SUBCELLULAR LOCATION:|$)',
            cleaned_str,
            re.DOTALL
        )

        if not location_sections:
            return None

        # STEP 3: Remove Notes from each section individually
        cleaned_sections = []
        for section in location_sections:
            # Remove Note sections from this specific section only
            section_cleaned = re.sub(r'\s*Note=.*?$', '', section, flags=re.DOTALL)
            section_cleaned = section_cleaned.strip()
            if section_cleaned:
                cleaned_sections.append(section_cleaned)

        if not cleaned_sections:
            return None

        # STEP 4: Check if any section has form information (broader than just "Isoform")
        has_forms = any(re.search(r'\[[^]]+\]:', section) for section in cleaned_sections)

        if has_forms or len(cleaned_sections) > 1:
            # Multiple sections or named forms - use the multi-section parser
            return self._parse_multiple_isoform_sections(cleaned_sections)
        else:
            # Single section with no form name - use simple parser
            combined_text = " ".join(cleaned_sections)
            return self._parse_simple_locations(combined_text)

    def _parse_multiple_isoform_sections(self, sections: List[str]) -> Optional[str]:
        """
        Parse multiple sections that may contain isoform/form information.
        Each section is a complete SUBCELLULAR LOCATION entry.

        Handles cases like:
        - [Isoform 1]: locations...
        - [Neutral ceramidase]: locations...
        - [Neutral ceramidase soluble form]: locations...
        - Direct locations without a form name
        """
        isoform_entries = []

        for section in sections:
            section = section.strip()
            if not section:
                continue

            # Remove any leading semicolons or periods from splitting
            section = re.sub(r'^[;.]\s*', '', section)

            # Check if this section starts with a form identifier like [Something]:
            form_match = re.match(r'\[([^]]+)\]:\s*(.*)', section, re.DOTALL)
            if form_match:
                # This section has a named form (isoform, enzyme variant, etc.)
                form_name = form_match.group(1)
                location_text = form_match.group(2)
            else:
                # This section has no form identifier, just locations
                # Don't add any prefix - just use the locations directly
                form_name = None
                location_text = section

            # Extract locations for this section
            locations = self._extract_locations_from_text(location_text)

            if locations:
                # Remove duplicates and sort
                unique_locations = sorted(set(locations))

                if form_name:
                    # Named form: "Isoform 1: Location1, Location2"
                    isoform_entry = f"{form_name}: {', '.join(unique_locations)}"
                else:
                    # Unnamed section: just "Location1, Location2" (no prefix)
                    isoform_entry = ', '.join(unique_locations)

                isoform_entries.append(isoform_entry)

        return " | ".join(isoform_entries) if isoform_entries else None

    def _parse_simple_locations(self, text: str) -> Optional[str]:
        """
        Parse locations without isoform information.
        Returns format: "Location1 | Location2 | Location3"
        """
        locations = self._extract_locations_from_text(text)
        if not locations:
            return None

        # Remove duplicates and sort
        unique_locations = sorted(set(locations))
        return " | ".join(unique_locations)

    def _extract_locations_from_text(self, text: str) -> List[str]:
        """
        Extract clean location names from a text section.
        Handles complex punctuation patterns correctly.
        """
        if not text:
            return []

        # Clean up the text first
        text = text.strip()

        # Remove trailing semicolons and periods that might interfere
        text = re.sub(r'[;.]\s*$', '', text)

        # Protect common abbreviations from period splitting
        text = text.replace('e.g.', 'EGSUB').replace('i.e.', 'IESUB')

        # Split by periods followed by whitespace
        # This handles cases like "Location1. Location2" and "Location1.Location2"
        statements = re.split(r'\.\s*', text)

        # Restore abbreviations
        statements = [stmt.replace('EGSUB', 'e.g.').replace('IESUB', 'i.e.') for stmt in statements]

        locations = []

        for statement in statements:
            statement = statement.strip()
            if not statement:
                continue

            # Take the part before the first semicolon (main location)
            main_location = statement.split(';')[0].strip()

            if main_location:
                # Clean the location name
                cleaned_location = self._clean_location_name(main_location)
                if cleaned_location:
                    locations.append(cleaned_location)

        return locations

    def _clean_location_name(self, location: str) -> Optional[str]:
        """
        Clean individual location names, removing descriptive terms.
        """
        if not location:
            return None

        location = location.strip()

        # Skip empty or very short strings
        if len(location) < 2:
            return None

        # Skip strings that are mostly punctuation or numbers
        if re.match(r'^[^a-zA-Z]*$', location):
            return None

        # Remove trailing punctuation
        location = re.sub(r'[.,:;]+$', '', location).strip()

        # Skip descriptive terms that aren't actual locations
        skip_patterns = [
            r'^Single-pass membrane protein$',
            r'^Multi-pass membrane protein$',
            r'^Peripheral membrane protein$',
            r'^Single-pass type II membrane protein$',
            r'^membrane protein$',
            r'^Note$',
            r'^\d+$'
        ]

        for pattern in skip_patterns:
            if re.match(pattern, location, re.IGNORECASE):
                return None

        return location
