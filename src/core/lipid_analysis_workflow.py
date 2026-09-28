import logging
import os
import json
from typing import List, Dict, Tuple, Union, Set, Optional
import pandas as pd
from dataclasses import dataclass, field
import re
from pathlib import Path

from src.core.lipid_preprocessor import LipidPreprocessor
from src.core.lipid_processor import LipidProcessor
from src.data.data_access import DatabaseLipidDataAccess
from src.matching.lipid_matcher import LipidMatcher
from src.parsing.lipid_parser import LipidParser, LipidParserFactory, LipidComponentCache
from src.models.models import ReactionValue, OutputContext, SingleLipidSearchContext, SingleEnzymeSearchContext
from src.output.output_generators import (create_output_generators, TranslateOutputGenerator,
                                          SingleLipidSearchOutputGenerator, SingleEnzymeSearchOutputGenerator)
from src.data.config_loader import ConfigLoader
from src.core.reaction_processor import ReactionProcessor
from src.utils.lipid_translator import LipidTranslator
from src.core.lipid_network_manager import LipidNetworkManager
from src.utils.translation_mapper import TranslationMapper
from src.core.lipid_preprocessor import LipidPreprocessor, PreprocessingResult
from src.utils.text_utils import split_reaction_text
from src.utils.resolver import MoleculeResolver

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class LipidAnalysisWorkflow:
    def __init__(self, config: Union[str, ConfigLoader], output_dir: str, debug_mode: bool = False):
        # Initialize config and basic settings
        self.config_loader = config if isinstance(config, ConfigLoader) else ConfigLoader(config)
        self.output_dir = output_dir
        self.debug_mode = debug_mode

        # Initialize data access
        self.data_access = DatabaseLipidDataAccess(self.config_loader)

        # Initialize parser components
        lipid_component_cache = LipidComponentCache()
        lipid_parser_factory = LipidParserFactory(self.config_loader)
        self.lipid_parser = LipidParser(
            factory=lipid_parser_factory,
            components_cache=lipid_component_cache
        )

        # Initialize translator
        self.lipid_translator = LipidTranslator(self.config_loader, self.lipid_parser)
        self.translation_mapper = None

        # Path to the directory of this file (i.e., src/core/)
        current_dir = Path(__file__).resolve().parent
        # Go up two levels to project root
        project_root = current_dir.parent.parent
        resolutions_path = project_root / "db_updates" / "api_fallback_resolutions.yaml"

        self.molecule_resolver = MoleculeResolver(
            data_access=self.data_access,
            lipid_parser=self.lipid_parser,
            lipid_translator=self.lipid_translator,
            resolutions_filepath=resolutions_path
        )

        # Initialize preprocessor
        self.lipid_preprocessor = LipidPreprocessor(self.lipid_parser)

        # Initialize matcher and processors
        self.lipid_matcher = LipidMatcher(self.lipid_parser, self.config_loader)
        self.reaction_processor = ReactionProcessor(self.lipid_matcher, self.lipid_parser, self.data_access)
        self.lipid_processor = LipidProcessor(self.data_access, self.reaction_processor)

        # Initialize output generators
        self.output_generators = create_output_generators(self.config_loader, debug_mode=self.debug_mode)
        self.lipid_network_manager = None

        output_paths = self.config_loader.get_output_paths()
        self.single_lipid_search_generator = SingleLipidSearchOutputGenerator(
            output_file=output_paths.get('single_lipid_report', 'single_lipid_search_report.tsv')
        )
        self.single_enzyme_search_generator = SingleEnzymeSearchOutputGenerator(
            output_file=output_paths.get('single_enzyme_report', 'single_enzyme_search_report.tsv')
        )

    def _generate_single_lipid_search_output(self, search_results: Dict[str, Dict], base_name: str) -> None:
        """Orchestrates the output generation for a single lipid search."""
        logger.info(f"Generating single lipid search output with base name: {base_name}")

        # We still need to prefetch Rhea IDs for the report.
        all_reactions = []
        for result in search_results.values():
            # Check for the correct keys: 'forward_reactions' and 'reverse_reactions'
            if 'forward_reactions' in result:
                all_reactions.extend(result['forward_reactions'])
            if 'reverse_reactions' in result:
                all_reactions.extend(result['reverse_reactions'])

        # Create ReactionValue objects to reuse the existing prefetcher
        reaction_values = [ReactionValue(r.rhea_id, r) for r in all_reactions]

        # --- ADD DEBUGGING HERE ---
        if self.debug_mode:
            generalized_ids_found = [rv.reaction.generalized_from_ids for rv in reaction_values if
                                     rv.reaction.generalized_from_ids]
            logger.debug(f"Found {len(generalized_ids_found)} reactions with generalized_from_ids field.")
            if generalized_ids_found:
                logger.debug(f"Example generalized_from_ids: {generalized_ids_found[:5]}")
        # ---------------------------

        specific_rhea_map = self._prefetch_specific_rhea_ids(reaction_values)

        # --- ADD MORE DEBUGGING HERE ---
        if self.debug_mode:
            logger.debug(f"Prefetched specific Rhea map contains {len(specific_rhea_map)} entries.")
            if specific_rhea_map:
                # Log a few example entries from the map
                logger.debug(f"Example specific_rhea_map entries: {dict(list(specific_rhea_map.items())[:5])}")
        # ------------------------------

        # Create the context and generate the file.
        context = SingleLipidSearchContext(
            search_results=search_results,
            base_name=base_name,
            output_dir=self.output_dir,
            debug_mode=self.debug_mode,
            specific_rhea_map=specific_rhea_map
        )
        self.single_lipid_search_generator.generate(context)

    def _generate_single_enzyme_search_output(self, search_results: Dict, base_name: str) -> None:
        """Orchestrates the output generation for a single enzyme search."""
        logger.info(f"Generating single enzyme search output with base name: {base_name}")

        context = SingleEnzymeSearchContext(
            search_results=search_results,
            base_name=base_name,
            output_dir=self.output_dir,
            debug_mode=self.debug_mode
        )
        self.single_enzyme_search_generator.generate(context)

    def search_products_by_lipid(self, lipids: Union[str, List[str]], organism_id: int,
                                 base_name: str = "lipid_search") -> Dict[str, Dict]:
        """
        Performs an independent, unconstrained search for each lipid provided.
        MODIFIED: Now accepts a base_name and generates an output file.
        """
        logger.info(f"Executing single lipid search for {lipids} in organism {organism_id}.")

        lipids_to_process = [lipids] if isinstance(lipids, str) else lipids
        final_results = {}

        for lipid_name in lipids_to_process:
            try:
                cleaned_lipid = self.lipid_translator.clean_lipid_abbreviation(lipid_name)
                reactions = self.lipid_processor.find_products_for_single_lipid(cleaned_lipid, organism_id)
                products = set(p for r in reactions for p in r.products)

                final_results[lipid_name] = {
                    'cleaned_input': cleaned_lipid, 'reactions': reactions,
                    'products': sorted(list(products)), 'reaction_count': len(reactions),
                    'product_count': len(products)
                }
            except Exception as e:
                logger.error(f"Single lipid search failed for '{lipid_name}': {e}", exc_info=True)
                final_results[lipid_name] = {'error': str(e)}

        # Generate the output file after all searches are complete.
        self._generate_single_lipid_search_output(final_results, base_name)

        return final_results

    def analyze_single_lipid(self, lipids: Union[str, List[str]], organism_id: int,
                             base_name: str = "lipid_analysis") -> Dict[str, Dict]:
        """
        Performs a comprehensive forward and reverse analysis for each lipid provided.
        MODIFIED: Now includes translation for each input lipid.
        """
        logger.info(f"Executing unified analysis for {lipids} in organism {organism_id}.")

        lipids_to_process = [lipids] if isinstance(lipids, str) else lipids
        final_results = {}

        for lipid_name in lipids_to_process:
            try:
                # --- NEW TRANSLATION STEP ---
                # We process each lipid individually to get its set of translated names.
                # We pass `translate=True` to our existing helper.
                preprocessing_result, translation_tracker = self._preprocess_and_translate([lipid_name], translate=True)

                # Create the mapper for this specific lipid's translations.
                mapper = TranslationMapper(translation_tracker)


                # `all_valid_lipids` will now contain the original lipid plus all its translated forms.
                lipids_for_this_search = preprocessing_result.all_valid_lipids
                logger.debug(f"Input '{lipid_name}' translated to: {lipids_for_this_search}")
                # --------------------------

                # We will aggregate the results from all translated forms of this one lipid.
                all_forward_reactions = {}
                all_reverse_reactions = {}

                for translated_lipid in lipids_for_this_search:
                    # Perform the analysis on each translated variant.
                    analysis = self.lipid_processor.analyze_single_lipid(translated_lipid, organism_id)

                    # Merge the results, ensuring no duplicates.
                    for f_reac in analysis.get('forward_reactions', []):
                        all_forward_reactions[f_reac.get_key()] = f_reac
                    for r_reac in analysis.get('reverse_reactions', []):
                        all_reverse_reactions[r_reac.get_key()] = r_reac

                # --- Step 3: Reverse-translate the final reactions ---
                # We need the original name as a single-element set for the mapper.
                original_input_set = {lipid_name}

                final_forward = []
                for reaction in all_forward_reactions.values():
                    # The mapper returns a list of possible reverse-translated reactions.
                    # We take the first one as it's sufficient for this context.
                    rev_translated = mapper.reverse_translate_reaction(
                        reaction, original_input_set, validate_products=False
                    )
                    if rev_translated:
                        final_forward.append(rev_translated[0])

                final_reverse = []
                for reaction in all_reverse_reactions.values():
                    rev_translated = mapper.reverse_translate_reaction(
                        reaction, original_input_set, validate_products=False
                    )
                    if rev_translated:
                        final_reverse.append(rev_translated[0])

                # Store the final, aggregated results under the original input name.
                final_results[lipid_name] = {
                    "forward_reactions": final_forward,
                    "reverse_reactions": final_reverse,
                }

            except Exception as e:
                logger.error(f"Unified analysis failed for '{lipid_name}': {e}", exc_info=True)
                final_results[lipid_name] = {'error': str(e)}

        self._generate_single_lipid_search_output(final_results, base_name)
        return final_results

    def search_reactions_by_enzyme(self, enzyme_identifier: str, organism_id: int, base_name: str = "enzyme_search") -> \
    Dict[str, Union[str, int, bool, List]]:
        """
        Performs a single enzyme search with robust parsing, abbreviation, and de-duplication.
        MODIFIED: No longer filters out common metabolites.
        """
        logger.info(f"Executing single enzyme search for '{enzyme_identifier}' in organism {organism_id}.")

        try:
            enzyme_reactions = self.data_access.get_reactions_by_enzyme(enzyme_identifier, organism_id)

            if not enzyme_reactions:
                output_payload = {'enzyme_identifier': enzyme_identifier, 'found': False, 'reactions': []}
                self._generate_single_enzyme_search_output(output_payload, base_name)
                return output_payload

            aggregated_reactions = {}

            for reaction_info in enzyme_reactions:
                reaction_text = reaction_info.get('reaction_text', '')
                rhea_id = reaction_info.get('rhea_id')  # Get the rhea_id for context

                parsed_components, is_reversible = split_reaction_text(reaction_text)

                if not parsed_components:
                    abbreviated_text = f"COULD_NOT_PARSE: {reaction_text}"
                else:
                    sl_reactant_names, sl_product_names = parsed_components
                    arrow = " <=> " if is_reversible else " => "

                    reactant_classes = [self.molecule_resolver.resolve_to_class(r, rhea_id) for r in sl_reactant_names]
                    product_classes = [self.molecule_resolver.resolve_to_class(p, rhea_id) for p in sl_product_names]

                    abbreviated_text = f"{' + '.join(sorted(reactant_classes))} {arrow} {' + '.join(sorted(product_classes))}"

                # --- NEW EVIDENCE AGGREGATION LOGIC ---
                if abbreviated_text not in aggregated_reactions:
                    # First time we've seen this abbreviated reaction. Initialize it.
                    aggregated_reactions[abbreviated_text] = {
                        'abbreviated_reaction_text': abbreviated_text,
                        'original_reaction_text': reaction_text,  # Keep the first one as an example
                        'rhea_id_example': str(reaction_info['rhea_id']) if reaction_info.get('rhea_id') is not None else None,
                        'consolidation_count': 1
                    }
                else:
                    # We've seen this abbreviated reaction before. Aggregate evidence.
                    entry = aggregated_reactions[abbreviated_text]
                    entry['consolidation_count'] += 1

                    # If the current entry has no Rhea ID, but this new one does, adopt it.
                    if not entry['rhea_id_example'] and reaction_info.get('rhea_id'):
                        entry['rhea_id_example'] = str(reaction_info['rhea_id'])

            processed_reactions = list(aggregated_reactions.values())

            if not processed_reactions:
                output_payload = {
                    'enzyme_identifier': enzyme_identifier, 'found': False,
                    'reactions': [], 'reaction_count': 0
                }
                self._generate_single_enzyme_search_output(output_payload, base_name)
                return output_payload

            # Get enzyme info from the FIRST original reaction, as it's shared by all.
            first_original_reaction = enzyme_reactions[0]
            output_payload = {
                'enzyme_identifier': enzyme_identifier,
                'enzyme_name': first_original_reaction['enzyme_name'],
                'uniprot_id': first_original_reaction['uniprot_id'],
                'found': True,
                'reactions': processed_reactions,  # Pass the new aggregated list
                'reaction_count': len(processed_reactions)
            }

            self._generate_single_enzyme_search_output(output_payload, base_name)
            return output_payload

        except Exception as e:
            logger.error(f"Single enzyme search failed for '{enzyme_identifier}': {e}", exc_info=True)
            raise

    def process_lipid_list(self, lipids: List[str], organism_id: int) -> Dict[str, Dict[Tuple, ReactionValue]]:
        """
        Process a list of lipids for a given organism.

        Args:
            lipids (List[str]): A list of lipids to process.
            organism_id (int): The ID of the organism for which to process the lipids.

        Returns:
            Dict[str, Dict[Tuple, ReactionValue]]: A dictionary containing the processed
            reaction dictionaries (regular and super reactions).

        Raises:
            Exception: If an error occurs during lipid processing.
        """
        logger.info(f"Starting lipid analysis workflow for organism {organism_id}")
        try:
            return self.lipid_processor.process_lipids(lipids, organism_id)
        except Exception as error:
            logger.error(f"Error in process_lipid_list: {str(error)}")
            raise

    def generate_outputs(self, results: Dict[str, Union[Dict[Tuple, ReactionValue], Set[str]]], translate: bool,
                         base_name: str, headgroup_map: Dict[str, Set[str]]) -> None:
        """
        Generate output files based on the processed results.

        Args:
            results (Dict[str, Union[Dict[Tuple, ReactionValue], Set[str]]]): The processed results.
            translate (bool): Whether translation was performed.
            base_name (str): The base name for output files.

        Raises:
            Exception: If an error occurs during output generation.
        """
        logger.info("Generating outputs")

        all_reaction_values = list(results.get('reaction_dict', {}).values()) + list(
            results.get('super_reaction_dict', {}).values())
        specific_rhea_map = self._prefetch_specific_rhea_ids(all_reaction_values)

        # Add the map to the results dictionary to be passed down.
        results['specific_rhea_map'] = specific_rhea_map

        try:
            context = OutputContext(
                reaction_dict=results['reaction_dict'],
                super_reaction_dict=results['super_reaction_dict'],
                input_lipids=results['original_input_lipids'],
                enzyme_dict=results.get('enzyme_dict', {}),
                headgroup_map=headgroup_map,
                translation_tracker=results.get('translation_tracker'),
                specific_rhea_map=results.get('specific_rhea_map', {}),
                base_name=base_name,
                output_dir=self.output_dir,
                debug_mode=self.debug_mode
            )
            for generator in self.output_generators:
                if not isinstance(generator, TranslateOutputGenerator) or translate:
                    generator.generate(context)

        except Exception as error:
            logger.error(f"Error in generate_outputs: {str(error)}")
            raise

    def _prefetch_specific_rhea_ids(self, all_reaction_values: List[ReactionValue]) -> Dict[int, str]:
        """
        NEW HELPER: Queries the database once to get all Rhea IDs for the specific reactions
        referenced in the 'generalized_from_ids' field.
        """
        all_specific_ids = set()
        for rv in all_reaction_values:
            generalized_ids_str = rv.reaction.generalized_from_ids
            if generalized_ids_str:
                ids = generalized_ids_str.strip('|').split('|')
                all_specific_ids.update(int(i) for i in ids if i.isdigit())

        if not all_specific_ids:
            return {}

        query = "SELECT reaction_id, rhea_id FROM lipograph.reactions WHERE reaction_id = ANY(%s)"
        try:
            # Use this instance's already-existing cursor
            self.data_access.cursor.execute(query, (list(all_specific_ids),))
            rhea_map = {row['reaction_id']: row['rhea_id'] for row in self.data_access.cursor.fetchall() if
                        row.get('rhea_id')}
            return rhea_map
        except Exception as e:
            logger.error(f"Failed to prefetch specific Rhea IDs: {e}")
            return {}

    def _process_lipid_network(self, results: Dict[str, Union[Dict[Tuple, ReactionValue], Set[str]]], base_name: str) -> str:
        """
        Process the lipid network and save it as a JSON file.

        Args:
            results (Dict[str, Union[Dict[Tuple, ReactionValue], Set[str]]]): The processed results.
            base_name (str): The base name for the output file.

        Returns:
            str: The path to the saved network JSON file.

        Raises:
            Exception: If an error occurs during network processing.
        """
        try:
            json_data = self._initialize_network_json()

            self.lipid_network_manager = LipidNetworkManager(
                json_data,
                results['reaction_dict'],
                results['super_reaction_dict'],
                self.config_loader,
                self.lipid_parser
            )

            self.lipid_network_manager.process_network()

            network_file_path = os.path.join(self.output_dir, f"{base_name}_network.json")
            with open(network_file_path, 'w') as f:
                json.dump(self.lipid_network_manager.json_data, f, indent=2)

            logger.info(f"Network JSON saved to {network_file_path}")

            return network_file_path
        except Exception as error:
            logger.error(f"Error in _process_lipid_network: {str(error)}")
            raise

    @staticmethod
    def _initialize_network_json() -> Dict:
        """
        Initialize the JSON data structure for the network.

        Returns:
            Dict: The initialized JSON data structure.
        """
        return {
            "nodes": [],
            "edges": [],
            "layers": [
                {
                    "name": "Glycerophospholipids",
                    "position_x": -480,
                    "position_y": 0,
                    "position_z": 0,
                    "last_layer_scale": 1,
                    "rotation_x": 0,
                    "rotation_y": 0,
                    "rotation_z": 0,
                    "floor_current_color": "#777777",
                    "geometry_parameters_width": 2471
                },
                {
                    "name": "Enzymes",
                    "position_x": 0,
                    "position_y": 0,
                    "position_z": 0,
                    "last_layer_scale": 1,
                    "rotation_x": 0,
                    "rotation_y": 0,
                    "rotation_z": 0,
                    "floor_current_color": "#777777",
                    "geometry_parameters_width": 2471
                },
                {
                    "name": "Sphingolipids",
                    "position_x": 480,
                    "position_y": 0,
                    "position_z": 0,
                    "last_layer_scale": 1,
                    "rotation_x": 0,
                    "rotation_y": 0,
                    "rotation_z": 0,
                    "floor_current_color": "#777777",
                    "geometry_parameters_width": 2471
                }
            ],
            "universalLabelColor": "#ffffff",
            "direction": True,
            "edgeOpacityByWeight": True
        }

    def _clear_caches(self) -> None:
        """Clear all caches in the data access layer"""
        self.data_access.clear_caches()

    def _preprocess_and_translate(self, lipids: List[str], translate: bool = False) -> Tuple[
        PreprocessingResult, Optional[Dict[str, List[str]]]]:
        """
        Preprocess and optionally translate the lipids.
        """
        # Clean lipids
        cleaned_lipids = [self.lipid_translator.clean_lipid_abbreviation(lipid) for lipid in lipids]

        unique_cleaned_set = set(cleaned_lipids)

        processed_lipids = list(unique_cleaned_set)

        # Translate
        translation_tracker = self.lipid_translator.translate_lipids(processed_lipids)
        self.translation_mapper = TranslationMapper(translation_tracker)
        processed_lipids = list(set(name for names in translation_tracker.values() for name in names))

        # Use the existing preprocessor
        preprocessing_result = self.lipid_preprocessor.preprocess_lipids(processed_lipids)

        return preprocessing_result, translation_tracker

    def _prepare_results(self, reaction_results: Dict, preprocessing_result: PreprocessingResult,
                         translation_tracker: Optional[Dict[str, List[str]]],
                         original_lipids: Set[str]) -> Dict:  # <-- Add new parameter
        """
        Prepare the final results dictionary from processing outputs.
        MODIFIED: Now includes the original, untranslated lipid set for output generation.
        """
        results = reaction_results.copy()

        # This is the set of all translated names, used for internal processing
        results['processed_lipids'] = preprocessing_result.all_valid_lipids

        # --- KEY FIX: Store the original lipid names for output ---
        results['original_input_lipids'] = original_lipids
        # ---------------------------------------------------------

        results['unparseable_lipids'] = preprocessing_result.unparseable_lipids

        if translation_tracker:
            results['translation_tracker'] = translation_tracker

        return results

    def process_files(self, file_paths: List[str], organism_id: int, translate: bool, process_network: bool) -> None:
        """
        Process multiple files containing lipid data.

        Args:
            file_paths: List of paths to files containing lipid data
            organism_id: ID of the organism to process
            translate: Whether to translate lipids before processing
            process_network: Whether to process and save lipid network

        Returns:
            None

        Note:
            Results structure will contain:
            - reaction_dict: Dict[Tuple, ReactionValue]
            - super_reaction_dict: Dict[Tuple, ReactionValue]
            - enzyme_dict: Dict[str, EnzymeInfo]
            - input_lipids: Set[str]
            - unparseable_lipids: List[str]
            - translation_tracker (optional): Dict[str, List[str]]
        """
        """Process multiple files containing lipid data."""
        try:
            for file_path in file_paths:
                # Extract lipids from file
                raw_lipids = self._extract_lipids_from_file(file_path)

                # Clean the raw lipids first for consistency.
                original_lipids_set = {self.lipid_translator.clean_lipid_abbreviation(lipid) for lipid in raw_lipids}

                # Preprocess and translate lipids
                preprocessing_result, translation_tracker = self._preprocess_and_translate(raw_lipids, translate)

                # Update both lipid_matcher and reaction_processor with headgroup map
                self.lipid_matcher.set_headgroup_map(preprocessing_result.headgroup_map)
                self.reaction_processor.set_headgroup_map(preprocessing_result.headgroup_map)

                # Update translation mapper in lipid processor
                self.lipid_processor.set_translation_mapper(self.translation_mapper)

                # Process lipids
                reaction_results = self.process_lipid_list(
                    sorted(preprocessing_result.all_valid_lipids),
                    organism_id
                )

                # Prepare final results
                results = self._prepare_results(reaction_results, preprocessing_result, translation_tracker, original_lipids_set)

                # Generate outputs
                base_name = os.path.splitext(os.path.basename(file_path))[0]
                self.generate_outputs(results, translate, base_name, preprocessing_result.headgroup_map)

                if process_network:
                    self._process_lipid_network(results, base_name)

            logger.info("File processing completed successfully")
        except Exception as error:
            logger.error(f"File processing failed: {str(error)}")
            raise

    def _extract_lipids_from_file(self, file_path: str) -> List[str]:
        """
        Extract lipids from an Excel file.

        Args:
            file_path (str): Path to the Excel file.

        Returns:
            List[str]: List of extracted lipids.
        """
        df = pd.read_excel(file_path, header=None)
        return [lipid.strip() for lipid in df.iloc[0] if isinstance(lipid, str) and lipid.strip()]

    def _process_lipids(self, lipids: List[str]) -> Tuple[List[str], Dict[str, List[str]]]:
        """
        Process lipids, with optimal translation and cleaning

        Args:
            lipids (List[str]): Raw input lipids.
            translate (bool): Whether to translate the lipids.

        Returns:
            List[str]: Processed lipids.
        """
        # First cleaning pass on input lipids
        initial_cleaned_lipids = [self.lipid_translator.clean_lipid_abbreviation(lipid) for lipid in lipids]

        # Convert to set immediately for O(1) lookups
        unique_cleaned_set = set(initial_cleaned_lipids)

        # Translation is now mandatory
        # translate_lipids will clean the translated lipids internally
        translation_tracker = self.lipid_translator.translate_lipids(list(unique_cleaned_set))
        self.translation_mapper = TranslationMapper(translation_tracker)


        # Get all translated versions (already cleaned)
        final_lipids = set(
            name for names in translation_tracker.values()
            for name in names
        )

        return list(final_lipids), translation_tracker

    def run(self, lipids: List[str], organism_id: int, translate: bool = False) -> None:
        try:
            if translate:
                translation_tracker = self.lipid_translator.translate_lipids(lipids)
                translated_lipids = list(set([name for names in translation_tracker.values() for name in names]))
                results = self.process_lipid_list(translated_lipids, organism_id)
                results['translation_tracker'] = translation_tracker
                results['input_lipids'] = set(translated_lipids)
            else:
                results = self.process_lipid_list(lipids, organism_id)
                results['input_lipids'] = set(lipids)

            self.generate_outputs(results, translate, "test_run")
            logger.info("Lipid analysis completed successfully")
        except Exception as e:
            logger.error(f"Lipid analysis failed: {str(e)}")

# Example usage
def main():
    workflow = LipidAnalysisWorkflow('lipid_config.yaml', '../../')
    lipids = ["SM(d18:1-D3/14:0) +2", "Cer(d18:1/16:0", "PC(16:0/18:1)", "SM(d18:1/16:0)", "PC(24:0/20:0)"]
    organism_id = 9606

    workflow.run(lipids, organism_id, translate=False)

if __name__ == "__main__":
    main()




# TODO: handle Hex lipids in input lipids
