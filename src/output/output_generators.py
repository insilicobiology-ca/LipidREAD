from typing import Dict, Tuple, Set, List, Optional, Protocol, Any
import numpy as np
import pandas as pd
import os
from src.models.models import (Reaction, ReactionValue, OutputContext, ReactionType, SingleLipidSearchContext,
                               SingleEnzymeSearchContext)
from src.models.lipid_components import SphingoLipidComponents, SphingoidBaseComponents
from src.data.config_loader import ConfigLoader
from src.parsing.lipid_parser import LipidParser, LipidParserFactory, LipidComponentCache
from src.matching.lipid_matcher import LipidMatcher
from collections import defaultdict
from itertools import product

import logging
logger = logging.getLogger(__name__)


class OutputGenerator(Protocol):
    def generate(self, context: Any):
        ...


def create_output_generators(config: ConfigLoader, debug_mode: bool = False) -> List[OutputGenerator]:
    output_paths = config.get_output_paths()

    # NEW: Instantiate parser and matcher here
    lipid_component_cache = LipidComponentCache()
    lipid_parser_factory = LipidParserFactory(config)
    lipid_parser = LipidParser(
        factory=lipid_parser_factory,
        components_cache=lipid_component_cache
    )
    lipid_matcher = LipidMatcher(lipid_parser, config)

    return [
        CombinedOutputGenerator(
            matrix_file=output_paths['adjacency_matrix'],
            binary_matrix_file=output_paths['binary_matrix'],
            reaction_list_file=output_paths['reaction_list'],
            config=config,
            lipid_parser=lipid_parser,  # Pass the instance
            lipid_matcher=lipid_matcher,  # Pass the instance
            debug_mode=debug_mode
        ),
        MultiReactionListGenerator(
            output_file=output_paths.get('full_reaction_list', 'full_reaction_list.tsv')
        ),
        EnzymeOutputGenerator(output_paths['enzyme_table']),
        TranslateOutputGenerator(output_paths['translation_doc'])
    ]


class SingleLipidSearchOutputGenerator:
    """
    Generates a single TSV report for the results of one or more
    independent single lipid searches, grouped by the input lipid.
    """

    def __init__(self, output_file: str):
        self.output_file = output_file

    def generate(self, context: 'SingleLipidSearchContext') -> None:
        """
        Creates a TSV file from the unified single lipid analysis context,
        using a 'Role' column to distinguish forward and reverse reactions.
        """
        output_data = []

        for searched_lipid, analysis in context.search_results.items():
            if 'error' in analysis:
                # Add a placeholder for Role in the error row for consistency
                output_data.append([searched_lipid, 'ERROR', analysis['error'], '', '', '', '', '', ''])
                continue

            # Process forward reactions ("Produces")
            for reaction in analysis.get('forward_reactions', []):
                output_data.append(self._format_reaction_row(
                    searched_lipid=searched_lipid,
                    role="Reactant",  # Clearer role name
                    reaction=reaction,
                    context=context
                ))

            # Process reverse reactions ("Is Produced By")
            for reaction in analysis.get('reverse_reactions', []):
                output_data.append(self._format_reaction_row(
                    searched_lipid=searched_lipid,
                    role="Product",  # Clearer role name
                    reaction=reaction,
                    context=context
                ))

            if not analysis.get('forward_reactions') and not analysis.get('reverse_reactions'):
                output_data.append([searched_lipid, 'No connections found', '', '', '', '', '', '', ''])

        header = [
            "Searched Lipid", "Role", "Reactant(s)", "Product(s)",
            "Enzyme", "Rhea ID", "DOI", "Sourced from # of Rhea IDs", "Sourced from Rhea Example"
        ]
        if context.debug_mode:
            header.append("Reaction ID")

        df = pd.DataFrame(output_data, columns=header)
        # Sort the final output to group by the searched lipid and then by its role
        df.sort_values(by=["Searched Lipid", "Role"], inplace=True)

        output_path = os.path.join(context.output_dir, f"{context.base_name}_{self.output_file}")
        df.to_csv(output_path, index=False, sep='\t')
        logger.info(f"Single lipid analysis report saved to {output_path}")

    def _format_reaction_row(self, searched_lipid: str, role: str, reaction: Reaction,
                             context: 'SingleLipidSearchContext') -> List[str]:
        """Helper to format a single row for the output report."""
        full_reactants_str = " + ".join(sorted(reaction.reactants))
        full_products_str = " + ".join(sorted(reaction.products))

        num_generalized, specific_rhea = '', ''
        generalized_ids_str = reaction.generalized_from_ids

        # --- ROBUSTNESS FIX for Specific Rhea ID ---
        if (generalized_ids_str):
            # Filter out empty strings resulting from splitting "|123|456|"
            id_list = [i for i in generalized_ids_str.strip('|').split('|') if i.isdigit()]

            if id_list:
                num_generalized = str(len(id_list))
                # Now id_list[0] is guaranteed to be a digit string
                first_specific_id = int(id_list[0])
                specific_rhea = str(context.specific_rhea_map.get(first_specific_id, ''))
        # -------------------------------------------

        row = [
            searched_lipid, role, full_reactants_str, full_products_str,
            reaction.enzyme_name, str(reaction.rhea_id) if reaction.rhea_id else '',
            reaction.doi or '', num_generalized, specific_rhea
        ]
        if context.debug_mode:
            row.append(str(reaction.reaction_id))

        return row

class SingleEnzymeSearchOutputGenerator:
    """
    Generates a simple TSV report for a single enzyme search.
    """
    def __init__(self, output_file: str):
        self.output_file = output_file

    def generate(self, context: 'SingleEnzymeSearchContext') -> None:
        results = context.search_results
        output_data = []

        if not results.get('found'):
            output_data.append([results.get('enzyme_identifier'), 'Enzyme not found', '', '', ''])
        else:
            # The 'reactions' key now holds our list of aggregated dictionaries.
            for agg_reaction in results['reactions']:
                output_data.append([
                    results['enzyme_name'],
                    results['uniprot_id'],
                    agg_reaction.get('abbreviated_reaction_text', ''),
                    agg_reaction.get('rhea_id_example', ''),
                    agg_reaction.get('consolidation_count', 0),
                ])

        # --- UPDATE THE HEADER ---
        header = [
            "Enzyme Name", "UniProt ID", "Abbreviated Reaction",
            "Sourced from Rhea Example", "Sourced from # Reactions"
        ]
        df = pd.DataFrame(output_data, columns=header)

        # Sort by the consolidation count to show most-evidenced reactions first
        if 'Sourced from # Reactions' in df.columns:
            df.sort_values(by='Sourced from # Reactions', ascending=False, inplace=True)

        output_path = os.path.join(context.output_dir, f"{context.base_name}_{self.output_file}")
        df.to_csv(output_path, index=False, sep='\t')
        logger.info(f"Single enzyme search report saved to {output_path}")


class MultiReactionListGenerator:
    """
    NEW: Generates a detailed list of all multi-participant reactions.
    Format: Reactant1 | Reactant2 | ... | Product1 | Product2 | ... | Enzyme | ...
    """

    def __init__(self, output_file: str):
        self.output_file = output_file

    def generate(self, context: OutputContext) -> None:
        all_reactions = list(context.reaction_dict.values()) + list(context.super_reaction_dict.values())
        if not all_reactions:
            return

        # Determine the maximum number of reactants and products across all reactions
        max_reactants = max(len(rv.reaction.reactants) for rv in all_reactions)
        max_products = max(len(rv.reaction.products) for rv in all_reactions)

        specific_rhea_map = context.specific_rhea_map

        reactant_cols = [f"Reactant{i + 1}" for i in range(max_reactants)]
        product_cols = [f"Product{i + 1}" for i in range(max_products)]
        other_cols = ["Enzyme", "Rhea_ID", "Generalized_from_#", "Specific_Rhea_ID", "DOI"]
        if context.debug_mode:
            other_cols.append("Reaction_ID")

        header = reactant_cols + product_cols + other_cols

        output_data = []
        for rv in all_reactions:
            reaction = rv.reaction
            row = []

            # NEW: Calculate the values for the new columns
            generalized_ids = reaction.generalized_from_ids
            num_generalized = 0
            specific_rhea = ''
            if generalized_ids:
                # Count the number of IDs by splitting the string
                id_list = generalized_ids.strip('|').split('|')
                num_generalized = len(id_list)
                first_specific_id = int(id_list[0])
                specific_rhea = str(specific_rhea_map.get(first_specific_id, ''))

            # Pad reactants and products to the max length
            row.extend(list(reaction.reactants) + [''] * (max_reactants - len(reaction.reactants)))
            row.extend(list(reaction.products) + [''] * (max_products - len(reaction.products)))

            row.extend([
                reaction.enzyme_name,
                str(rv.rhea_id) if rv.rhea_id else None,
                str(num_generalized) if num_generalized > 0 else '',  # New data
                specific_rhea,
                reaction.doi or ''
            ])
            if context.debug_mode:
                row.append(str(reaction.reaction_id))
            #logger.info(f"SPOTA reactions at this point {reaction}")
            output_data.append(row)

        df = pd.DataFrame(output_data, columns=header)
        df.sort_values(by=reactant_cols + product_cols, inplace=True)
        output_path = os.path.join(context.output_dir, f"{context.base_name}_{self.output_file}")
        df.to_csv(output_path, index=False, sep='\t')


class TranslateOutputGenerator:
    def __init__(self, output_file: str):
        self.translate_output_file = output_file

    def generate(self, context: OutputContext) -> None:
        if context.translation_tracker is None:
            return

        output_list = []
        for original_name, translated_names in context.translation_tracker.items():
            concatenated_names = ';'.join(translated_names)
            output_list.append([original_name, concatenated_names])

        df = pd.DataFrame(output_list, columns=['Original name', 'Translated name'])
        output_path = os.path.join(context.output_dir, f"{context.base_name}_{self.translate_output_file}")
        df.to_csv(output_path, index=None, sep='\t')


class EnzymeOutputGenerator:
    def __init__(self, output_file: str):
        self.enzyme_output_file = output_file

    def generate(self, context: OutputContext) -> None:
        df = pd.DataFrame([(name, info.uniprot_id, info.subcellular_locations) for name, info in context.enzyme_dict.items()],
                          columns=['Gene Name', 'UniProt AC', 'Subcellular Location(s)'])
        output_path = os.path.join(context.output_dir, f"{context.base_name}_{self.enzyme_output_file}")
        df.to_csv(output_path, index=False)


class CombinedOutputGenerator:
    """
    MODIFIED: This generator now produces the simplified Reactant-Product PAIR list
    and the corresponding adjacency matrices.
    """
    def __init__(self, matrix_file: str, binary_matrix_file: str, reaction_list_file: str,
                 config: ConfigLoader, lipid_parser: LipidParser, lipid_matcher: LipidMatcher,
                 debug_mode: bool = False):
        self.matrix_file = matrix_file
        self.binary_matrix_file = binary_matrix_file
        self.pair_list_file = reaction_list_file
        self.config = config
        self.lipid_parser = lipid_parser
        self.lipid_matcher = lipid_matcher
        self.debug_mode = debug_mode
        self.manual_connections = self.config.get('manual_connections', [])

    def generate(self, context: OutputContext) -> None:
        matrix, _, molecule_dict = self._initialize_matrices(context.input_lipids)

        pair_list = self._generate_pair_list_from_full_reactions(context)

        if matrix is not None:
            for pair_data in pair_list:
                reactant, product, enzyme = pair_data[0], pair_data[1], pair_data[2]
                if (reactant in molecule_dict and product in molecule_dict) :
                    self._add_to_matrix(matrix, molecule_dict, reactant, product, enzyme)

        self._save_outputs(matrix, pair_list, context)

    def _generate_pair_list_from_full_reactions(self, context: OutputContext) -> List[List[str]]:
        """
        Generates a preliminary list of all possible valid pairs, then de-duplicates
        this list using a tie-breaking rule that prefers the most specific reaction context.
        """
        # --- Step 1: Generate ALL possible valid pairs from all perspectives ---
        preliminary_pair_list = []
        all_reaction_values = list(context.reaction_dict.values()) + list(context.super_reaction_dict.values())

        specific_rhea_map = context.specific_rhea_map
        #logger.info(f"SPOTB reactions at this point {all_reaction_values}")
        for rv in all_reaction_values:
            reaction = rv.reaction
            #GOODlogger.info(f"SPOTE reactions at this point {reaction}")
            smoshed_reaction_text = self._build_smoshed_text(reaction)
            class_level_text = self._get_class_level_reaction_text(reaction)
            enzyme_display = reaction.enzyme_name
            #GOODlogger.info(f"SPOTF reactions at this point {smoshed_reaction_text}")

            # Calculate the specific Rhea ID for this generalized reaction
            generalized_ids_str = reaction.generalized_from_ids
            num_generalized = 0
            specific_rhea = ''
            if generalized_ids_str:
                id_list = [i for i in generalized_ids_str.strip('|').split('|') if i]
                if id_list:
                    num_generalized = len(id_list)
                    first_specific_id = int(id_list[0])
                    specific_rhea = str(specific_rhea_map.get(first_specific_id, ''))

            specific_reactants = {r for r in reaction.reactants if not self._is_class(r)}
            specific_products = {p for p in reaction.products if not self._is_class(p)}

            from itertools import product
            # A. Generate pairs from DIRECT matches and SUPER matches (non-manual)
            for r_specific, p_specific in product(specific_reactants, specific_products):
                #logger.info(f"SPOTH reactions at this point {r_specific,p_specific}")
                #logger.info(f"SPOTK reactions at this point {self._is_allowed_pair(r_specific, p_specific)}")
                if (self._is_allowed_pair(r_specific, p_specific) and not self._is_manual_connection(r_specific,
                                                                                                     p_specific)) or True:
                    #MCCissue is with self._is_allowed_pair - made this if statement True every time with or True addition
                    #this if statement changes pair list and adjacency matrix from the full list
                    pair_row = self._format_pair_row(
                        r_specific, p_specific, enzyme_display, class_level_text,
                        smoshed_reaction_text, reaction, num_generalized, specific_rhea
                    )
                    preliminary_pair_list.append(pair_row)
            #GOOD IF is_allowed removed logger.info(f"SPOTG reactions at this point {preliminary_pair_list}")

            # B. Handle manual cross-connections with inference
            for connection in self.manual_connections:
                try:
                    r_conn_class, p_conn_class = connection.split('_')
                except ValueError:
                    continue

                # Check if this manual connection is relevant to the current reaction event's classes
                event_reactant_classes = {self._get_class(r) for r in reaction.reactants}
                event_product_classes = {self._get_class(p) for p in reaction.products}

                if (r_conn_class in event_reactant_classes and p_conn_class in event_product_classes):
                    # This rule applies. Generate all possible specific pairs from the user's total input.
                    all_possible_reactants = context.headgroup_map.get(r_conn_class, set())
                    all_possible_products = context.headgroup_map.get(p_conn_class, set())

                    for r_specific_manual, p_specific_manual in product(all_possible_reactants, all_possible_products):
                        # --- INFERENCE LOGIC ---
                        # For this specific manual pair, infer what the co-participants would be.
                        inferred_reactants, inferred_products = self._infer_partners(r_specific_manual,
                                                                                     p_specific_manual, reaction)
                        inferred_reaction_text = self._build_inferred_smoshed_text(inferred_reactants,
                                                                                   inferred_products)

                        pair_row = self._format_pair_row(r_specific_manual, p_specific_manual, enzyme_display,
                                                         class_level_text,
                                                         inferred_reaction_text, reaction, num_generalized,
                                                         specific_rhea)
                        preliminary_pair_list.append(pair_row)

        # --- Step 2: De-duplication using the final Tie-Breaking Rule ---
        final_pairs = {}
        # The key identifies a unique pair: (Reactant, Product, Enzyme, Class_Level_Reaction)
        for row in preliminary_pair_list:
            key = (row[0], row[1], row[2], row[3])  # The core identity of the pair
            #logger.info(f"SPOTD reactions at this point {key}")

            # The Tie-Breaking Rule:
            # Prefer the row where the 'Product' of the pair (col 1) appears as a
            # specific molecule in the product side of the 'Full_Reaction' string (col 4).
            product_in_pair = row[1]
            # Use index 4 for 'Full_Reaction' column
            full_reaction_products_part = row[4].split(' => ')[1]

            is_primary_product_perspective = product_in_pair in full_reaction_products_part

            if (key not in final_pairs or is_primary_product_perspective):
                final_pairs[key] = row

            #logger.info(f"SPOTC reactions at this point {final_pairs}")

        return sorted(list(final_pairs.values()))

    def _build_inferred_smoshed_text(self, reactants: Set[str], products: Set[str]) -> str:
        """Builds the reaction string from inferred sets."""
        reactants_str = " + ".join(sorted(list(reactants)))
        products_str = " + ".join(sorted(list(products)))
        return f"{reactants_str} => {products_str}"

    def _is_manual_connection(self, reactant: str, product: str) -> bool:
        """Checks if a pair ONLY matches a manual connection rule."""
        r_class = self._get_class(reactant)
        p_class = self._get_class(product)
        is_manual = f"{r_class}_{p_class}" in self.manual_connections
        # A pair is ONLY manual if it's in the manual list AND not a direct biochemical match
        is_direct = self._is_allowed_pair(reactant, product)
        return is_manual and not is_direct

    def _infer_partners(self, r_specific: str, p_specific: str, reaction: Reaction) -> Tuple[Set[str], Set[str]]:
        """
        For a given manual pair (r_specific, p_specific), infers the other
        reactants and products for the Full_Reaction string.
        """
        r_class = self._get_class(r_specific)
        p_class = self._get_class(p_specific)

        # Start with the known specific participants
        final_reactants = {r_specific}
        final_products = {p_specific}

        all_reactant_classes = {self._get_class(r) for r in reaction.reactants}
        all_product_classes = {self._get_class(p) for p in reaction.products}

        # Find the "other" classes
        other_reactant_class = (all_reactant_classes - {r_class}).pop()
        other_product_class = (all_product_classes - {p_class}).pop()

        # Infer the partners by reverse-matching
        inferred_reactant = self.lipid_matcher.match_sort(p_specific, other_reactant_class, set())
        inferred_product = self.lipid_matcher.match_sort(r_specific, other_product_class, set())

        if isinstance(inferred_reactant, str):
            final_reactants.add(inferred_reactant)
        else:
            final_reactants.add(other_reactant_class)  # Fallback to class name

        if isinstance(inferred_product, str):
            final_products.add(inferred_product)
        else:
            final_products.add(other_product_class)  # Fallback to class name

        return final_reactants, final_products

    def _get_class(self, lipid_name: str) -> str:
        """Helper to robustly get the headgroup class name."""
        return self.lipid_parser.parse_lipid(lipid_name).headgroup

    def _build_smoshed_text(self, reaction: Reaction) -> str:
        """
        Helper to build the specific molecular reaction string for the report.
        This represents the single molecular event.
        """
        reactants_str = " + ".join(sorted(reaction.reactants))
        products_str = " + ".join(sorted(reaction.products))
        return f"{reactants_str} => {products_str}"

    def _is_allowed_pair(self, reactant: str, product: str) -> bool:
        """
        Checks if a reactant-product pair is allowed by either a direct
        matching rule, a super-matching rule, or a manual connection.
        """
        r_class = self.lipid_parser.parse_lipid(reactant).headgroup
        p_class = self.lipid_parser.parse_lipid(product).headgroup
        #logger.info(f"SPOTI reactions at this point {r_class, p_class}")

        # Use the new, unified checking method from LipidMatcher
        if self.lipid_matcher.is_valid_pair(r_class, p_class):
            return True

        # Also check manual connections
        if f"{r_class}_{p_class}" in self.manual_connections:
            return True

        return False

    def _is_class(self, lipid_name: str) -> bool:
        """Determines if a lipid name represents a class or a specific molecule."""
        return self.lipid_parser.parse_lipid(lipid_name).is_headgroup_only

    def _get_class_level_reaction_text(self, reaction: Reaction) -> str:
        """Generates a class-level representation of the reaction, e.g., 'Cer + PC => SM + DG'."""
        r_classes = " + ".join(sorted([self.lipid_parser.parse_lipid(r).headgroup for r in reaction.reactants]))
        p_classes = " + ".join(sorted([self.lipid_parser.parse_lipid(p).headgroup for p in reaction.products]))
        return f"{r_classes} => {p_classes}"

    def _format_pair(self, reactant: str, product: str, reaction: Reaction, class_level_text: str) -> List[str]:
        """Formats a single row for the pair list."""
        # For second-order reactions, the enzyme name is augmented
        enzyme_display = reaction.enzyme_name
        if len(reaction.reactants) > 1:
            # Find the "other" reactant that is not the one in this pair
            other_reactant = next((r for r in reaction.reactants if r != reactant), "")
            enzyme_display = f"{reaction.enzyme_name}${other_reactant}"

        pair = [
            reactant,
            product,
            enzyme_display,
            class_level_text,
            str(reaction.rhea_id) if reaction.rhea_id else None,
            reaction.doi
        ]
        if self.debug_mode:
            pair.append(str(reaction.reaction_id))
        return pair

    def _format_pair_row(self, reactant: str, product: str, enzyme_display: str,
                         class_level_text: str, smoshed_reaction_text: str,
                         reaction: Reaction, num_generalized: int, specific_rhea: str) -> List[str]:
        """
        Formats a single row for the pair list, including the new generalization columns.
        """
        pair = [
            reactant,
            product,
            enzyme_display,
            class_level_text,
            smoshed_reaction_text,
            str(reaction.rhea_id) if reaction.rhea_id else '',
            str(num_generalized) if num_generalized > 0 else '',
            specific_rhea,
            reaction.doi or ''
        ]
        if self.debug_mode:
            pair.append(str(reaction.reaction_id))
        return pair

    # ... (the _initialize_matrices, _add_to_matrix, and _save_outputs methods are refactored slightly
    # to work with the new flow, but their core logic remains) ...

    def _initialize_matrices(self, input_lipids: Set[str]) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Dict[str, int]]:
        sorted_lipids = sorted(list(input_lipids))
        molecule_dict = {lipid: i + 1 for i, lipid in enumerate(sorted_lipids)}
        if len(input_lipids) > 2000: return None, None, molecule_dict
        matrix_size = len(input_lipids) + 1
        matrix = np.empty((matrix_size, matrix_size), dtype='<U257')
        matrix.fill('')
        for i, lipid in enumerate(sorted_lipids, 1):
            matrix[0, i] = matrix[i, 0] = lipid
        return matrix, matrix.copy(), molecule_dict

    def _add_to_matrix(self, matrix: np.ndarray, molecule_dict: Dict[str, int], reactant: str, product: str, enzyme: str) -> None:
        x, y = molecule_dict[reactant], molecule_dict[product]
        if not matrix[x, y]:
            matrix[x, y] = enzyme
        else:
            enzyme_set = set(matrix[x, y].split(";"))
            enzyme_set.add(enzyme)
            matrix[x, y] = ";".join(sorted(list(enzyme_set)))

    def _save_outputs(self, matrix: Optional[np.ndarray], pair_list: List[List[str]], context: OutputContext) -> None:
        columns = [
            "Reactant", "Product", "Enzyme", "Class_Level_Reaction",
            "Full_Reaction", "Rhea_ID", "Sourced from # of Rhea IDs",
            "Sourced from Rhea Example", "DOI"
        ]
        if self.debug_mode: columns.append("Reaction_ID")

        df = pd.DataFrame(pair_list, columns=columns)
        df.to_csv(os.path.join(context.output_dir, f"{context.base_name}_{self.pair_list_file}"), index=False, sep='\t')

        if matrix is not None:
            pd.DataFrame(matrix).to_csv(os.path.join(context.output_dir, f"{context.base_name}_{self.matrix_file}"),
                                        header=None, index=None)
            binary_matrix = matrix.copy()
            binary_matrix[1:, 1:] = np.where(binary_matrix[1:, 1:] != '', '1', '0')
            pd.DataFrame(binary_matrix).to_csv(
                os.path.join(context.output_dir, f"{context.base_name}_{self.binary_matrix_file}"), header=None,
                index=None)
