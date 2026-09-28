from typing import NamedTuple, List, Optional, Dict, Set, Union
from enum import Enum
from dataclasses import dataclass, field


class ReactionType(Enum):
    REGULAR = "REGULAR"
    SUPER = "SUPER"


@dataclass
class EnzymeInfo:
    enzyme_name: str
    uniprot_id: str
    subcellular_locations: Optional[str] = None

class MatchResult(Enum):
    """
    Enum class that is used when trying to match SUPER reactions, used in 2 cases:
    1. (NO_MATCH) Used when reactant and product do not share chain length info, need to grab 2nd reactant for matching
    2. (SUPER_MATCH) Used when matching requires components from both reactants, thus triggers super_match_sort
    """
    NO_MATCH = 1
    SUPER_MATCH = 2

class ChainPosition(Enum):
    SN1 = 'sn1'
    SN2 = 'sn2'

    def __str__(self):
        return self.value


class Reaction(NamedTuple):
    reaction_id: int
    reaction_enzyme_id: int
    enzyme_name: str
    uniprot_id: str
    reactants: List[str]  # List of reactant(s)
    products: List[str]  # Product abbreviation
    rhea_id: Union[str, int]
    doi: Optional[str]
    reaction_type: ReactionType
    subcellular_locations: Optional[str] = None
    generalized_from_ids: Optional[str] = None

    def get_key(self):
        sorted_reactants = tuple(sorted(self.reactants))
        sorted_products = tuple(sorted(self.products))
        return (sorted_reactants, sorted_products, self.enzyme_name, self.doi)

    @property
    def key_with_rhea(self):
        return self.get_key() + (self.rhea_id,)


class ReactionValue(NamedTuple):
    rhea_id: Optional[Union[int, str]]
    reaction: Reaction

    def update_rhea_id(self, new_rhea_id: Optional[Union[int, str]]) -> 'ReactionValue':
        rhea_ids_set = set()

        # Helper function to safely add rhea_id to the set
        def add_to_set(id_value):
            if id_value is not None:
                if isinstance(id_value, int):
                    rhea_ids_set.add(str(id_value))
                elif isinstance(id_value, str):
                    rhea_ids_set.update(id_value.split(';'))

        # Add existing rhea_id(s) to the set
        add_to_set(self.rhea_id)

        # Add new rhea_id to the set
        add_to_set(new_rhea_id)

        # Join unique ids into a string
        updated_rhea_id = ';'.join(sorted(rhea_ids_set)) if rhea_ids_set else None

        return ReactionValue(updated_rhea_id, self.reaction)

    def merge_with(self, other_value: 'ReactionValue') -> 'ReactionValue':
        """
        NEW, ROBUST MERGE LOGIC: Merges this ReactionValue with another,
        intelligently combining Rhea IDs and `generalized_from_ids`.
        """
        # --- 1. Combine Rhea IDs ---
        # This logic correctly handles merging multiple IDs.
        rhea_ids_set = set()

        def add_rhea_to_set(id_value):
            if id_value is not None and str(id_value).strip() != '':
                # Handle semicolon-separated strings of Rhea IDs
                if isinstance(id_value, str) and ';' in id_value:
                    rhea_ids_set.update(filter(None, id_value.split(';')))
                else:
                    rhea_ids_set.add(str(id_value))

        add_rhea_to_set(self.rhea_id)
        add_rhea_to_set(self.reaction.rhea_id)
        add_rhea_to_set(other_value.rhea_id)
        add_rhea_to_set(other_value.reaction.rhea_id)

        updated_rhea_id = ';'.join(sorted(list(rhea_ids_set))) if rhea_ids_set else None

        # --- 2. Combine `generalized_from_ids` (this logic was already correct) ---
        generalized_ids_set = set()

        def add_gen_to_set(id_value):
            if id_value:
                generalized_ids_set.update(filter(None, id_value.strip('|').split('|')))

        add_gen_to_set(self.reaction.generalized_from_ids)
        add_gen_to_set(other_value.reaction.generalized_from_ids)

        updated_generalized_ids = f"|{'|'.join(sorted(list(generalized_ids_set)))}|" if generalized_ids_set else None

        # --- 3. Return the new, correctly merged ReactionValue ---
        # The new Reaction object should also contain the merged Rhea ID for consistency.
        updated_reaction = self.reaction._replace(
            generalized_from_ids=updated_generalized_ids,
            rhea_id=updated_rhea_id  # Also update the nested Rhea ID
        )

        return ReactionValue(rhea_id=updated_rhea_id, reaction=updated_reaction)


@dataclass
class OutputContext:
    reaction_dict: Dict
    super_reaction_dict: Dict
    input_lipids: Set[str]
    enzyme_dict: Dict
    headgroup_map: Dict[str, Set[str]]
    translation_tracker: Optional[Dict[str, list[str]]] = None
    specific_rhea_map: Dict[int, str] = field(default_factory=dict)
    base_name: str = ""
    output_dir: str = ""
    debug_mode: bool = False




@dataclass
class SingleLipidSearchContext:
    """Context object for single lipid search results."""
    # This structure mirrors the output of LipidAnalysisWorkflow.search_products_by_lipid
    search_results: Dict[str, Dict]
    base_name: str
    output_dir: str
    debug_mode: bool = False
    specific_rhea_map: Dict[int, str] = field(default_factory=dict)


@dataclass
class SingleEnzymeSearchContext:
    """Context object for single enzyme search results."""
    # This structure mirrors the output of LipidAnalysisWorkflow.search_reactions_by_enzyme
    search_results: Dict[str, Union[str, int, bool, List]]
    base_name: str
    output_dir: str
    debug_mode: bool = False