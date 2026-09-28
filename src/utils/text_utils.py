import re
import logging
from typing import List, Tuple, Optional
import itertools # Import the itertools module

# Get a logger for this module
logger = logging.getLogger(__name__)


def split_reaction_text(reaction_text: str) -> Tuple[Optional[List[List[str]]], bool]:
    """
    Splits a full reaction text into reactant and product name lists.
    Handles various arrow types and returns a reversibility flag.
    This is a standalone utility function.
    """
    if not isinstance(reaction_text, str):
        logger.warning(
            f"Invalid input to split_reaction_text: expected a string, got {type(reaction_text)}"
        )
        return (None, False)

    reactants_str, products_str = "", ""
    reversible = False

    arrow_map = {" <=> ": True, " => ": False, " <= ": False, " + + ": True}

    found_arrow = None
    for arrow in arrow_map:
        if arrow in reaction_text:
            found_arrow = arrow
            reversible = arrow_map[arrow]
            parts = reaction_text.split(arrow, 1)
            if len(parts) != 2:  # Safety check
                continue
            if found_arrow == " <= ":
                products_str, reactants_str = parts
            else:
                reactants_str, products_str = parts
            break

    if not found_arrow:
        logger.warning(
            f"Could not find a valid arrow in reaction text: '{reaction_text}'"
        )
        return (None, False)

    def clean_participant(p_str: str) -> str:
        cleaned = re.sub(r"^\d+\s+", "", p_str.strip())
        if cleaned.lower().startswith("a "):
            cleaned = cleaned[2:]
        return cleaned

    reactants = [clean_participant(r) for r in reactants_str.split(" + ") if r.strip()]
    products = [clean_participant(p) for p in products_str.split(" + ") if p.strip()]

    if not reactants or not products:
        logger.warning(
            f"Could not parse valid reactants or products from '{reaction_text}'"
        )
        return (None, False)

    return ([reactants, products], reversible)

def generate_reaction_combinations(reactant_ids: List[int],
                                   product_ids: List[int],
                                   reversible: bool
                                   ) -> List[Tuple[int, int]]:
    """
    Generate all possible combinations of (reactant_id, product_id) pairs.
    If reversible, also generates (product_id, reactant_id) pairs.

    Args:
        reactant_ids (List[int]): List of reactant molecule IDs.
        product_ids (List[int]): List of product molecule IDs.
        reversible (bool): Flag indicating if the reaction is reversible.

    Returns:
        List[Tuple[int, int]]: A list of (reactant_id, product_id) tuples.
    """
    if not reactant_ids or not product_ids:
        return []

    # Generate forward pairs: (each reactant, each product)
    forward_combinations = list(itertools.product(reactant_ids, product_ids))

    all_combinations = forward_combinations

    if reversible:
        # Generate reverse pairs: (each original product as reactant, each original reactant as product)
        # This correctly represents the reverse reaction direction.
        reversed_combinations = list(itertools.product(product_ids, reactant_ids))
        all_combinations.extend(reversed_combinations)

    # Remove duplicate pairs if any (e.g., if R={A}, P={A} and reversible)
    # Though for distinct R and P lists, product usually won't create duplicates here.
    # If lists can have duplicates, set(all_combinations) might be needed,
    # but itertools.product handles unique elements from input lists fine.
    # The main source of duplicates would be if reversible adds identical pairs to forward.
    # Example: R=[1], P=[2], rev=True -> (1,2), (2,1)
    # Example: R=[1], P=[1], rev=True -> (1,1), (1,1) -> needs set for [(1,1)]
    # However, reaction_pairs table unique constraint should handle DB-side duplicates.
    # For clarity, let's remove Python-side duplicates if they arise from this logic.
    # Using list(dict.fromkeys(all_combinations)) is an ordered way to get unique tuples.

    unique_combinations = list(dict.fromkeys(all_combinations))

    logger.debug(
        f"Generated {len(unique_combinations)} unique reaction ID pairs for R_IDs:{reactant_ids}, P_IDs:{product_ids}, Rev:{reversible}.")
    return unique_combinations
