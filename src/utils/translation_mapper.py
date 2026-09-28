from src.models.models import Reaction
from typing import List, Dict, Tuple, Optional, Set
from itertools import product


class TranslationMapper:
    def __init__(self, translation_tracker: Dict[str, List[str]]):
        self.translation_tracker = translation_tracker
        # Create reverse mapping where values are sets of original names
        self.reverse_mapping: Dict[str, Set[str]] = {}

        for original, translations in translation_tracker.items():
            for translated in translations:
                if translated not in self.reverse_mapping:
                    self.reverse_mapping[translated] = set()
                self.reverse_mapping[translated].add(original)

    def get_original_names(self, translated_name: str) -> Set[str]:
        """
        Get all possible original names for a translated form.
        Returns the translated name itself if no mapping exists.
        """
        return self.reverse_mapping.get(translated_name, {translated_name})

    # def reverse_translate_reaction(self, reaction: Reaction) -> Reaction:
    #     return Reaction(
    #         reaction_id=reaction.reaction_id,
    #         reaction_enzyme_id=reaction.reaction_enzyme_id,
    #         enzyme_name=reaction.enzyme_name,
    #         uniprot_id=reaction.uniprot_id,
    #         reactants=[self.get_original_name(r) for r in reaction.reactants],
    #         product=self.get_original_name(reaction.product),
    #         rhea_id=reaction.rhea_id,
    #         doi=reaction.doi,
    #         reaction_type=reaction.reaction_type
    #     )
    def reverse_translate_reaction(self, reaction: Reaction, input_lipids: Set[str],
                                   validate_products: bool = True) -> List[Reaction]:
        """
        MODIFIED: Creates all possible original reaction forms for multi-product reactions.
        Filters results to only include combinations where at least one of the specific,
        non-class-level products exists in the input lipids.
        """
        # Get all possible original forms for each reactant and product
        reactant_options = [self.get_original_names(r) for r in reaction.reactants]
        product_options = [self.get_original_names(p) for p in reaction.products]

        valid_reactions = []

        # Generate all possible combinations of reactants and products
        for reactant_combo in product(*reactant_options):
            for product_combo in product(*product_options):

                # The validation logic has already been performed in LipidProcessor.
                # Here, we are just creating all possible valid name permutations.
                # The original input_lipids set contains the non-translated names, which is what we need.

                # Check if at least one of the reverse-translated products is in the original user input list.
                # This ensures we don't create reaction variants that lead to products the user never provided.
                if validate_products and not any(p in input_lipids for p in product_combo):
                    continue

                # Create new reaction with this combination
                valid_reactions.append(Reaction(
                    reaction_id=reaction.reaction_id,
                    reaction_enzyme_id=reaction.reaction_enzyme_id,
                    enzyme_name=reaction.enzyme_name,
                    uniprot_id=reaction.uniprot_id,
                    reactants=list(reactant_combo),
                    products=list(product_combo),  # MODIFIED: Use the product combo
                    rhea_id=reaction.rhea_id,
                    doi=reaction.doi,
                    reaction_type=reaction.reaction_type,
                    subcellular_locations=reaction.subcellular_locations,
                    generalized_from_ids=reaction.generalized_from_ids
                ))

        # Return a list of unique reactions, as some combinations might be identical
        # (e.g., if a lipid has no translation)
        unique_reactions = list({r.get_key(): r for r in valid_reactions}.values())
        return unique_reactions