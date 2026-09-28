from functools import lru_cache
from abc import ABC, abstractmethod
from typing import List, Set, Tuple, Optional, Dict
import psycopg  # Change from psycopg2 to psycopg
from psycopg.rows import dict_row  # Replaces DictCursor
from src.models.models import Reaction, ReactionType
from src.data.config_loader import ConfigLoader
import logging
from collections import defaultdict
import pandas as pd
from pandas import json_normalize
import csv

logger = logging.getLogger(__name__)


class LipidDataAccess(ABC):
    @abstractmethod
    def get_molecule_info(self, abbreviation: str) -> List[Tuple[int, str, str]]:
        pass

    @abstractmethod
    def get_reactions(self, molecule_id: int, organism_id: int) -> Dict[int, List[Reaction]]:
        # MODIFIED: Return type is now a dictionary grouping reactions by enzyme event.
        pass

    def get_synthesis_reactions_for_molecule(self, molecule_id: int, organism_id: int) -> Dict[int, List[Reaction]]:
        pass

    @abstractmethod
    def get_family_molecule(self, sl_lipid_class: str) -> Optional[Tuple[int, str, str]]:
        pass

    @abstractmethod
    def get_second_reactant(self, first_reactant_molecule_id: int, reaction_enzyme_id: int) -> Optional[Tuple]:
        pass

    @abstractmethod
    def get_all_participants(self, reaction_enzyme_id: int) -> Optional[Tuple[Set[str], Set[str]]]:
        # NEW: Abstract method for the new helper.
        pass

    @abstractmethod
    def get_super_reaction_ids(self) -> Set[str]:
        pass

    @abstractmethod
    def clear_caches(self) -> None:
        pass

class DatabaseLipidDataAccess(LipidDataAccess):
    def __init__(self, config_loader: ConfigLoader):
        self.config = config_loader
        self.conn = None
        self.cursor = None
        self._connect_to_database()
        self._super_reaction_ids: Set[int] = set(self.config.get_super_reaction_ids())
        self.query_cache: Dict[Tuple[int, int], Optional[Tuple]] = {}

    def _connect_to_database(self):
        try:
            self.conn = psycopg.connect(
                dbname=self.config.get_db_name(),
                user=self.config.get_db_user(),
                password=self.config.get_db_password(),
                host=self.config.get_db_host(),
                port=self.config.get_db_port(),
                options=f"-c search_path={self.config.get_db_schema()}"
            )
            self.conn.row_factory = psycopg.rows.dict_row
            self.cursor = self.conn.cursor()
        except psycopg.Error as e:
            logger.error(f"Failed to connect to the database: {str(e)}")
            raise

    @lru_cache(maxsize=5000)
    def get_molecule_info(self, identifier: str) -> List[Tuple[int, str, str]]:
        query = """
                    SELECT molecule_id, molecule_name, sl_lipid_class 
                    FROM molecules 
                    WHERE cleaned_abbreviation ILIKE %s

                    UNION

                    SELECT molecule_id, molecule_name, sl_lipid_class 
                    FROM molecules 
                    WHERE molecule_name ILIKE %s
                """
        # We pass the identifier twice, once for each part of the query.
        # ILIKE provides a case-insensitive match for the full name.
        params = (identifier, identifier)

        self.cursor.execute(query, params)

        rows = self.cursor.fetchall()
        return [(row['molecule_id'], row['molecule_name'], row['sl_lipid_class']) for row in rows]

    @lru_cache(maxsize=5000)
    def get_molecule_info_by_synonym(self, synonym: str) -> List[Tuple[int, str, str]]:
        """
        Retrieves molecule info by searching the molecule_synonyms field.
        Returns a list of tuples to be consistent with get_molecule_info.
        """
        query = """
                SELECT molecule_id, molecule_name, sl_lipid_class
                FROM molecules 
                WHERE molecule_synonyms ILIKE %s
            """
        search_pattern = f"%|{synonym}|%"

        self.cursor.execute(query, (search_pattern,))
        rows = self.cursor.fetchall()
        # Ensure the return type is identical to get_molecule_info
        return [(row['molecule_id'], row['molecule_name'], row['sl_lipid_class']) for row in rows]

    def get_full_molecule_record_by_slm_id(self, slm_id: str) -> Optional[Dict]:
        """
        Retrieves a full molecule record dictionary using its SwissLipids ID.
        """
        query = "SELECT * FROM molecules WHERE swisslipids_id = %s"
        self.cursor.execute(query, (slm_id,))
        return self.cursor.fetchone()  # Returns a single dict or None

    @lru_cache(maxsize=3000)
    def get_reactions(self, molecule_id: int, organism_id: int) -> Dict[int, List[Reaction]]:
        """
        MODIFIED: Retrieves reactions grouped by reaction_enzyme_id.
        Each group contains a list of reaction pairs that constitute a single biochemical event.
        """
        query = """SELECT r.reaction_id, re.reaction_enzyme_id, e.enzyme_name, 
                          mr.cleaned_abbreviation as reactant_abbreviation,
                          mp.cleaned_abbreviation as product_abbreviation,
                          r.rhea_id, r.doi, e.uniprot_id,
                          e.subcellular_location,
                          r.generalized_from_reaction_ids -- NEWLY ADDED COLUMN
                       FROM reaction_pairs AS rp
                       JOIN reaction_enzyme AS re ON rp.reaction_enzyme_id = re.reaction_enzyme_id
                       JOIN reactions AS r on re.reaction_id = r.reaction_id
                       JOIN molecules AS mr ON rp.reactant_molecule_id = mr.molecule_id
                       JOIN molecules AS mp ON rp.product_molecule_id = mp.molecule_id
                       JOIN enzymes AS e ON re.enzyme_id = e.enzyme_id
                       WHERE rp.reactant_molecule_id = %s AND e.organism_id = %s"""
        self.cursor.execute(query, (molecule_id, organism_id))

        grouped_reactions = defaultdict(list)
        for row in self.cursor.fetchall():
            reaction_enzyme_id = row['reaction_enzyme_id']
            reaction_type = ReactionType.SUPER if row[
                                                      'reaction_id'] in self._super_reaction_ids else ReactionType.REGULAR

            # This object represents a single reactant->product pair within a larger event.
            reaction_pair = Reaction(
                reaction_id=row['reaction_id'],
                reaction_enzyme_id=reaction_enzyme_id,
                enzyme_name=row['enzyme_name'],
                uniprot_id=row['uniprot_id'],
                reactants=[row['reactant_abbreviation']],
                products=[row['product_abbreviation']],  # Using new list format
                rhea_id=row['rhea_id'],
                doi=row['doi'],
                reaction_type=reaction_type,
                subcellular_locations = row['subcellular_location'],
                generalized_from_ids=row['generalized_from_reaction_ids']
            )
            grouped_reactions[reaction_enzyme_id].append(reaction_pair)

        return dict(grouped_reactions)

    @lru_cache(maxsize=3000)  # Caching is beneficial here too
    def get_synthesis_reactions_for_molecule(self, product_molecule_id: int, organism_id: int) -> Dict[
        int, List[Reaction]]:
        """
        Retrieves all reaction events that produce a specific molecule, grouped by reaction_enzyme_id.
        This is the reverse of get_reactions.

        Args:
            product_molecule_id: The molecule_id of the lipid in the product role.
            organism_id: The NCBI taxonomy ID of the organism.

        Returns:
            A dictionary where keys are reaction_enzyme_id and values are lists of
            Reaction objects representing the pairs in that event.
        """
        query = """
                SELECT r.reaction_id, re.reaction_enzyme_id, e.enzyme_name, 
                       mr.cleaned_abbreviation as reactant_abbreviation,
                       mp.cleaned_abbreviation as product_abbreviation,
                       r.rhea_id, r.doi, e.uniprot_id,
                       e.subcellular_location,
                       r.generalized_from_reaction_ids
                FROM reaction_pairs AS rp
                JOIN reaction_enzyme AS re ON rp.reaction_enzyme_id = re.reaction_enzyme_id
                JOIN reactions AS r on re.reaction_id = r.reaction_id
                JOIN molecules AS mr ON rp.reactant_molecule_id = mr.molecule_id
                JOIN molecules AS mp ON rp.product_molecule_id = mp.molecule_id
                JOIN enzymes AS e ON re.enzyme_id = e.enzyme_id
                WHERE rp.product_molecule_id = %s AND e.organism_id = %s
            """
        self.cursor.execute(query, (product_molecule_id, organism_id))

        grouped_reactions = defaultdict(list)
        for row in self.cursor.fetchall():
            reaction_enzyme_id = row['reaction_enzyme_id']
            reaction_type = ReactionType.SUPER if row[
                                                      'reaction_id'] in self._super_reaction_ids else ReactionType.REGULAR

            reaction_pair = Reaction(
                reaction_id=row['reaction_id'],
                reaction_enzyme_id=reaction_enzyme_id,
                enzyme_name=row['enzyme_name'],
                uniprot_id=row['uniprot_id'],
                reactants=[row['reactant_abbreviation']],
                products=[row['product_abbreviation']],
                rhea_id=row['rhea_id'],
                doi=row['doi'],
                reaction_type=reaction_type,
                subcellular_locations=row['subcellular_location'],
                generalized_from_ids=row['generalized_from_reaction_ids']
            )
            grouped_reactions[reaction_enzyme_id].append(reaction_pair)

        return dict(grouped_reactions)

    def get_all_participants(self, reaction_enzyme_id: int) -> Optional[Tuple[Set[str], Set[str]]]:
        """
        NEW: For a given reaction_enzyme_id, retrieves the complete set of
        reactant and product abbreviations from the reaction_pairs table.
        This is used for pre-specified molecular reactions.
        """
        if reaction_enzyme_id is None:
            return None

        try:
            # Query 1: Get all unique reactants for the event
            reactants_query = """
                SELECT DISTINCT mr.cleaned_abbreviation 
                FROM reaction_pairs rp 
                JOIN molecules mr ON rp.reactant_molecule_id = mr.molecule_id 
                WHERE rp.reaction_enzyme_id = %s
            """
            self.cursor.execute(reactants_query, (reaction_enzyme_id,))
            all_reactants = {row['cleaned_abbreviation'] for row in self.cursor.fetchall()}

            # Query 2: Get all unique products for the event
            products_query = """
                SELECT DISTINCT mp.cleaned_abbreviation 
                FROM reaction_pairs rp 
                JOIN molecules mp ON rp.product_molecule_id = mp.molecule_id 
                WHERE rp.reaction_enzyme_id = %s
            """
            self.cursor.execute(products_query, (reaction_enzyme_id,))
            all_products = {row['cleaned_abbreviation'] for row in self.cursor.fetchall()}

            if not all_reactants or not all_products:
                logger.warning(f"Could not find complete participants for REID {reaction_enzyme_id}")
                return None

            return all_reactants, all_products

        except Exception as e:
            logger.error(f"Failed to get all participants for REID {reaction_enzyme_id}: {e}")
            return None

    @lru_cache(maxsize=3000)
    def get_family_molecule(self, sl_lipid_class: str) -> Optional[Tuple[int, str, str]]:
        query = """SELECT molecule_id, molecule_name, sl_lipid_class 
                   FROM molecules 
                   WHERE swisslipids_id = %s
                   LIMIT 1"""
        self.cursor.execute(query, [sl_lipid_class])
        row = self.cursor.fetchone()
        return (row['molecule_id'], row['molecule_name'], row['sl_lipid_class']) if row else None


    def get_second_reactant(self, first_reactant_molecule_id: int, reaction_enzyme_id: int) -> Optional[Tuple]:
        cache_key = (first_reactant_molecule_id, reaction_enzyme_id)
        if cache_key in self.query_cache:
            return self.query_cache[cache_key]

        query = """SELECT rp.reactant_molecule_id, mr.molecule_name, mr.cleaned_abbreviation, e.enzyme_name, r.reaction_id
                    FROM reaction_pairs AS rp
                    JOIN reaction_enzyme AS re ON rp.reaction_enzyme_id = re.reaction_enzyme_id
                    JOIN reactions AS r on re.reaction_id = r.reaction_id
                    JOIN molecules AS mr ON rp.reactant_molecule_id = mr.molecule_id
                    JOIN molecules AS mp ON rp.product_molecule_id = mp.molecule_id
                    JOIN enzymes AS e ON re.enzyme_id = e.enzyme_id
                    WHERE rp.reactant_molecule_id != %s AND re.reaction_enzyme_id = %s"""
        try:
            self.cursor.execute(query, [first_reactant_molecule_id, reaction_enzyme_id])
            result = self.cursor.fetchone()
            self.query_cache[cache_key] = result
            return result
        except Exception as e:
            logger.error(f"Database error in get_second_reactant: {str(e)}")
            return None

    def get_super_reaction_ids(self) -> Set[int]:
        return self._super_reaction_ids

    def get_reactions_by_enzyme(self, enzyme_identifier: str, organism_id: int) -> List[Dict[str, str]]:
        """
        Retrieves all reaction texts associated with a given enzyme for a specific organism.
        The identifier can be a UniProt AC or an enzyme name (case-insensitive).

        Args:
            enzyme_identifier: The UniProt AC or enzyme name to search for.
            organism_id: The ID of the organism.

        Returns:
            A list of dictionaries, each containing reaction and enzyme details.
        """
        query = """
            SELECT DISTINCT r.reaction_text, r.rhea_id, e.enzyme_name, e.uniprot_id
            FROM lipograph.enzymes e
            JOIN lipograph.reaction_enzyme re ON e.enzyme_id = re.enzyme_id
            JOIN lipograph.reactions r ON re.reaction_id = r.reaction_id
            WHERE e.organism_id = %s 
              AND (e.uniprot_id = %s OR e.enzyme_name = %s)
            ORDER BY r.reaction_text
        """
        params = (organism_id, enzyme_identifier, enzyme_identifier)

        self.cursor.execute(query, params)

        # The row_factory=dict_row ensures each result is a dict.
        # We can just return the list of fetched results directly.
        results = self.cursor.fetchall()
        logger.info(
            f"Found {len(results)} reactions for enzyme identifier '{enzyme_identifier}' in organism {organism_id}.")

        return results

    def __del__(self):
        """Close resources when object is deleted"""
        if hasattr(self, 'cursor') and self.cursor:
            self.cursor.close()
        if hasattr(self, 'conn') and self.conn:
            self.conn.close()

    def clear_caches(self):
        """Clear all method caches"""
        self.get_molecule_info.cache_clear()
        self.get_reactions.cache_clear()
        self.get_family_molecule.cache_clear()
        self.query_cache.clear()
