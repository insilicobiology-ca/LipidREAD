import logging
import re, yaml
from functools import lru_cache
from typing import Optional, Dict, Any, List, Tuple
from pathlib import Path

from src.data.data_access import DatabaseLipidDataAccess
from src.parsing.lipid_parser import LipidParser
from src.utils.lipid_translator import LipidTranslator

logger = logging.getLogger(__name__)

CONTEXT_RESOLUTIONS_KEY = "context_specific_resolutions"
GLOBAL_NAME_RESOLUTIONS_KEY = "global_name_resolutions"
NO_MATCH_MARKER = "__NO_MATCH_CONFIRMED__"

class MoleculeResolver:
    """
    A dedicated service for resolving chemical names from reaction texts
    into specific molecule entries from the database.
    """

    def __init__(
        self,
        data_access: DatabaseLipidDataAccess,
        lipid_parser: LipidParser,
        lipid_translator: LipidTranslator,
        resolutions_filepath: Optional[Path] = None
    ):
        self.data_access = data_access
        self.lipid_parser = lipid_parser
        self.lipid_translator = lipid_translator

        self.context_resolutions, self.global_name_resolutions = self._load_resolutions_from_yaml(resolutions_filepath)

    def _load_resolutions_from_yaml(self, yaml_path_str: Optional[str]) -> Tuple[Dict[str, str], Dict[str, str]]:
        """Loads manual resolutions from the given filepath."""
        context_res, global_res = {}, {}

        if not yaml_path_str:
            logger.info("No resolutions file path provided. No manual resolutions loaded.")
            return context_res, global_res

        yaml_path = Path(yaml_path_str)
        if not yaml_path.exists():
            logger.warning(f"Resolution file not found at '{yaml_path}'. No manual resolutions loaded.")
            return context_res, global_res

        logger.info(f"Loading manual resolutions from '{yaml_path}'...")
        try:
            with open(yaml_path, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
                if data:
                    context_res = data.get("context_specific_resolutions", {})
                    global_res = data.get("global_name_resolutions", {})
                    logger.info(
                        f"Loaded {len(context_res)} context-specific and {len(global_res)} global name resolutions.")
        except Exception as e:
            logger.error(f"Failed to load or parse resolution YAML file '{yaml_path}': {e}.")

        return context_res, global_res

    def _check_manual_resolutions(self, name: str, rhea_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """
        Checks for a manually curated resolution in the loaded YAML file.
        Returns the molecule info dictionary if a match is found, None otherwise.
        Can return a special marker to indicate a confirmed non-match.
        """
        # 1. Check for a context-specific match (name::Rhea_ID)
        if rhea_id:
            context_key = f"{name}::{rhea_id}"
            if context_key in self.context_resolutions:
                resolution = self.context_resolutions[context_key]
                logger.debug(f"Resolver: Found CONTEXT manual resolution for '{context_key}' -> '{resolution}'")
                if resolution == NO_MATCH_MARKER: return {'status': NO_MATCH_MARKER}
                return self.data_access.get_full_molecule_record_by_slm_id(resolution)

        # 2. Check for a global match (name)
        if name in self.global_name_resolutions:
            resolution = self.global_name_resolutions[name]
            logger.debug(f"Resolver: Found GLOBAL manual resolution for '{name}' -> '{resolution}'")
            if resolution == NO_MATCH_MARKER: return {'status': NO_MATCH_MARKER}
            return self.data_access.get_full_molecule_record_by_slm_id(resolution)

        return None

    @lru_cache(maxsize=4096)
    def resolve_to_class(self, participant_name: str, rhea_id: Optional[str] = None) -> str:
        """
        Takes a full chemical name and resolves it to an abbreviated lipid class.
        If it cannot be resolved, it returns the cleaned original name.
        """
        if not participant_name or not participant_name.strip():
            return ""

        name_to_resolve = self._clean_participant_name(participant_name)

        molecule_info = self._resolve_name_to_molecule_info(name_to_resolve, rhea_id)

        if molecule_info:
            abbreviation = molecule_info.get("cleaned_abbreviation")
            if abbreviation:
                try:
                    components = self.lipid_parser.parse_lipid(abbreviation)
                    # Return the base headgroup, stripping positional info like _sn1
                    return components.headgroup
                except (ValueError, TypeError):
                    return (
                        abbreviation  # Fallback to the abbreviation if it's unparsable
                    )

        return name_to_resolve

    def _clean_participant_name(self, name: str) -> str:
        """Removes stoichiometric coefficients and common prefixes."""
        cleaned = re.sub(r"^\d+\s+", "", name.strip())
        if cleaned.lower().startswith(("a ", "an ")):
            cleaned = cleaned[2:].strip()
        return cleaned

    def _resolve_name_to_molecule_info(self, name: str, rhea_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """
        The core multi-stage resolver logic.

        MODIFIED:
        1. Adds a special rule to handle all "-CoA" names by assigning them to the generic Fatty Acyl-CoA class.
        2. Adds a fallback to look up the parent class abbreviation if the direct molecule is missing one.
        """

        # --- NEW: Rule 0 - Check manual YAML resolutions first ---
        manual_resolution = self._check_manual_resolutions(name, rhea_id)
        if manual_resolution:
            if manual_resolution.get('status') == NO_MATCH_MARKER:
                return None  # This is a confirmed non-match, stop processing.
            return manual_resolution  # Return the full record from the manual lookup.


        # Rule 1: Handle -CoA (This is a name-based override, so it stays first)
        if name.upper().endswith(('-COA', ' COA')):
            logger.debug(f"Resolved '{name}' as a generic Fatty Acyl-CoA using special rule.")
            return {'molecule_id': 55, 'molecule_name': 'Fatty Acyl-CoA', 'sl_lipid_class': 'SLM:000000249',
                    'cleaned_abbreviation': 'FA-CoA'}

        # --- Step 1: Find the molecule and its best available abbreviation ---
        molecule_infos = self.data_access.get_molecule_info(name)
        if not molecule_infos:
            molecule_infos = self.data_access.get_molecule_info_by_synonym(name)

        if not molecule_infos:
            return None

        mol_id, mol_name, sl_lipid_class = molecule_infos[0]

        # Find the best abbreviation: either its own or its parent's.
        abbreviation_query = "SELECT cleaned_abbreviation FROM molecules WHERE molecule_id = %s"
        self.data_access.cursor.execute(abbreviation_query, (mol_id,))
        abbreviation_result = self.data_access.cursor.fetchone()
        final_abbr = abbreviation_result['cleaned_abbreviation'] if abbreviation_result else None

        if not final_abbr and sl_lipid_class:
            # Fallback to parent if direct abbreviation is missing
            logger.debug(
                f"Molecule ID {mol_id} ('{mol_name}') is missing an abbreviation. Attempting fallback to parent class '{sl_lipid_class}'.")
            parent_query = "SELECT cleaned_abbreviation FROM molecules WHERE swisslipids_id = %s"
            self.data_access.cursor.execute(parent_query, (sl_lipid_class,))
            parent_info = self.data_access.cursor.fetchone()
            if parent_info and parent_info.get('cleaned_abbreviation'):
                final_abbr = parent_info['cleaned_abbreviation']

        # --- Step 2: Now, with the definitive abbreviation (or None), apply rules ---

        # Rule 2: Check for SFA/MUFA/PUFA abbreviations
        if final_abbr and final_abbr.upper() in ('SFA', 'MUFA', 'PUFA'):
            logger.debug(f"Resolved '{name}' (abbr: {final_abbr}) to generic Fatty Acid class.")
            return {'molecule_id': 613, 'molecule_name': 'Fatty Acids', 'sl_lipid_class': 'SLM:000000984',
                    'cleaned_abbreviation': 'FA'}

        # Rule 3: Check if abbreviation has chain info (species-level)
        if final_abbr and re.search(r'\(\d+:\d+', final_abbr) and sl_lipid_class:
            logger.debug(f"Abbreviation '{final_abbr}' is species-level. Resolving to parent class '{sl_lipid_class}'.")
            parent_query = "SELECT * FROM molecules WHERE swisslipids_id = %s"
            self.data_access.cursor.execute(parent_query, (sl_lipid_class,))
            parent_info = self.data_access.cursor.fetchone()
            if parent_info:
                # --- KEY FIX: Perform the SFA/MUFA/PUFA check AGAIN on the PARENT's abbreviation ---
                parent_abbr = parent_info.get('cleaned_abbreviation')
                if parent_abbr and parent_abbr.upper() in ('SFA', 'MUFA', 'PUFA'):
                    logger.debug(f"Parent class is '{parent_abbr}'. Resolving to generic Fatty Acid class.")
                    return {'molecule_id': 613, 'molecule_name': 'Fatty Acids', 'sl_lipid_class': 'SLM:000000984',
                            'cleaned_abbreviation': 'FA'}
                # ---------------------------------------------------------------------------------

                # If the parent is not SFA/MUFA/PUFA, return the parent's full record.
                return dict(parent_info)

        # If no special rules applied, return the molecule's info with the best abbreviation we found.
        return {
            'molecule_id': mol_id,
            'molecule_name': mol_name,
            'sl_lipid_class': sl_lipid_class,
            'cleaned_abbreviation': final_abbr  # This can be None if no abbreviation was ever found
        }
