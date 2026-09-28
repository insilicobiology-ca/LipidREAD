import re
import psycopg
from psycopg.rows import dict_row
from typing import List, Dict, Optional, Set
from pathlib import Path
from src.data.config_loader import ConfigLoader
from src.models.lipid_components import LysoGlycerophospholipidComponents
from src.parsing.lipid_parser import LipidParser
import logging

logger = logging.getLogger(__name__)


class LipidTranslator:
    def __init__(self, config_loader: ConfigLoader, lipid_parser: LipidParser):
        self.config = config_loader
        self.lipid_parser = lipid_parser
        self.conn: Optional[psycopg.Connection] = None
        self.cursor: Optional[psycopg.Cursor] = None
        self._connect_to_db()

    # Constants for Hex translations
    HEX_PREFIX = "Hex"
    BETA_GLC_PREFIX = "beta-Glc"
    BETA_GAL_PREFIX = "beta-Gal"
    HEX2_PREFIX = "Hex2"
    LAC_PREFIX = "Lac"
    GALA_PREFIX = "Gala"
    DG_PREFIX = "DG"
    ZERO_CHAIN = "/0:0"

    def _connect_to_db(self) -> None:
        try:
            self.conn = psycopg.connect(
                dbname=self.config.get_db_name(),
                user=self.config.get_db_user(),
                password=self.config.get_db_password(),
                host=self.config.get_db_host(),
                port=self.config.get_db_port(),
                options=f"-c search_path=LipidTranslator",
            )
            self.cursor = self.conn.cursor(row_factory=dict_row)
            logger.info("Successfully connected to LipidTranslator database")
        except psycopg.Error as e:
            logger.error(f"Failed to connect to LipidTranslator database: {e}")
            self.conn = None
            self.cursor = None
            raise ConnectionError(f"LipidTranslator DB connection failed: {e}") from e

    def _execute_query(self, query: str, params: Dict[str, str]) -> List[str]:
        """Execute query and return list of names."""
        if not self.cursor or not self.conn or self.conn.closed:
            logger.error("Database connection not available")
            return []

        try:
            self.cursor.execute(query, params)
            results = self.cursor.fetchall()
            return [row["NAME"] for row in results if row and "NAME" in row]
        except psycopg.Error as e:
            logger.error(f"Database query failed: {e}")
            return []

    def translate_lipid(self, lipid_name: str) -> List[str]:
        """Translate a single lipid name using database lookup."""
        query = """
            SELECT LT."NAME" from "LipidTranslator".lipids AS LT 
            LEFT JOIN "LipidTranslator".hmdb AS LH on LT."HMDB_ID" = LH."HMDB_ID" 
            LEFT JOIN "LipidTranslator".pubchem LP on LT."PUBCHEM_CID" = LP."PUBCHEM_CID" 
            WHERE LT."LM_ID" ILIKE %(lipid)s 
            OR LT."NAME" ILIKE %(lipid)s 
            OR LT."SYSTEMATIC_NAME" ILIKE %(lipid)s 
            OR LT."EXACT_MASS" ILIKE %(lipid)s 
            OR LT."PUBCHEM_CID" ILIKE %(lipid)s 
            OR LT."HMDB_ID" ILIKE %(lipid)s 
            OR LT."CHEBI_ID" ILIKE %(lipid)s 
            OR LT."SWISSLIPIDS_ID" ILIKE %(lipid)s 
            OR LT."ABBREVIATION" ILIKE %(lipid)s 
            OR LT."SYNONYMS" ILIKE %(synonym)s 
            OR LT."FORMULA" ILIKE %(lipid)s
        """

        params = {"lipid": f"{lipid_name}", "synonym": f"%|{lipid_name}|%"}

        return self._execute_query(query, params)

    def translate_lipids(self, lipid_names: List[str]) -> Dict[str, List[str]]:
        """
        Translate a list of lipid names, handling Hex and Hex2 translations.
        Returns a dictionary mapping original names to their translations.
        """
        translation_tracker = {}

        for lipid in lipid_names:
            all_translations = set()

            # Get database translations
            try:
                translated_names = self.translate_lipid(lipid)
                if translated_names:
                    for name in translated_names:
                        cleaned_name = self.clean_lipid_abbreviation(name)
                        try:
                            components = self.lipid_parser.parse_lipid(cleaned_name)
                            if isinstance(
                                components, LysoGlycerophospholipidComponents
                            ):
                                cleaned_name = str(components)
                            all_translations.add(cleaned_name)
                        except ValueError:
                            logger.warning(
                                f"Could not parse lipid name '{cleaned_name}' - skipping"
                            )
            except ConnectionError:
                logger.warning(
                    "Database not available - skipping database translations"
                )

            # Get hex translations
            hex_translations = self._handle_hex_translations(lipid)
            if hex_translations:
                all_translations.update(
                    self.clean_lipid_abbreviation(name) for name in hex_translations
                )

            # Always include original lipid name
            all_translations.add(lipid)

            translation_tracker[lipid] = sorted(list(all_translations))

        return translation_tracker

    def _handle_hex_translations(self, lipid: str) -> Set[str]:
        """Handle both Hex and Hex2 translations for a given lipid."""
        translations = set()

        if lipid.startswith(self.HEX2_PREFIX):
            lac_version = lipid.replace(self.HEX2_PREFIX, self.LAC_PREFIX, 1)
            gala_version = lipid.replace(self.HEX2_PREFIX, self.GALA_PREFIX, 1)
            translations.update([lac_version, gala_version])
        elif lipid.startswith(self.HEX_PREFIX):
            glc_version = lipid.replace(self.HEX_PREFIX, self.BETA_GLC_PREFIX, 1)
            gal_version = lipid.replace(self.HEX_PREFIX, self.BETA_GAL_PREFIX, 1)
            translations.update([glc_version, gal_version])

        return translations

    def _is_dg_lipid(self, lipid: str) -> bool:
        """Check if the given lipid is a DG (Diacylglycerol) lipid."""
        return lipid.startswith(self.DG_PREFIX)

    def _clean_dg_lipid(self, lipid: str) -> str:
        """Remove the redundant '0:0' chain from DG lipids if present."""
        if not self._is_dg_lipid(lipid):
            return lipid
        if self.ZERO_CHAIN in lipid:
            return lipid.replace(self.ZERO_CHAIN, "")
        return lipid

    def clean_lipid_abbreviation(self, lipid_abbreviation: str) -> str:
        """
        Clean a lipid abbreviation by removing specific patterns and standardizing hydroxyl group notation.

        Args:
            lipid_abbreviation: The original lipid abbreviation.

        Returns:
            The cleaned lipid abbreviation.
        """
        if not isinstance(lipid_abbreviation, str):
            return str(lipid_abbreviation)

        # Remove double bond positions, 'iso-' prefix, [iso#], and [rac]
        cleaned_abbreviation = re.sub(
            r"\(\d+(E|Z)(,\d+(E|Z))*\)|iso-|\[iso\d\]|\[rac\]", "", lipid_abbreviation
        )

        # Standardize hydroxyl group notation to '(OH)'
        cleaned_abbreviation = re.sub(
            r"\((\d+)OH\)|-(\d+)OH", "(OH)", cleaned_abbreviation
        )

        cleaned_abbreviation = self._clean_dg_lipid(cleaned_abbreviation)

        return cleaned_abbreviation.strip()

    def close_connection(self) -> None:
        """Close the database connection."""
        if self.cursor:
            self.cursor.close()
            self.cursor = None
        if self.conn and not self.conn.closed:
            self.conn.close()
            self.conn = None

        try:
            logger.info("LipidTranslator database connection closed.")
        except Exception:
            pass
        logger.info("LipidTranslator database connection closed")

    def __del__(self):
        """Cleanup method - explicitly close connections."""
        self.close_connection()
