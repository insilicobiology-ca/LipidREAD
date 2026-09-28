import pandas as pd
import psycopg  # Change from psycopg2 to psycopg
from psycopg.rows import dict_row  # Replaces DictCursor
import os
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any
import io
import re
import itertools
import collections
import requests
from thefuzz import fuzz
import json
import sys
import codecs
import functools
from functools import lru_cache
from typing import NamedTuple
import yaml

# Assuming ConfigLoader and potentially LipidTranslator are needed
from src.data.config_loader import ConfigLoader
from src.utils.text_utils import split_reaction_text, generate_reaction_combinations


def _to_int_or_none(val_str: Any) -> Optional[int]:
    if pd.isna(val_str) or str(val_str).strip().lower() in ['', 'none', 'null']:
        return None
    try:
        return int(float(str(val_str).strip()))  # float first handles "123.0"
    except ValueError as e:
        logging.getLogger(__name__).error(f"Invalid integer value for update: '{val_str}'. Error: {e}")
        raise ValueError(f"Invalid integer value for update: '{val_str}'")


def _process_string_value_for_update(val_str: Any) -> Optional[str]:
    if pd.isna(val_str):
        return None
    s_val = str(val_str).strip()
    if s_val.lower() in ['null', 'none', '']:  # If blank or literal "NULL"/"None" (case-insensitive)
        return None
    return str(val_str)  # Return original string (not stripped version if only "null" or "none" matched)
    # This allows strings like "None" (the word) to be stored if not intended as NULL marker.
    # Let's refine: if str(val_str).strip().lower() is 'null' or 'none', it's NULL.
    # if str(val_str).strip() is '', it's also NULL for non-text or optional text.
    # For a text field, if user provides "  ", it should be "  ".
    # The goal is to map reviewer's intent (blanking means NULL) correctly.


# Refined string processor:
def _process_approved_string_value(val_from_csv: Any) -> Optional[str]:
    if pd.isna(val_from_csv):  # Handles NaN from pandas
        return None

    val_str = str(val_from_csv)  # Work with string representation

    # If reviewer explicitly typed "NULL" or "None" (case-insensitive) or left it empty after stripping
    if val_str.strip().lower() in ['null', 'none', '']:
        return None

    # Otherwise, return the string as is (it might have leading/trailing spaces if intended by reviewer)
    return val_str

class APIMatch(NamedTuple):
    score: int
    slm_id: str
    name: str


# Configuration for field updates
# This should ideally be a class attribute or loaded from config.
# For now, defining it here for clarity within the class context.
FIELD_UPDATE_CONFIG = {
    "molecule": {
        "table_name": "lipograph.molecules",
        "pk_column_db": "swisslipids_id",  # DB column name for PK
        "pk_column_csv": "entity_key",  # CSV column name for PK value
        "pk_converter": _process_approved_string_value,  # SwissLipids IDs are strings
        "fields": {
            # field_name_csv: {"db_column": db_col_name, "converter": converter_func}
            "molecule_name": {"db_column": "molecule_name", "converter": _process_approved_string_value},
            "abbreviation": {"db_column": "abbreviation", "converter": _process_approved_string_value},
            "sl_lipid_class": {"db_column": "sl_lipid_class", "converter": _process_approved_string_value},
            # If 'cleaned_abbreviation' were to be updatable via this file (it shouldn't directly, it's derived),
            # it would also need an entry. Its update should be triggered by 'abbreviation' change.
        }
    },
    "enzyme": {
        "table_name": "lipograph.enzymes",
        "pk_column_db": "uniprot_id",
        "pk_column_csv": "entity_key",
        "pk_converter": _process_approved_string_value,  # UniProt IDs are strings
        "fields": {
            "enzyme_name": {"db_column": "enzyme_name", "converter": _process_approved_string_value},
            "organism_id": {"db_column": "organism_id", "converter": _to_int_or_none},
            "swisslipids_p_id": {"db_column": "swisslipids_p_id", "converter": _process_approved_string_value}
        }
    },
    "reaction": {
        "table_name": "lipograph.reactions",
        "pk_column_db": "reaction_id",  # PK in DB
        "pk_column_csv": "entity_key",  # PK from CSV
        "pk_converter": _to_int_or_none,  # Reaction ID is int
        "fields": {
            "rhea_id": {"db_column": "rhea_id", "converter": _to_int_or_none},
            "doi": {"db_column": "doi", "converter": _process_approved_string_value}
        }
    }
}

# --- Add these constants at the top of the file ---
YAML_RESOLUTIONS_FILE_NAME = "api_fallback_resolutions.yaml"
CONTEXT_RESOLUTIONS_KEY = "context_specific_resolutions"
GLOBAL_NAME_RESOLUTIONS_KEY = "global_name_resolutions"
NO_MATCH_MARKER = "__NO_MATCH_CONFIRMED__"


class LipidDatabaseUpdater:
    """
        Handles the process of updating the LipidCRED database
        from SwissLipids TSV files.
    """
    FIELD_UPDATE_MAPPINGS = FIELD_UPDATE_CONFIG
    # --- TEMPORARY: For collecting unresolved names ---
    _temp_unresolved_molecule_names = collections.defaultdict(int)

    # Store name and count of occurrences
    # --- END TEMPORARY ---

    def __init__(self, config_loader: ConfigLoader):
        # ...
        self.config = config_loader
        self.output_dir = Path(self.config.get_update_output_dir())
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

        from src.utils.lipid_translator import LipidTranslator
        from src.parsing.lipid_parser import LipidParser, LipidParserFactory, LipidComponentCache
        # Ensure these are initialized correctly based on your project structure
        # Assuming LipidParser and its dependencies are correctly set up if needed by LipidTranslator
        try:
            parser_factory = LipidParserFactory(self.config)
            cache = LipidComponentCache()
            lipid_parser_instance = LipidParser(factory=parser_factory, components_cache=cache)
            self.lipid_translator = LipidTranslator(self.config, lipid_parser_instance)
        except Exception as e:
            self.logger.error(f"Error initializing LipidTranslator or its dependencies: {e}", exc_info=True)
            # Depending on how critical LipidTranslator is for other parts, you might re-raise or handle.
            # For _apply_approved_field_updates, it's not directly used, but good to ensure init path is fine.
            self.lipid_translator = None # Ensure it's defined even if init fails, for other methods.

        self.family_molecule_id_map = self._load_family_molecule_id_map(config_loader)
        self.logger.info(f"Loaded {len(self.family_molecule_id_map)} entries into family_molecule_id_map.")

        self.conn: Optional[psycopg.Connection] = None
        self.cursor: Optional[psycopg.Cursor] = None
        self.logger.info(f"LipidDatabaseUpdater instance created. Output directory: {self.output_dir}")
        LipidDatabaseUpdater._temp_unresolved_molecule_names.clear()

        # ## NEW: Load the resolution file ##
        self.context_resolutions, self.global_name_resolutions = self._load_resolutions_from_yaml()

        # ## NEW: Instance variable for unresolved names ##
        self.unresolved_names_this_run = collections.defaultdict(int)

        # ## NEW: Dictionary to hold summary of new additions ##
        self.new_additions_summary = collections.defaultdict(list)

        self.abbreviation_warnings = []

        self.api_slp_cache = {}

    # ## NEW HELPER METHOD ##
    def _load_resolutions_from_yaml(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        """Loads manual resolutions from the YAML file."""
        yaml_path = self.output_dir / YAML_RESOLUTIONS_FILE_NAME
        context_resolutions, global_name_resolutions = {}, {}

        if not yaml_path.exists():
            self.logger.info(f"Resolution file not found at '{yaml_path}'. Starting with empty resolutions.")
            return context_resolutions, global_name_resolutions

        self.logger.info(f"Loading manual resolutions from '{yaml_path}'...")
        try:
            with open(yaml_path, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
                if data:
                    context_resolutions = data.get(CONTEXT_RESOLUTIONS_KEY, {})
                    global_name_resolutions = data.get(GLOBAL_NAME_RESOLUTIONS_KEY, {})
                    self.logger.info(
                        f"Loaded {len(context_resolutions)} context-specific and {len(global_name_resolutions)} global name resolutions.")
        except Exception as e:
            self.logger.error(
                f"Failed to load or parse resolution YAML file: {e}. Continuing with empty resolutions.")

        return context_resolutions, global_name_resolutions

    def _connect_to_database(self):
        """Establishes and returns a database connection and cursor."""
        if self.conn is None or self.conn.closed:
            try:
                conn = psycopg.connect(
                    dbname=self.config.get_db_name(),
                    user=self.config.get_db_user(),
                    password=self.config.get_db_password(),
                    host=self.config.get_db_host(),
                    port=self.config.get_db_port(),
                    options=f"-c search_path={self.config.get_db_schema()}"
                )
                # Use dict_row factory if you want dict-like row access
                # cur = conn.cursor(row_factory=dict_row)
                cur = conn.cursor()
                self.logger.info("Database connection established successfully.")
                return conn, cur
            except psycopg.Error as e:
                self.logger.error(f"Failed to connect to the database: {e}")
                raise ConnectionError(f"Database connection failed: {e}") from e
        return self.conn, self.cursor

    def __enter__(self):
        """Context manager entry: connect to the database."""
        self.logger.debug("Entering context: connecting to database.")
        if self.conn is None or self.conn.closed:
            try:
                self.conn = psycopg.connect(
                    dbname=self.config.get_db_name(),
                    user=self.config.get_db_user(),
                    password=self.config.get_db_password(),
                    host=self.config.get_db_host(),
                    port=self.config.get_db_port(),
                    options=f"-c search_path={self.config.get_db_schema()}"
                )
                # Use dict_row factory if you want dict-like row access
                # self.cursor = self.conn.cursor(row_factory=dict_row)
                self.cursor = self.conn.cursor()
                self.logger.info("Database connection established via context manager.")
            except psycopg.Error as e:
                self.logger.error(f"Failed to connect to the database in __enter__: {e}")
                if self.cursor: self.cursor.close()
                if self.conn: self.conn.close()
                self.cursor, self.conn = None, None
                raise ConnectionError(f"Database connection failed in __enter__: {e}") from e
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logger.debug("Exiting context: closing database connection (psycopg3).")
        # --- TEMPORARY: Print unresolved names at the end ---
        if LipidDatabaseUpdater._temp_unresolved_molecule_names:
            self.logger.info("=" * 50)
            self.logger.info("POTENTIAL METABOLITES / UNRESOLVED MOLECULE NAMES (and occurrence count):")
            # Sort by name for consistency, or by count for frequency
            sorted_unresolved = sorted(LipidDatabaseUpdater._temp_unresolved_molecule_names.items(),
                                       key=lambda item: item[0])
            for name, count in sorted_unresolved:
                self.logger.info(f"- \"{name}\"  (occurs {count} times)")
            self.logger.info(
                "Consider adding these to 'ignore_metabolites_in_pairs' in config.yaml if they are not lipids.")
            self.logger.info("=" * 50)
            LipidDatabaseUpdater._temp_unresolved_molecule_names.clear()  # Clear for next potential run in same session (not typical for CLI)
        # --- END TEMPORARY ---
        if self.cursor:
            self.cursor.close()
            self.logger.debug("Database cursor closed.")
        if self.conn:
            if exc_type is not None: # An exception occurred
                try:
                    # For psycopg3, check transaction status with conn.info.transaction_status
                    # psycopg.pq.TransactionStatus.INTRANS, INERROR, UNKNOWN
                    from psycopg.pq import TransactionStatus
                    if not self.conn.closed and self.conn.info.transaction_status in (TransactionStatus.INTRANS, TransactionStatus.INERROR):
                        self.logger.warning(f"Exception occurred in 'with' block ({exc_type}), rolling back transaction.")
                        self.conn.rollback()
                except Exception as rb_err:
                    self.logger.error(f"Error during rollback on exception: {rb_err}")
            self.conn.close()
            self.logger.info("Database connection closed via context manager (psycopg3).")
        self.cursor, self.conn = None, None
        return False

    def _load_family_molecule_id_map(self, config_loader: ConfigLoader) -> Dict[str, int]:
        family_map_from_config = config_loader.get_family_molecule_ids()
        if family_map_from_config:
            return {str(k): int(v) for k, v in family_map_from_config.items()} # Ensure types

    def report_missing_abbreviations(self, output_file: str):
        """
        Queries the production database for molecules missing an abbreviation
        or a cleaned_abbreviation and saves the report to a file.
        """
        self.logger.info("Generating report on molecules with missing abbreviations...")

        query = """
            SELECT swisslipids_id, molecule_name, abbreviation, cleaned_abbreviation
            FROM lipograph.molecules
            WHERE abbreviation IS NULL OR cleaned_abbreviation IS NULL
               OR TRIM(abbreviation) = '' OR TRIM(cleaned_abbreviation) = '';
        """
        df = self._query_to_df(query)

        if df.empty:
            self.logger.info("No molecules found with missing abbreviations. Excellent!")
            return

        output_path = Path(output_file)
        df.to_csv(output_path, sep='\t', index=False)
        self.logger.warning(f"Found {len(df)} molecules with missing/empty abbreviations. "
                            f"Report saved to: {output_path}")

    def _apply_approved_field_updates(self, approved_updates_file: Path) -> int:
        """
        Applies approved field updates from the specified CSV file to production tables.
        This method expects to be run within an active database transaction.

        Args:
            approved_updates_file (Path): Path to the 'approved_delta_updated_entries.csv' file.

        Returns:
            int: The number of successful update statements executed.
        """
        if not approved_updates_file.exists():
            self.logger.info(f"No approved field updates file found: {approved_updates_file}. Skipping.")
            return 0

        self.logger.info(f"Applying approved field updates from: {approved_updates_file}")
        try:
            # Read with dtype=str to preserve original string forms, na_filter=False to read empty strings as empty.
            # Or use na_values to specify what should be NaN, then converters handle NaN.
            # Let's use default NaN handling and let converters manage pd.isna().
            updates_df = pd.read_csv(approved_updates_file, sep='\t', dtype=str, keep_default_na=True,
                                     na_values=['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN',
                                                '-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>',
                                                'N/A', 'NA', 'NULL', 'NaN', 'n/a', 'nan', 'null'])
        except Exception as e:
            self.logger.error(f"Failed to read approved updates file {approved_updates_file}: {e}")
            raise  # Propagate error to roll back transaction

        if updates_df.empty:
            self.logger.info("Approved field updates file is empty. No updates to apply.")
            return 0

        update_count = 0
        for index, row in updates_df.iterrows():
            try:
                entity_type = row.get('entity_type')
                field_name_csv = row.get('field_name')  # Field name as defined in the CSV/delta logic

                # Value to update with (potentially modified by reviewer)
                new_value_from_csv = row.get('new_value_from_current_sl')

                if pd.isna(entity_type) or pd.isna(field_name_csv):
                    self.logger.warning(f"Skipping row {index + 2} in {approved_updates_file} due to missing "
                                        f"'entity_type' or 'field_name'. Row: {row.to_dict()}")
                    continue

                table_config = self.FIELD_UPDATE_MAPPINGS.get(entity_type)
                if not table_config:
                    self.logger.warning(f"No update configuration found for entity_type '{entity_type}' "
                                        f"(row {index + 2}). Skipping update. Row: {row.to_dict()}")
                    continue

                field_config = table_config["fields"].get(field_name_csv)
                if not field_config:
                    self.logger.warning(f"No field configuration for field_name '{field_name_csv}' "
                                        f"under entity_type '{entity_type}' (row {index + 2}). Skipping. Row: {row.to_dict()}")
                    continue

                # Get and convert Primary Key value for the WHERE clause
                pk_csv_col_name = table_config["pk_column_csv"]
                entity_key_csv_val = row.get(pk_csv_col_name)
                if pd.isna(entity_key_csv_val):
                    self.logger.warning(f"Missing entity_key ('{pk_csv_col_name}') for update in row {index + 2}. "
                                        f"Skipping. Row: {row.to_dict()}")
                    continue

                pk_converter = table_config["pk_converter"]
                db_pk_value = pk_converter(entity_key_csv_val)
                if db_pk_value is None and not (
                        pk_converter == _to_int_or_none or pk_converter == _process_approved_string_value):  # Check if None is legitimate for this PK type
                    # Most PKs are not nullable. If converter yields None, it's likely an issue.
                    self.logger.error(f"Primary key value '{entity_key_csv_val}' for {entity_type}.{field_name_csv} "
                                      f"converted to None, which is usually invalid for a PK. Skipping update for row {index + 2}.")
                    continue

                # Get and convert the new value for the SET clause
                value_converter = field_config["converter"]
                db_new_value = value_converter(new_value_from_csv)

                db_table_name = table_config["table_name"]
                db_column_to_update = field_config["db_column"]
                db_pk_column = table_config["pk_column_db"]

                # Construct and execute UPDATE statement
                # Ensure column names are quoted if they might contain special characters or are keywords
                # For simplicity, assuming well-behaved column names from config. If not, quote them: f'"{db_column_to_update}"'
                update_sql = f"""
                    UPDATE {db_table_name}
                    SET "{db_column_to_update}" = %s
                    WHERE "{db_pk_column}" = %s;
                """

                self.logger.debug(f"Executing update: {update_sql} with params ({db_new_value}, {db_pk_value})")
                self.cursor.execute(update_sql, (db_new_value, db_pk_value))

                if self.cursor.rowcount == 0:
                    self.logger.warning(
                        f"Update for {entity_type} with key '{db_pk_value}' and field '{db_column_to_update}' "
                        f"did not affect any rows. The entity might have been deleted or key changed. "
                        f"Row {index + 2}: {row.to_dict()}")
                elif self.cursor.rowcount == 1:
                    update_count += 1
                    self.logger.info(
                        f"Successfully updated {entity_type}.{db_column_to_update} for PK {db_pk_value} to '{db_new_value}'.")
                else:  # Should not happen for PK-based update
                    self.logger.warning(
                        f"Update for {entity_type} with key '{db_pk_value}' and field '{db_column_to_update}' "
                        f"affected {self.cursor.rowcount} rows. Expected 1. Row {index + 2}: {row.to_dict()}")
                    update_count += self.cursor.rowcount


            except ValueError as ve:  # Catch conversion errors specifically
                self.logger.error(f"Data conversion error processing update for row {index + 2} "
                                  f"in {approved_updates_file}: {ve}. Row: {row.to_dict()}",
                                  exc_info=False)  # Keep log clean from full tb for conv.
                # Continue to next row or re-raise to stop all updates?
                # For now, let's make it strict: a bad row in approved file stops the process.
                raise ValueError(
                    f"Error in approved_delta_updated_entries.csv, row {index + 2}: {ve}. Cannot apply updates.") from ve
            except Exception as e:
                self.logger.error(f"Error processing update for row {index + 2} "
                                  f"in {approved_updates_file}: {e}. Row: {row.to_dict()}", exc_info=True)
                raise  # Re-raise to ensure transaction rollback

        self.logger.info(f"Applied {update_count} field updates from {approved_updates_file}.")
        return update_count

    def _load_tsv_to_staging(self,
                             tsv_file_path: str,
                             table_name: str,
                             staging_table_columns_ordered: List[str],
                             file_specific_encoding: str,
                             tsv_column_rename_map: Optional[Dict[str, str]] = None,
                             text_replacements_in_cols: Optional[Dict[str, List[Tuple[str, str]]]] = None):
        self.logger.info(
            f"Preparing to load '{tsv_file_path}' into staging table '{table_name}' (file encoding: {file_specific_encoding}) using psycopg3.")
        if self.conn is None or self.conn.closed or self.cursor is None or self.cursor.closed:
            self.logger.error("Database connection not available in _load_tsv_to_staging.")
            raise ConnectionError("Database connection not established.")

        try:
            self.cursor.execute("SET synchronous_commit = OFF;")
            # 1. Preprocess the TSV data into a BytesIO buffer (UTF-8 encoded)
            processed_data_bytes_buffer = self._preprocess_tsv_for_copy(
                tsv_file_path,
                desired_columns_ordered=staging_table_columns_ordered,
                file_encoding=file_specific_encoding,
                rename_map=tsv_column_rename_map,
                text_replacements=text_replacements_in_cols
            )

            # 2. Truncate the staging table
            self.logger.debug(f"Truncating staging table '{table_name}'...")
            self.cursor.execute(f"TRUNCATE TABLE {table_name};")
            self.logger.info(f"Staging table '{table_name}' truncated.")

            # 3. Load data using new psycopg3 COPY interface
            # The `COPY ... ENCODING 'UTF8'` tells PostgreSQL to expect UTF-8 bytes from the client.
            # Our `processed_data_bytes_buffer` now provides exactly that.
            copy_sql = f"COPY {table_name} FROM STDIN WITH (FORMAT CSV, DELIMITER E'\\t', HEADER TRUE, NULL '', ENCODING 'UTF8', FREEZE TRUE);"

            self.logger.debug(f"Executing COPY command for '{table_name}' with psycopg3...")
            # `cursor.copy()` is a context manager itself for the COPY operation
            with self.cursor.copy(copy_sql) as copy:
                # Read from the BytesIO buffer in chunks and write to the COPY operation
                # This is more memory efficient for very large buffers than copy.write(buffer.readall())
                chunk_size = 1024 * 1024  # Adjust as needed
                while True:
                    chunk = processed_data_bytes_buffer.read(chunk_size)
                    if not chunk:
                        break
                    copy.write(chunk)

            # copy.rowcount might be available here if needed, after the 'with copy:' block
            self.logger.info(
                f"COPY operation to '{table_name}' completed (row count from copy not directly fetched here).")

            # 4. Commit the transaction (which includes TRUNCATE and COPY)
            self.conn.commit()

            self.cursor.execute("SET synchronous_commit = ON;")

            self.cursor.execute(f"SELECT COUNT(*) FROM {table_name};")
            row_count = self.cursor.fetchone()
            actual_row_count = row_count[0] if row_count else 0
            self.logger.info(f"Successfully loaded {actual_row_count} rows into '{table_name}'. Transaction committed.")

        except psycopg.Error as e:  # Catch psycopg.Error for psycopg3
            if self.conn and not self.conn.closed:
                self.conn.rollback()
                try:
                    self.cursor.execute("SET synchronous_commit = ON;")
                except:
                    pass  # Ignore errors when restoring settings
                self.logger.error(
                    f"Database error during loading into '{table_name}' (psycopg3): {e.diag.message_primary if e.diag else e}",
                    exc_info=True)
            raise
        except Exception as e:
            if self.conn and not self.conn.closed:
                self.conn.rollback()
                try:
                    self.cursor.execute("SET synchronous_commit = ON;")
                except:
                    pass  # Ignore errors when restoring settings
                self.logger.error(f"An unexpected error occurred while loading '{tsv_file_path}' into '{table_name}': {e}",
                                  exc_info=True)
            raise

    def _preprocess_tsv_for_copy(self,
                                 tsv_file_path: str,
                                 desired_columns_ordered: List[str],
                                 file_encoding: str = 'utf-8',
                                 rename_map: Optional[Dict[str, str]] = None,
                                 text_replacements: Optional[Dict[str, List[Tuple[str, str]]]] = None
                                 ) -> io.BytesIO:
        self.logger.info(
            f"Preprocessing TSV file: '{tsv_file_path}' for columns: {desired_columns_ordered}. "
            f"Attempting with specified encoding: '{file_encoding}'.")

        if not Path(tsv_file_path).exists():
            self.logger.error(f"TSV file not found for preprocessing: {tsv_file_path}")
            raise FileNotFoundError(f"TSV file not found: {tsv_file_path}")

        common_read_csv_args = {
            "sep": '\t', "dtype": str, "keep_default_na": False,
            "na_values": [''], "low_memory": False
        }

        df = None
        try:
            df = pd.read_csv(tsv_file_path, encoding=file_encoding, **common_read_csv_args)
        except UnicodeDecodeError:
            if file_encoding == 'utf-8':  # Only fallback to 'replace' if primary was utf-8
                self.logger.warning(f"Strict UTF-8 failed for {tsv_file_path}, retrying with errors='replace'")
                df = pd.read_csv(tsv_file_path, encoding='utf-8', encoding_errors='replace', **common_read_csv_args)
            else:  # If specified encoding (e.g. latin-1) failed, re-raise
                raise
        if df is None: raise ValueError("Failed to read TSV into DataFrame.")  # Should not happen if above raises

        if rename_map: df.rename(columns=rename_map, inplace=True)
        if text_replacements:
            for col_name, replacements in text_replacements.items():
                if col_name in df.columns:
                    for pattern, rep_str in replacements:
                        df[col_name] = df[col_name].astype(str).str.replace(pattern, rep_str, regex=False)
        missing_cols = [col for col in desired_columns_ordered if col not in df.columns]
        if missing_cols: raise ValueError(f"Missing columns after processing: {missing_cols}")
        df_processed = df[desired_columns_ordered]
        # End of placeholder for robust df_processed creation

        # IMPORTANT: For psycopg3's copy.write(), it expects bytes.
        # So, we write the DataFrame to a BytesIO buffer, encoded as UTF-8.
        bytes_buffer = io.BytesIO()
        # The header from df_processed.to_csv will be UTF-8 bytes.
        # The data itself (Python strings) will be encoded to UTF-8 bytes by to_csv.
        df_processed.to_csv(bytes_buffer, sep='\t', index=False, header=True, na_rep='', encoding='utf-8')
        bytes_buffer.seek(0)  # Reset buffer position to the beginning for reading

        self.logger.info(f"Successfully preprocessed '{tsv_file_path}'. BytesIO buffer ready for COPY.")
        return bytes_buffer

    def ingest_swisslipids_data(self,
                                lipids_file_path: str,
                                enzymes_file_path: str,
                                lipids2uniprot_file_path: Optional[str] = None):
        self.logger.info("Starting ingestion of SwissLipids data into staging tables with preprocessing.")

        staging_lipids_table = self.config.get_staging_lipids_table()
        staging_enzymes_table = self.config.get_staging_enzymes_table()

        lipids_staging_table_cols_ordered = [  # As defined before
            "Lipid ID", "Level", "Name", "Abbreviation*", "Synonyms*",
            "Lipid class*", "Parent", "Components*", "SMILES (pH7.3)",
            "Formula (pH7.3)", "Exact Mass (neutral form)"
        ]
        # Use config or hardcode confirmed encodings
        lipids_file_actual_encoding = 'latin-1'  # latin-1 because lipids.tsv is problematic

        enzymes_staging_table_cols_ordered = [  # As defined before
            "SwissLipids ID", "UniProtKB AC(s)", "Gene name", "Protein taxon",
            "Taxon scientific name", "Rhea ID", "Reaction text", "Evidence tag ID"
        ]
        enzymes_file_actual_encoding = 'utf-8'  # utf-8 for enzymes.tsv
        enzymes_text_replacements = {"Reaction text": [("=&gt;", "=>")]}

        try:
            self.logger.info(
                f"Processing lipids file: {lipids_file_path} with file-specific encoding: {lipids_file_actual_encoding}")
            self._load_tsv_to_staging(
                lipids_file_path,
                staging_lipids_table,
                staging_table_columns_ordered=lipids_staging_table_cols_ordered,
                file_specific_encoding=lipids_file_actual_encoding,  # Pass correct encoding
                text_replacements_in_cols=None  # No replacements for lipids.tsv for now
            )

            self.logger.info(
                f"Processing enzymes file: {enzymes_file_path} with file-specific encoding: {enzymes_file_actual_encoding}")
            self._load_tsv_to_staging(
                enzymes_file_path,
                staging_enzymes_table,
                staging_table_columns_ordered=enzymes_staging_table_cols_ordered,
                file_specific_encoding=enzymes_file_actual_encoding,  # Pass correct encoding
                text_replacements_in_cols=enzymes_text_replacements
            )
            # ... (handle lipids2uniprot) ...
            self.logger.info("All SwissLipids data successfully ingested into staging tables.")
        except Exception as e:
            self.logger.critical(f"Failed to ingest all SwissLipids data: {e}")
            raise

    def _generate_manual_entry_templates(self, run_output_dir: Path) -> None:
        """Generates blank template files for manual additions."""
        self.logger.info("Generating manual entry templates...")

        # --- Manual Reactions Template ---
        manual_reactions_template_path = run_output_dir / "template_manual_add_reactions.tsv"

        headers = [
            "InputReactionText_CleanedAbbr",
            # e.g., "Cer(d18:1/16:0) + PC(16:0/18:1) <=> SM(d18:1/16:0) + DAG(16:0/18:1)"
            "Enzyme_UniProtKB_AC",  # e.g., "P12345" or "P12345|Q67890"
            "Organism_ID",  # e.g., "9606" (NCBI Taxonomy ID)
            "Organism_ScientificName",  # e.g., "Homo sapiens" (Optional, for clarity)
            "Enzyme_GeneName",  # e.g., "SMS1" (Optional)
            "Rhea_ID",  # (Optional)
            "DOI"  # (Optional)
        ]

        # Create an empty DataFrame with these headers to ensure file creation
        pd.DataFrame(columns=headers).to_csv(manual_reactions_template_path, sep='\t', index=False, encoding='utf-8')

        self.logger.info(f"Generated manual reaction template: {manual_reactions_template_path}")
        # Add calls here to generate other templates (e.g., for manual molecules) if needed in the future.

    def _perform_diff_and_generate_reports(self,
                                           run_output_dir_path: Path,  # Path object for the output dir of this run
                                           prev_lipids_snapshot_df: Optional[pd.DataFrame] = None,
                                           prev_enzymes_snapshot_df: Optional[pd.DataFrame] = None
                                           ) -> Dict[str, Any]:  # Returns summary
        """
        Performs the core difference detection (2-way or 3-way based on snapshot availability)
        and generates all delta CSV files and a summary report.
        """
        self.logger.info("Performing difference detection and generating reports...")

        # --- Find all differences ---
        # Pass snapshot DataFrames to the find methods. They will handle if they are None.
        new_organisms_df = self._find_new_organisms(prev_enzymes_snapshot_df=prev_enzymes_snapshot_df)
        new_molecules_df = self._find_new_molecules(prev_lipids_snapshot_df=prev_lipids_snapshot_df)
        new_reactions_df = self._find_new_reaction_definitions(prev_enzymes_snapshot_df=prev_enzymes_snapshot_df)
        new_enzymes_df = self._find_new_enzyme_definitions(prev_enzymes_snapshot_df=prev_enzymes_snapshot_df)
        new_reaction_links_df = self._find_new_reaction_enzyme_links(prev_enzymes_snapshot_df)

        updated_entries_df = self._find_updated_entries(
            prev_lipids_snapshot_df=prev_lipids_snapshot_df,
            prev_enzymes_snapshot_df=prev_enzymes_snapshot_df
        )

        # --- Generate delta CSV files (using the existing _generate_single_delta_file) ---
        # Note: _generate_single_delta_file needs to be adapted to take run_output_dir_path
        # instead of constructing it from self.output_dir / timestamp internally.

        # Helper to adapt _generate_single_delta_file call for now
        def save_delta(df, filename_suffix):
            self._generate_single_delta_file(df, filename_suffix, run_output_dir_path)

        save_delta(new_organisms_df, "delta_new_organisms.tsv")
        save_delta(new_molecules_df, "delta_new_molecules.tsv")
        save_delta(new_reactions_df, "delta_new_reaction_definitions.tsv")
        save_delta(new_enzymes_df, "delta_new_enzyme_definitions.tsv")
        save_delta(new_reaction_links_df, "delta_new_reaction_links.tsv")
        save_delta(updated_entries_df, "delta_updated_entries.tsv")

        # --- Prepare and generate summary report ---
        summary = {
            'run_timestamp_dir': str(run_output_dir_path.name),  # The YYYYMMDD_HHMMSS_test part
            'new_organisms': len(new_organisms_df),
            'new_molecules': len(new_molecules_df),
            'new_reaction_definitions': len(new_reactions_df),
            'new_enzyme_definitions': len(new_enzymes_df),
            'new_reaction_links': len(new_reaction_links_df),
            'updated_entries': len(updated_entries_df)
        }
        # _generate_summary_report also needs to take run_output_dir_path
        self._generate_summary_report(summary, "update_summary.txt", run_output_dir_path)

        self.logger.info("Delta file generation complete based on current diff logic.")
        return summary

    def generate_delta_files_for_review(self,
                                        current_data_run_timestamp: Optional[str] = None,
                                        previous_snapshot_lipids_file_path: Optional[str] = None,
                                        previous_snapshot_enzymes_file_path: Optional[str] = None
                                        ) -> Dict[str, Any]:
        """
        Orchestrates finding all new/updated entities and generating CSV delta files
        and a summary report. This is the main method called by the 'generate' mode.

        Args:
            current_data_run_timestamp (Optional[str]): A string representing the "version" or
                                                timestamp for organizing output files for this run.
                                                If None, uses current datetime.
            previous_snapshot_lipids_file_path (Optional[str]): Path to previous lipids snapshot TSV.
            previous_snapshot_enzymes_file_path (Optional[str]): Path to previous enzymes snapshot TSV.
        Returns:
            Dict[str, Any]: A summary of the findings.
        """
        self.logger.info("Starting generation of delta files for review...")

        run_timestamp_str = current_data_run_timestamp if current_data_run_timestamp else datetime.now().strftime(
            "%Y%m%d_%H%M%S")
        run_output_dir = self.output_dir / run_timestamp_str  # e.g., db_updates/20250517_103000/
        run_output_dir.mkdir(parents=True, exist_ok=True)

        # 1. Generate manual entry templates
        self._generate_manual_entry_templates(run_output_dir)

        # 2. Load previous snapshots into DataFrames if files exist
        #    Define desired columns for snapshots (should match what was archived)
        lipids_snapshot_cols = [
            "Lipid ID", "Level", "Name", "Abbreviation*", "Synonyms*",
            "Lipid class*", "Parent", "Components*", "SMILES (pH7.3)",
            "Formula (pH7.3)", "Exact Mass (neutral form)"
        ]
        enzymes_snapshot_cols = [
            "SwissLipids ID", "UniProtKB AC(s)", "Gene name", "Protein taxon",
            "Taxon scientific name", "Rhea ID", "Reaction text", "Evidence tag ID"
        ]

        prev_lipids_df = None
        if previous_snapshot_lipids_file_path and Path(previous_snapshot_lipids_file_path).exists():
            self.logger.info(f"Loading previous lipids snapshot from: {previous_snapshot_lipids_file_path}")
            # Assuming snapshots are UTF-8 and clean, or use a specific encoding if known
            buffer = self._preprocess_tsv_for_copy(previous_snapshot_lipids_file_path, lipids_snapshot_cols,
                                                   file_encoding='utf-8')
            prev_lipids_df = pd.read_csv(buffer, sep='\t', dtype=str, keep_default_na=False, na_values=[''])
            self.logger.info(f"Loaded {len(prev_lipids_df)} records from previous lipids snapshot.")
        else:
            self.logger.info("No previous lipids snapshot file provided or found.")

        prev_enzymes_df = None
        if previous_snapshot_enzymes_file_path and Path(previous_snapshot_enzymes_file_path).exists():
            self.logger.info(f"Loading previous enzymes snapshot from: {previous_snapshot_enzymes_file_path}")
            buffer = self._preprocess_tsv_for_copy(previous_snapshot_enzymes_file_path, enzymes_snapshot_cols,
                                                   file_encoding='utf-8')
            prev_enzymes_df = pd.read_csv(buffer, sep='\t', dtype=str, keep_default_na=False, na_values=[''])
            self.logger.info(f"Loaded {len(prev_enzymes_df)} records from previous enzymes snapshot.")
        else:
            self.logger.info("No previous enzymes snapshot file provided or found.")

        # 3. Call the core diffing and report generation logic
        summary = self._perform_diff_and_generate_reports(
            run_output_dir_path=run_output_dir,
            prev_lipids_snapshot_df=prev_lipids_df,
            prev_enzymes_snapshot_df=prev_enzymes_df
        )

        self.logger.info("Delta file generation process completed.")
        return summary

    def _query_to_df(self, query: str, params: Optional[Tuple] = None) -> pd.DataFrame:
        """
        Executes a SQL query and returns results as a Pandas DataFrame.
        Uses the column names from the cursor description.
        """
        self.logger.debug(f"Executing query for DataFrame: {query} with params: {params}")
        try:
            self.cursor.execute(query, params or ())

            # Check if the query returned any rows before fetching description
            if self.cursor.description is None:  # For queries like INSERT/UPDATE without RETURNING
                self.logger.debug("Query did not return rows (e.g., DML without RETURNING). Returning empty DataFrame.")
                return pd.DataFrame()

            colnames = [desc[0] for desc in self.cursor.description]
            results = self.cursor.fetchall()
            df = pd.DataFrame(results, columns=colnames)
            self.logger.debug(f"Query returned {len(df)} rows.")
            return df
        except psycopg.Error as e:
            self.conn.rollback()  # Rollback if a query fails within a transaction context
            self.logger.error(f"Database query failed: {e}\nQuery: {query}\nParams: {params}")
            # Depending on context, you might want to raise the error or return an empty DataFrame
            raise  # Or return pd.DataFrame() if that's preferable for the calling code
        except Exception as e:
            self.logger.error(f"Unexpected error executing query: {e}\nQuery: {query}\nParams: {params}")
            raise

    def _dataframe_to_temp_db_table(self, df: pd.DataFrame, temp_table_name: str) -> bool:
        """
        Creates and populates a temporary PostgreSQL table from a pandas DataFrame.
        The temporary table is session-specific and will be dropped automatically
        when the database session ends.

        Args:
            df (pd.DataFrame): The DataFrame to load.
            temp_table_name (str): The desired name for the temporary table (must be SQL-safe).

        Returns:
            bool: True if the table was created and populated, False otherwise.
        """
        if df is None or df.empty:
            self.logger.info(f"DataFrame for temporary table '{temp_table_name}' is empty. Skipping creation.")
            # It's important that downstream SQL queries can handle the temp table not existing
            # (e.g., by using LEFT JOIN and checking for NULLs).
            # Alternatively, create an empty temp table with the expected schema if queries demand its presence.
            # For now, we'll just skip creation and return False.
            return False

        if self.conn is None or self.conn.closed or self.cursor is None or self.cursor.closed:
            self.logger.error(f"Database connection not available for creating temp table '{temp_table_name}'.")
            raise ConnectionError("Database connection not established for temp table creation.")

        self.logger.info(f"Loading DataFrame into temporary table: {temp_table_name} ({len(df)} rows)")

        # Ensure temp_table_name is safe (e.g., no SQL injection if it were dynamic, though here it's controlled)
        # For temporary tables, they are session-scoped, so name clashes are less of an issue
        # but good practice to ensure it's a valid identifier. We don't need IF NOT EXISTS for CREATE TEMP.

        try:
            # 1. Drop if it somehow exists from a failed previous step in the same session (unlikely but safe)
            self.cursor.execute(f"DROP TABLE IF EXISTS {temp_table_name};")

            # 2. Create the temporary table based on DataFrame columns.
            #    All columns will be TEXT for simplicity of loading and comparison.
            #    Ensure column names from DataFrame are SQL-safe (e.g., quoted if they have spaces/special chars).
            #    Pandas column names from read_csv(..., sep='\t') should be fine if from typical TSV headers.
            cols_sql_definitions = []
            for col_name in df.columns:
                # Basic quoting for safety if column names have spaces or are keywords
                safe_col_name = f'"{col_name}"'
                cols_sql_definitions.append(f'{safe_col_name} TEXT')

            create_table_sql = f"CREATE TEMPORARY TABLE {temp_table_name} ({', '.join(cols_sql_definitions)});"
            self.logger.debug(f"Creating temporary table with SQL: {create_table_sql}")
            self.cursor.execute(create_table_sql)

            # 3. Prepare data buffer (BytesIO with UTF-8 encoded TSV content)
            #    This is the same as what _preprocess_tsv_for_copy produces, but from an in-memory df.
            bytes_buffer = io.BytesIO()
            df.to_csv(bytes_buffer, sep='\t', index=False, header=True, na_rep='', encoding='utf-8')
            bytes_buffer.seek(0)

            # 4. COPY data into the temporary table
            #    The column order in the buffer (from df.columns) must match the order in CREATE TEMPORARY TABLE.
            #    Explicitly listing columns in COPY is safest.
            copy_columns_sql_list = [f'"{col}"' for col in df.columns]
            copy_sql = f"COPY {temp_table_name} ({', '.join(copy_columns_sql_list)}) FROM STDIN WITH (FORMAT CSV, DELIMITER E'\\t', HEADER TRUE, NULL '', ENCODING 'UTF8');"

            self.logger.debug(f"Executing COPY to temporary table {temp_table_name}...")
            with self.cursor.copy(copy_sql) as copy:
                chunk_size = 8192
                while True:
                    chunk = bytes_buffer.read(chunk_size)
                    if not chunk:
                        break
                    copy.write(chunk)

            # Temporary tables are typically part of the current transaction.
            # No separate commit is strictly needed *for the temp table creation itself*
            # if it's part of a larger operation that will be committed/rolled back.
            # However, for safety and clarity if this method is called standalone:
            # self.conn.commit() # Or manage transactions at a higher level.
            # For now, assume higher level manages the overall transaction.

            self.logger.info(
                f"Successfully created and populated temporary table {temp_table_name} with {len(df)} rows.")
            return True
        except psycopg.Error as e:
            self.logger.error(
                f"Database error creating/populating temporary table {temp_table_name}: {e.diag.message_primary if e.diag else e}",
                exc_info=True)
            # self.conn.rollback() # Rollback if this was part of a larger ongoing transaction
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error with temporary table {temp_table_name}: {e}", exc_info=True)
            # self.conn.rollback()
            raise

    def _find_new_organisms(self, prev_enzymes_snapshot_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        self.logger.info("Finding new organisms (3-way capable)...")
        staging_enzymes_table = self.config.get_staging_enzymes_table()

        temp_snapshot_table_name = f"temp_prev_org_snap_{datetime.now().strftime('%Y%m%d%H%M%S%f')}"
        snapshot_table_created_and_used = False  # Flag to track if temp table was made and should be dropped

        select_clause_base = """
            SELECT DISTINCT
                se."Protein taxon" AS organism_id, 
                se."Taxon scientific name" AS taxon_scientific_name
        """
        from_clause = f" FROM {staging_enzymes_table} se"
        join_production_clause = " LEFT JOIN lipograph.organism prod_o ON CAST(NULLIF(TRIM(se.\"Protein taxon\"), '') AS INTEGER) = prod_o.organism_id"

        # Common WHERE clause: must not be in production and must have valid taxon info
        where_clause_base = """
            WHERE prod_o.organism_id IS NULL 
              AND NULLIF(TRIM(se."Protein taxon"), '') IS NOT NULL 
              AND NULLIF(TRIM(se."Taxon scientific name"), '') IS NOT NULL
        """

        if prev_enzymes_snapshot_df is not None and not prev_enzymes_snapshot_df.empty:
            # We only need 'Protein taxon' from the enzymes snapshot for organisms
            relevant_snapshot_cols = ['Protein taxon']
            if all(col in prev_enzymes_snapshot_df.columns for col in relevant_snapshot_cols):
                org_snapshot_df = prev_enzymes_snapshot_df[relevant_snapshot_cols].copy()
                org_snapshot_df.rename(columns={'Protein taxon': 'snapshot_organism_id'},
                                       inplace=True)  # Avoid name clash
                org_snapshot_df.drop_duplicates(subset=['snapshot_organism_id'], inplace=True)

                if not org_snapshot_df.empty:
                    try:
                        if self._dataframe_to_temp_db_table(org_snapshot_df, temp_snapshot_table_name):
                            snapshot_table_created_and_used = True
                            self.logger.info(
                                f"Using previous enzymes snapshot (for organisms) via temp table: {temp_snapshot_table_name}.")

                            select_clause_extra = """,
                                CASE
                                    WHEN prev_snap.snapshot_organism_id IS NULL THEN 'NewToSwissLipids_And_Production'
                                    ELSE 'ReappearedInSL_Or_NewToProduction' 
                                END AS novelty_status,
                                prev_snap.snapshot_organism_id IS NULL AS is_new_since_last_sl_snapshot
                            """
                            join_snapshot_clause = f" LEFT JOIN {temp_snapshot_table_name} prev_snap ON se.\"Protein taxon\" = prev_snap.snapshot_organism_id"
                            query = select_clause_base + select_clause_extra + from_clause + join_production_clause + join_snapshot_clause + where_clause_base
                        else:  # Temp table creation failed
                            self.logger.warning(
                                "Failed to create temp table for organism snapshot. Falling back to 2-way diff.")
                            select_clause_extra = """, 'NewToProduction_SnapshotError' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                            query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
                    except Exception as e:
                        self.logger.error(
                            f"Error processing organism snapshot for 3-way diff: {e}. Falling back to 2-way.")
                        select_clause_extra = """, 'NewToProduction_SnapshotProcessingError' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                        query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
                else:  # Snapshot relevant columns were empty
                    self.logger.info("Organism snapshot data was empty after filtering. Performing 2-way diff.")
                    select_clause_extra = """, 'NewToProduction_SnapshotEmpty' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                    query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
            else:  # Snapshot DF didn't have the required columns
                self.logger.warning(
                    "Previous enzymes snapshot DataFrame missing 'Protein taxon'. Performing 2-way diff for organisms.")
                select_clause_extra = """, 'NewToProduction_SnapshotColMissing' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
        else:  # No snapshot DataFrame provided
            self.logger.info("No previous enzymes snapshot provided. Performing 2-way diff for new organisms.")
            select_clause_extra = """, 'NewToProduction_SnapshotUnavailable' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
            query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base

        df = self._query_to_df(query)

        if snapshot_table_created_and_used:
            try:
                self.cursor.execute(f"DROP TABLE IF EXISTS {temp_snapshot_table_name};")
                self.logger.info(f"Dropped temporary snapshot table: {temp_snapshot_table_name}")
            except psycopg.Error as e:
                self.logger.warning(f"Could not drop temporary table {temp_snapshot_table_name}: {e}")

        self.logger.info(f"Found {len(df)} potential new organisms.")
        return df

    def _find_new_molecules(self, prev_lipids_snapshot_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        self.logger.info("Finding new molecules (3-way capable)...")
        staging_lipids_table = self.config.get_staging_lipids_table()

        temp_snapshot_table_name = f"temp_prev_mol_snap_{datetime.now().strftime('%Y%m%d%H%M%S%f')}"
        snapshot_table_created_and_used = False

        select_clause_base = """
            SELECT
                sl."Name" AS molecule_name,
                sl."Lipid ID" AS swisslipids_id,
                sl."Abbreviation*" AS abbreviation,
                sl."Abbreviation*" AS cleaned_abbreviation, -- Initially same as raw, will be processed in Python
                sl."Synonyms*" AS molecule_synonyms,
                sl."Lipid class*" AS sl_lipid_class
        """
        from_clause = f" FROM {staging_lipids_table} sl"
        join_production_clause = " LEFT JOIN lipograph.molecules prod_m ON sl.\"Lipid ID\" = prod_m.swisslipids_id"
        where_clause_base = " WHERE prod_m.swisslipids_id IS NULL AND sl.\"Lipid ID\" IS NOT NULL"

        if prev_lipids_snapshot_df is not None and not prev_lipids_snapshot_df.empty:
            # Ensure the snapshot has the key column 'Lipid ID'
            if '"Lipid ID"' not in prev_lipids_snapshot_df.columns and 'Lipid ID' not in prev_lipids_snapshot_df.columns:
                self.logger.warning(
                    "Previous lipids snapshot DataFrame missing 'Lipid ID' key column. Performing 2-way diff.")
                select_clause_extra = """, 'NewToProduction_SnapshotKeyColMissing' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
            else:
                # If column name has quotes from df creation, handle it. Assume it's 'Lipid ID' if no quotes.
                snapshot_key_col = '"Lipid ID"' if '"Lipid ID"' in prev_lipids_snapshot_df.columns else 'Lipid ID'
                # Create a minimal snapshot df for joining
                mol_snapshot_df = prev_lipids_snapshot_df[[snapshot_key_col]].copy()
                mol_snapshot_df.rename(columns={snapshot_key_col: 'snapshot_lipid_id'}, inplace=True)
                mol_snapshot_df.drop_duplicates(subset=['snapshot_lipid_id'], inplace=True)

                if not mol_snapshot_df.empty:
                    try:
                        if self._dataframe_to_temp_db_table(mol_snapshot_df, temp_snapshot_table_name):
                            snapshot_table_created_and_used = True
                            self.logger.info(
                                f"Using previous lipids snapshot via temp table: {temp_snapshot_table_name} for 3-way diff.")
                            select_clause_extra = """,
                                CASE
                                    WHEN prev_snap.snapshot_lipid_id IS NULL THEN 'NewToSwissLipids_And_Production'
                                    ELSE 'ReappearedInSL_Or_NewToProduction'
                                END AS novelty_status,
                                prev_snap.snapshot_lipid_id IS NULL AS is_new_since_last_sl_snapshot
                            """
                            join_snapshot_clause = f" LEFT JOIN {temp_snapshot_table_name} prev_snap ON sl.\"Lipid ID\" = prev_snap.snapshot_lipid_id"
                            query = select_clause_base + select_clause_extra + from_clause + join_production_clause + join_snapshot_clause + where_clause_base
                        else:  # Temp table creation failed
                            self.logger.warning(
                                "Failed to create temp table for lipids snapshot. Falling back to 2-way diff.")
                            select_clause_extra = """, 'NewToProduction_SnapshotError' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                            query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
                    except Exception as e:
                        self.logger.error(
                            f"Error processing lipids snapshot for 3-way diff: {e}. Falling back to 2-way.")
                        select_clause_extra = """, 'NewToProduction_SnapshotProcessingError' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                        query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
                else:  # Snapshot relevant columns were empty
                    self.logger.info("Lipids snapshot data was empty after filtering. Performing 2-way diff.")
                    select_clause_extra = """, 'NewToProduction_SnapshotEmpty' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                    query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
        else:  # No snapshot DataFrame provided
            self.logger.info("No previous lipids snapshot. Performing 2-way diff for new molecules.")
            select_clause_extra = """, 'NewToProduction_SnapshotUnavailable' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
            query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base

        df = self._query_to_df(query)

        if snapshot_table_created_and_used:
            try:
                self.cursor.execute(f"DROP TABLE IF EXISTS {temp_snapshot_table_name};")
                self.logger.info(f"Dropped temporary snapshot table: {temp_snapshot_table_name}")
            except psycopg.Error as e:
                self.logger.warning(f"Could not drop temporary table {temp_snapshot_table_name}: {e}")

        # Post-process 'cleaned_abbreviation' (as before, using your LipidTranslator logic)
        if not df.empty:
            self.logger.debug(f"Applying Python-side cleaning to 'abbreviation' for {len(df)} new molecules.")
            # The 'abbreviation' column contains the raw abbreviation from SwissLipids.
            # The 'cleaned_abbreviation' column was initially a copy, now we overwrite it with the cleaned version.
            df['cleaned_abbreviation'] = df['abbreviation'].apply(
                lambda x: self.lipid_translator.clean_lipid_abbreviation(x) if pd.notna(x) and str(
                    x).strip() != '' else None
            )

        self.logger.info(f"Found {len(df)} potential new molecules.")
        return df

    def _find_new_reaction_definitions(self, prev_enzymes_snapshot_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        self.logger.info("Finding new reaction definitions (3-way capable)...")
        staging_enzymes_table = self.config.get_staging_enzymes_table()

        temp_snapshot_table_name = f"temp_prev_rxn_snap_{datetime.now().strftime('%Y%m%d%H%M%S%f')}"
        snapshot_table_created_and_used = False

        select_clause_base = """
            SELECT DISTINCT
                se."Reaction text" AS reaction_text,
                CAST(NULLIF(TRIM(se."Rhea ID"), '') AS INTEGER) AS rhea_id,
                NULL AS doi
        """
        from_clause = f" FROM {staging_enzymes_table} se"
        join_production_clause = """
            LEFT JOIN lipograph.reactions prod_r
                ON se."Reaction text" = prod_r.reaction_text
                AND (CAST(NULLIF(TRIM(se."Rhea ID"), '') AS INTEGER) IS NOT DISTINCT FROM prod_r.rhea_id)
        """  # Key for reaction is (text, rhea_id)
        where_clause_base = " WHERE prod_r.reaction_id IS NULL AND se.\"Reaction text\" IS NOT NULL"

        if prev_enzymes_snapshot_df is not None and not prev_enzymes_snapshot_df.empty:
            # Need "Reaction text" and "Rhea ID" from snapshot
            snapshot_key_cols = ['Reaction text', 'Rhea ID']  # Assuming these are the names in the snapshot df
            if all(col in prev_enzymes_snapshot_df.columns for col in snapshot_key_cols):
                rxn_snapshot_df = prev_enzymes_snapshot_df[snapshot_key_cols].copy()
                rxn_snapshot_df.rename(
                    columns={'Reaction text': 'snapshot_reaction_text', 'Rhea ID': 'snapshot_rhea_id'},
                    inplace=True)
                rxn_snapshot_df['snapshot_rhea_id'] = pd.to_numeric(
                    rxn_snapshot_df['snapshot_rhea_id'].str.strip().replace('', None), errors='coerce').astype(
                    'Int64')  # Handle NULLs and convert
                rxn_snapshot_df.drop_duplicates(subset=['snapshot_reaction_text', 'snapshot_rhea_id'], inplace=True)

                if not rxn_snapshot_df.empty:
                    try:
                        if self._dataframe_to_temp_db_table(rxn_snapshot_df, temp_snapshot_table_name):
                            snapshot_table_created_and_used = True
                            self.logger.info(
                                f"Using previous enzymes snapshot (for reactions) via temp table: {temp_snapshot_table_name}.")
                            select_clause_extra = """,
                                CASE
                                    WHEN prev_snap.snapshot_reaction_text IS NULL THEN 'NewToSwissLipids_And_Production'
                                    ELSE 'ReappearedInSL_Or_NewToProduction'
                                END AS novelty_status,
                                prev_snap.snapshot_reaction_text IS NULL AS is_new_since_last_sl_snapshot
                            """
                            join_snapshot_clause = f"""
                                LEFT JOIN {temp_snapshot_table_name} prev_snap
                                    ON se."Reaction text" = prev_snap.snapshot_reaction_text
                                    AND (CAST(NULLIF(TRIM(se."Rhea ID"), '') AS INTEGER) IS NOT DISTINCT FROM CAST(NULLIF(TRIM(prev_snap.snapshot_rhea_id), '') AS INTEGER))
                            """
                            query = select_clause_base + select_clause_extra + from_clause + join_production_clause + join_snapshot_clause + where_clause_base
                        else:  # Temp table creation failed
                            self.logger.warning(
                                "Failed to create temp table for reaction snapshot. Falling back to 2-way diff.")
                            select_clause_extra = """, 'NewToProduction_SnapshotError' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                            query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
                    except Exception as e:
                        self.logger.error(
                            f"Error processing reaction snapshot for 3-way diff: {e}. Falling back to 2-way.")
                        select_clause_extra = """, 'NewToProduction_SnapshotProcessingError' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                        query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
                else:  # Snapshot relevant columns were empty
                    self.logger.info("Reaction snapshot data was empty after filtering. Performing 2-way diff.")
                    select_clause_extra = """, 'NewToProduction_SnapshotEmpty' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                    query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
            else:  # Snapshot DF didn't have the required columns
                self.logger.warning(
                    "Previous enzymes snapshot DataFrame missing reaction key columns. Performing 2-way diff for reactions.")
                select_clause_extra = """, 'NewToProduction_SnapshotColMissing' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
        else:  # No snapshot DataFrame provided
            self.logger.info("No previous enzymes snapshot. Performing 2-way diff for new reaction definitions.")
            select_clause_extra = """, 'NewToProduction_SnapshotUnavailable' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
            query = select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base

        df = self._query_to_df(query)

        if snapshot_table_created_and_used:
            try:
                self.cursor.execute(f"DROP TABLE IF EXISTS {temp_snapshot_table_name};")
                self.logger.info(f"Dropped temporary snapshot table: {temp_snapshot_table_name}")
            except psycopg.Error as e:
                self.logger.warning(f"Could not drop temporary table {temp_snapshot_table_name}: {e}")

        self.logger.info(f"Found {len(df)} potential new reaction definitions.")
        return df

    def _find_new_enzyme_definitions(self, prev_enzymes_snapshot_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        self.logger.info("Finding new enzyme definitions (3-way capable)...")
        staging_enzymes_table = self.config.get_staging_enzymes_table()

        temp_snapshot_table_name = f"temp_prev_enz_snap_{datetime.now().strftime('%Y%m%d%H%M%S%f')}"
        snapshot_table_created_and_used = False

        # Base query for unnesting current staging enzymes
        unnested_staging_cte = f"""
            WITH UnnestedStagingEnzymes AS (
                SELECT
                    se."Gene name" AS enzyme_name,
                    TRIM(unnest(string_to_array(se."UniProtKB AC(s)", '|'))) AS uniprot_id,
                    CAST(NULLIF(TRIM(se."Protein taxon"), '') AS INTEGER) AS organism_id,
                    se."SwissLipids ID" AS swisslipids_p_id -- This is the SLP: ID from enzymes.tsv
                FROM {staging_enzymes_table} se
                WHERE NULLIF(TRIM(se."UniProtKB AC(s)"), '') IS NOT NULL
            )
        """
        select_clause_base = """
            SELECT DISTINCT
                use.enzyme_name,
                use.uniprot_id,
                use.organism_id,
                use.swisslipids_p_id
        """
        from_clause = " FROM UnnestedStagingEnzymes use"
        join_production_clause = " LEFT JOIN lipograph.enzymes prod_e ON use.uniprot_id = prod_e.uniprot_id"
        where_clause_base = " WHERE prod_e.enzyme_id IS NULL AND use.uniprot_id IS NOT NULL"  # New if not in prod by UniProt ID

        if prev_enzymes_snapshot_df is not None and not prev_enzymes_snapshot_df.empty:
            # Need "UniProtKB AC(s)" from snapshot, then unnest it similarly
            if '"UniProtKB AC(s)"' not in prev_enzymes_snapshot_df.columns and 'UniProtKB AC(s)' not in prev_enzymes_snapshot_df.columns:
                self.logger.warning(
                    "Previous enzymes snapshot DataFrame missing 'UniProtKB AC(s)'. Performing 2-way diff for enzymes.")
                select_clause_extra = """, 'NewToProduction_SnapshotKeyColMissing' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                query = unnested_staging_cte + select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
            else:
                # Create a DataFrame of unique UniProt IDs from the snapshot
                snapshot_uniprot_ids = set()
                key_col_name = '"UniProtKB AC(s)"' if '"UniProtKB AC(s)"' in prev_enzymes_snapshot_df.columns else 'UniProtKB AC(s)'
                for ac_string in prev_enzymes_snapshot_df[key_col_name].dropna():
                    for ac in ac_string.split('|'):
                        snapshot_uniprot_ids.add(ac.strip())

                if snapshot_uniprot_ids:
                    enz_snapshot_df = pd.DataFrame(list(snapshot_uniprot_ids), columns=['snapshot_uniprot_id'])
                    try:
                        if self._dataframe_to_temp_db_table(enz_snapshot_df, temp_snapshot_table_name):
                            snapshot_table_created_and_used = True
                            self.logger.info(
                                f"Using previous enzymes snapshot (for UniProt IDs) via temp table: {temp_snapshot_table_name}.")
                            select_clause_extra = """,
                                CASE
                                    WHEN prev_snap.snapshot_uniprot_id IS NULL THEN 'NewToSwissLipids_And_Production'
                                    ELSE 'ReappearedInSL_Or_NewToProduction'
                                END AS novelty_status,
                                prev_snap.snapshot_uniprot_id IS NULL AS is_new_since_last_sl_snapshot
                            """
                            join_snapshot_clause = f" LEFT JOIN {temp_snapshot_table_name} prev_snap ON use.uniprot_id = prev_snap.snapshot_uniprot_id"
                            query = unnested_staging_cte + select_clause_base + select_clause_extra + from_clause + join_production_clause + join_snapshot_clause + where_clause_base
                        else:  # Temp table creation failed
                            self.logger.warning(
                                "Failed to create temp table for enzyme snapshot. Falling back to 2-way diff.")
                            select_clause_extra = """, 'NewToProduction_SnapshotError' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                            query = unnested_staging_cte + select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
                    except Exception as e:
                        self.logger.error(
                            f"Error processing enzyme snapshot for 3-way diff: {e}. Falling back to 2-way.")
                        select_clause_extra = """, 'NewToProduction_SnapshotProcessingError' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                        query = unnested_staging_cte + select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
                else:  # Snapshot relevant columns were empty
                    self.logger.info(
                        "Enzyme snapshot data (UniProt IDs) was empty after processing. Performing 2-way diff.")
                    select_clause_extra = """, 'NewToProduction_SnapshotEmpty' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
                    query = unnested_staging_cte + select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base
        else:  # No snapshot DataFrame provided
            self.logger.info("No previous enzymes snapshot. Performing 2-way diff for new enzyme definitions.")
            select_clause_extra = """, 'NewToProduction_SnapshotUnavailable' AS novelty_status, TRUE AS is_new_since_last_sl_snapshot """
            query = unnested_staging_cte + select_clause_base + select_clause_extra + from_clause + join_production_clause + where_clause_base

        df = self._query_to_df(query)

        if snapshot_table_created_and_used:
            try:
                self.cursor.execute(f"DROP TABLE IF EXISTS {temp_snapshot_table_name};")
                self.logger.info(f"Dropped temporary snapshot table: {temp_snapshot_table_name}")
            except psycopg.Error as e:
                self.logger.warning(f"Could not drop temporary table {temp_snapshot_table_name}: {e}")

        self.logger.info(f"Found {len(df)} potential new enzyme definitions.")
        return df

    def _find_new_reaction_enzyme_links(self, prev_enzymes_snapshot_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        Finds new reaction-enzyme associations by comparing the current staging table
        to the previous snapshot. Returns a DataFrame of the full rows from staging
        that represent new links.
        """
        self.logger.info("Finding new reaction-enzyme links...")
        staging_enzymes_table = self.config.get_staging_enzymes_table()

        # Load the current staging data into a DataFrame
        current_staging_df = self._query_to_df(f"SELECT * FROM {staging_enzymes_table};")

        # If there's no snapshot, then ALL links are considered "new" for processing.
        if prev_enzymes_snapshot_df is None or prev_enzymes_snapshot_df.empty:
            self.logger.warning(
                "No previous enzymes snapshot found. All staging entries will be considered for linking.")
            # Add a column for clarity in the delta file
            current_staging_df['link_novelty_status'] = 'NewOrUnchanged_NoSnapshot'
            return current_staging_df

        self.logger.info("Comparing current staging enzymes to snapshot to find new links...")
        try:
            # To robustly find rows in A that are not in B, we use a left-only merge.
            # We must treat all columns as strings to avoid dtype issues with merge.
            # Key for a "link" is the entire row.
            cols = current_staging_df.columns.tolist()
            merged_df = pd.merge(
                current_staging_df.astype(str),
                prev_enzymes_snapshot_df.astype(str),
                on=cols,
                how='left',
                indicator=True
            )

            # The delta is the set of rows that only exist in the "left" DataFrame (current staging)
            delta_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])

            self.logger.info(f"Found {len(delta_df)} new or changed reaction-enzyme link rows.")
            delta_df['link_novelty_status'] = 'NewSinceLastSnapshot'
            return delta_df

        except Exception as e:
            self.logger.error(
                f"Failed to diff snapshots for links. Defaulting to processing all staging entries. Error: {e}")
            current_staging_df['link_novelty_status'] = 'NewOrUnchanged_DiffError'
            return current_staging_df

    def _find_updated_entries(self,
                              prev_lipids_snapshot_df: Optional[pd.DataFrame] = None,
                              prev_enzymes_snapshot_df: Optional[pd.DataFrame] = None
                              ) -> pd.DataFrame:
        self.logger.info("Finding updated entries (3-way comparison)...")
        updated_entries_list = []

        # --- Staging Table Names ---
        staging_lipids_table = self.config.get_staging_lipids_table()
        staging_enzymes_table = self.config.get_staging_enzymes_table()

        # --- Handle Lipids Snapshot ---
        temp_prev_lipids_table_name = f"temp_prev_lipids_upd_{datetime.now().strftime('%Y%m%d%H%M%S%f')}"
        prev_lipids_table_created_and_used = False
        if prev_lipids_snapshot_df is not None and not prev_lipids_snapshot_df.empty:
            required_lipid_cols = ["Lipid ID", "Name", "Abbreviation*", "Lipid class*"]
            if all(col in prev_lipids_snapshot_df.columns for col in required_lipid_cols):
                try:
                    # Pass only relevant columns to avoid loading huge snapshot into temp table if not needed
                    if self._dataframe_to_temp_db_table(prev_lipids_snapshot_df[required_lipid_cols],
                                                        temp_prev_lipids_table_name):
                        prev_lipids_table_created_and_used = True
                        self.logger.info(
                            f"Using previous lipids snapshot via temp table: {temp_prev_lipids_table_name} for updated entries.")
                except Exception as e:
                    self.logger.error(
                        f"Failed to create temp table for previous lipids snapshot (updates): {e}. Lipid updates will use 2-way logic.")
            else:
                self.logger.warning(
                    "Previous lipids snapshot DataFrame missing one or more required columns for comparison. Lipid updates will use 2-way logic.")

        # --- Molecule Updates ---
        molecule_fields_to_compare = {
            "molecule_name": {"sl_col": 'sl."Name"', "prod_col": 'prod_m.molecule_name', "snapshot_df_col": 'Name'},
            "abbreviation": {"sl_col": 'sl."Abbreviation*"', "prod_col": 'prod_m.abbreviation',
                             "snapshot_df_col": 'Abbreviation*'},
            "sl_lipid_class": {"sl_col": 'sl."Lipid class*"', "prod_col": 'prod_m.sl_lipid_class',
                               "snapshot_df_col": 'Lipid class*'}
        }
        # The _dataframe_to_temp_db_table quotes column names if created from df.columns.
        # So, prev_sl_col needs to reference "Original DataFrame Column Name".

        for field_key, field_config in molecule_fields_to_compare.items():
            sl_col = field_config["sl_col"]
            prod_col = field_config["prod_col"]
            # The snapshot_df_col is the plain column name from the pandas DataFrame
            snapshot_df_col_name = field_config["snapshot_df_col"]

            where_field_differs = f"{sl_col} IS DISTINCT FROM {prod_col}"

            # This logic is now clean, simple, and correct. No more parsing/stripping.
            prev_sl_col_sql = f'prev_sl."{snapshot_df_col_name}"' if prev_lipids_table_created_and_used else "NULL"
            select_prev_sl_val_sql = prev_sl_col_sql

            # --- REVISED LOGIC ---
            manual_override_case_sql = "FALSE"
            if prev_lipids_table_created_and_used:
                # A manual override exists if production differs from the snapshot,
                # AND the new SL data has NOT changed since the snapshot.
                manual_override_case_sql = f"""
                                CASE WHEN ({prod_col} IS DISTINCT FROM {prev_sl_col_sql}) AND ({sl_col} IS NOT DISTINCT FROM {prev_sl_col_sql})
                                     THEN TRUE ELSE FALSE
                                END
                            """

            # A change in SL exists if the new SL data is different from the snapshot.
            change_in_sl_case_sql = f"({sl_col} IS DISTINCT FROM {prev_sl_col_sql})"
            # --- END REVISED LOGIC ---

            join_prev_lipids_sql = ""
            if prev_lipids_table_created_and_used:
                join_prev_lipids_sql = f'LEFT JOIN {temp_prev_lipids_table_name} prev_sl ON sl."Lipid ID" = prev_sl."Lipid ID"'

            query_mol_field_update = f"""
                        SELECT
                            prod_m.swisslipids_id AS entity_key,
                            'molecule' AS entity_type,
                            '{field_key}' AS field_name,
                            {prod_col} AS old_value_in_production,
                            {sl_col} AS new_value_from_current_sl,
                            {select_prev_sl_val_sql} AS previous_value_from_snapshot_sl,
                            {manual_override_case_sql} AS is_manual_override_in_production,
                            {change_in_sl_case_sql} AS is_change_in_sl_since_snapshot
                        FROM {staging_lipids_table} sl
                        JOIN lipograph.molecules prod_m ON sl."Lipid ID" = prod_m.swisslipids_id
                        {join_prev_lipids_sql}
                        WHERE {where_field_differs};
                    """
            df_field_updates = self._query_to_df(query_mol_field_update)
            if not df_field_updates.empty:
                updated_entries_list.append(df_field_updates)

        if prev_lipids_table_created_and_used:
            try:
                self.cursor.execute(f"DROP TABLE IF EXISTS {temp_prev_lipids_table_name};")
                self.logger.info(f"Dropped temporary lipids snapshot table: {temp_prev_lipids_table_name}")
            except psycopg.Error as e:
                self.logger.warning(f"Could not drop {temp_prev_lipids_table_name}: {e}")

        # --- Enzyme Updates ---
        temp_prev_enzymes_table_name = f"temp_prev_enz_upd_{datetime.now().strftime('%Y%m%d%H%M%S%f')}"
        prev_enzymes_table_created_and_used = False
        if prev_enzymes_snapshot_df is not None and not prev_enzymes_snapshot_df.empty:
            required_enzyme_cols = ["UniProtKB AC(s)", "Gene name", "Protein taxon", "SwissLipids ID"]
            if all(col in prev_enzymes_snapshot_df.columns for col in required_enzyme_cols):
                try:
                    if self._dataframe_to_temp_db_table(prev_enzymes_snapshot_df, temp_prev_enzymes_table_name):
                        prev_enzymes_table_created_and_used = True
                        self.logger.info(
                            f"Using previous enzymes snapshot via temp table: {temp_prev_enzymes_table_name} for updated entries.")
                except Exception as e:
                    self.logger.error(
                        f"Failed to create temp table for previous enzymes snapshot (updates): {e}. Enzyme updates will use 2-way logic.")
            else:
                self.logger.warning(
                    "Previous enzymes snapshot DataFrame missing one or more required columns for comparison. Enzyme updates will use 2-way logic.")

        enzyme_fields_to_compare = {
            "enzyme_name": ('use.staging_gene_name', 'prod_e.enzyme_name', 'prev_use.snapshot_gene_name'),
            "organism_id": ('use.staging_organism_id', 'prod_e.organism_id', 'prev_use.snapshot_organism_id'),
            "swisslipids_p_id": (
            'use.staging_swisslipids_p_id', 'prod_e.swisslipids_p_id', 'prev_use.snapshot_swisslipids_p_id')
        }

        unnested_staging_enz_cte = f"""
            UnnestedStagingEnzymes AS (
                SELECT DISTINCT ON (uniprot_id)
                    uniprot_id,
                    staging_gene_name,
                    staging_organism_id,
                    staging_swisslipids_p_id
                FROM (
                    SELECT
                        TRIM(unnest(string_to_array(se."UniProtKB AC(s)", '|'))) AS uniprot_id,
                        se."Gene name" AS staging_gene_name,
                        CAST(NULLIF(TRIM(se."Protein taxon"), '') AS INTEGER) AS staging_organism_id,
                        se."SwissLipids ID" AS staging_swisslipids_p_id
                    FROM {staging_enzymes_table} se WHERE NULLIF(TRIM(se."UniProtKB AC(s)"), '') IS NOT NULL
                ) AS subquery
                ORDER BY uniprot_id
            )
        """
        unnested_prev_enz_cte_sql = ""
        join_prev_enzymes_sql = ""
        if prev_enzymes_table_created_and_used:
            unnested_prev_enz_cte_sql = f""",
            UnnestedPrevSnapshotEnzymes AS (
                SELECT DISTINCT ON (uniprot_id)
                    uniprot_id,
                    snapshot_gene_name,
                    snapshot_organism_id,
                    snapshot_swisslipids_p_id
                FROM (
                    SELECT
                        TRIM(unnest(string_to_array(prev_se."UniProtKB AC(s)", '|'))) AS uniprot_id,
                        prev_se."Gene name" AS snapshot_gene_name,
                        CAST(NULLIF(TRIM(prev_se."Protein taxon"), '') AS INTEGER) AS snapshot_organism_id,
                        prev_se."SwissLipids ID" AS snapshot_swisslipids_p_id
                    FROM {temp_prev_enzymes_table_name} prev_se WHERE NULLIF(TRIM(prev_se."UniProtKB AC(s)"), '') IS NOT NULL
                ) AS subquery
                ORDER BY uniprot_id
            )
            """
            join_prev_enzymes_sql = "LEFT JOIN UnnestedPrevSnapshotEnzymes prev_use ON use.uniprot_id = prev_use.uniprot_id"

        for field_key, (sl_col_expr, prod_col, prev_sl_col_expr_base) in enzyme_fields_to_compare.items():
            # ... (same logic as before for constructing enzyme field update SQL with CASTs) ...
            where_field_differs = f"{sl_col_expr} IS DISTINCT FROM {prod_col}"
            select_prev_sl_val_sql = prev_sl_col_expr_base if prev_enzymes_table_created_and_used else "NULL"
            prod_col_display = f"CAST({prod_col} AS TEXT)" if field_key == "organism_id" else prod_col
            sl_col_display = f"CAST({sl_col_expr} AS TEXT)" if field_key == "organism_id" else sl_col_expr
            select_prev_sl_val_display = f"CAST({select_prev_sl_val_sql} AS TEXT)" if field_key == "organism_id" and prev_enzymes_table_created_and_used else select_prev_sl_val_sql
            # --- REVISED LOGIC ---
            manual_override_case_sql = "FALSE"
            if prev_enzymes_table_created_and_used:
                prod_val_for_override = prod_col_display
                prev_val_for_override = f"CAST({prev_sl_col_expr_base} AS TEXT)" if field_key == "organism_id" else prev_sl_col_expr_base
                sl_val_for_override = sl_col_display

                manual_override_case_sql = f"""
                                CASE WHEN ({prod_val_for_override} IS DISTINCT FROM {prev_val_for_override}) AND ({sl_val_for_override} IS NOT DISTINCT FROM {prev_val_for_override})
                                     THEN TRUE ELSE FALSE
                                END
                            """

            sl_val_for_sl_change = sl_col_display
            prev_val_for_sl_change = f"CAST({prev_sl_col_expr_base} AS TEXT)" if field_key == "organism_id" and prev_enzymes_table_created_and_used else select_prev_sl_val_sql

            change_in_sl_case_sql = f"({sl_val_for_sl_change} IS DISTINCT FROM {prev_val_for_sl_change})"
            # --- END REVISED LOGIC ---

            query_enz_field_update = f"""
                WITH {unnested_staging_enz_cte}
                {unnested_prev_enz_cte_sql}
                SELECT
                    prod_e.uniprot_id AS entity_key, 'enzyme' AS entity_type, '{field_key}' AS field_name,
                    {prod_col_display} AS old_value_in_production, {sl_col_display} AS new_value_from_current_sl,
                    {select_prev_sl_val_display} AS previous_value_from_snapshot_sl,
                    {manual_override_case_sql} AS is_manual_override_in_production,
                    {change_in_sl_case_sql} AS is_change_in_sl_since_snapshot
                FROM UnnestedStagingEnzymes use
                JOIN lipograph.enzymes prod_e ON use.uniprot_id = prod_e.uniprot_id
                {join_prev_enzymes_sql}
                WHERE {where_field_differs};
            """
            df_field_updates = self._query_to_df(query_enz_field_update)
            if not df_field_updates.empty: updated_entries_list.append(df_field_updates)

        if prev_enzymes_table_created_and_used:
            try:
                self.cursor.execute(f"DROP TABLE IF EXISTS {temp_prev_enzymes_table_name};")
                self.logger.info(f"Dropped temporary enzymes snapshot table: {temp_prev_enzymes_table_name}")
            except psycopg.Error as e:
                self.logger.warning(f"Could not drop {temp_prev_enzymes_table_name}: {e}")

        # --- Reaction Definition Updates ---
        # Key for reaction: (reaction_text, rhea_id). We are primarily interested if Rhea ID changes for a given text.
        # Significant text changes usually mean a new reaction, handled by _find_new_reaction_definitions.
        temp_prev_rxn_table_name = f"temp_prev_rxn_upd_{datetime.now().strftime('%Y%m%d%H%M%S%f')}"
        prev_rxn_table_created_and_used = False
        if prev_enzymes_snapshot_df is not None and not prev_enzymes_snapshot_df.empty:
            # Snapshot needs "Reaction text" and "Rhea ID"
            required_rxn_cols = ["Reaction text", "Rhea ID"]
            if all(col in prev_enzymes_snapshot_df.columns for col in required_rxn_cols):
                try:
                    # Select only relevant columns for the reaction snapshot temp table
                    rxn_snapshot_subset_df = prev_enzymes_snapshot_df[required_rxn_cols].copy()
                    rxn_snapshot_subset_df.drop_duplicates(inplace=True)  # Avoid redundant rows in temp table
                    if self._dataframe_to_temp_db_table(rxn_snapshot_subset_df, temp_prev_rxn_table_name):
                        prev_rxn_table_created_and_used = True
                        self.logger.info(
                            f"Using previous reaction snapshot data via temp table: {temp_prev_rxn_table_name}.")
                except Exception as e:
                    self.logger.error(
                        f"Failed to create temp table for previous reaction snapshot (updates): {e}. Reaction updates will use 2-way logic.")
            else:
                self.logger.warning(
                    "Previous enzymes snapshot DataFrame missing reaction key columns for comparison. Reaction updates will use 2-way logic.")

        # Field to compare: rhea_id. Key: reaction_text.
        # Staging Rhea ID
        sl_rhea_col_expr = "CAST(NULLIF(TRIM(se.\"Rhea ID\"), '') AS INTEGER)"
        # Production Rhea ID
        prod_rhea_col_expr = "prod_r.rhea_id"
        # Previous Snapshot Rhea ID (from temp_prev_rxn_table_name, which has column "Rhea ID")
        # The _dataframe_to_temp_db_table creates columns like prev_rxn."Rhea ID"
        prev_sl_rhea_col_expr = 'CAST(NULLIF(TRIM(prev_rxn."Rhea ID"), \'\') AS INTEGER)' if prev_rxn_table_created_and_used else "NULL"

        where_rhea_differs = f"{sl_rhea_col_expr} IS DISTINCT FROM {prod_rhea_col_expr}"

        # For SELECT list, we want the text representation or NULL
        select_prev_sl_rhea_val_sql = 'TRIM(prev_rxn."Rhea ID")' if prev_rxn_table_created_and_used else "NULL"

        manual_override_rhea_sql = "FALSE"
        if prev_rxn_table_created_and_used:
            manual_override_rhea_sql = f"CASE WHEN {prev_sl_rhea_col_expr} IS NOT NULL AND {prod_rhea_col_expr} IS DISTINCT FROM {prev_sl_rhea_col_expr} THEN TRUE ELSE FALSE END"

        change_in_sl_rhea_sql = "TRUE"
        if prev_rxn_table_created_and_used:
            change_in_sl_rhea_sql = f"CASE WHEN {prev_sl_rhea_col_expr} IS NULL OR {sl_rhea_col_expr} IS DISTINCT FROM {prev_sl_rhea_col_expr} THEN TRUE ELSE FALSE END"

        join_prev_rxn_sql = ""
        if prev_rxn_table_created_and_used:
            join_prev_rxn_sql = f"""
                LEFT JOIN {temp_prev_rxn_table_name} prev_rxn 
                    ON se."Reaction text" = prev_rxn."Reaction text" 
                    -- Optional: AND (CAST(NULLIF(TRIM(se."Rhea ID"), '') AS INTEGER) IS NOT DISTINCT FROM CAST(NULLIF(TRIM(prev_rxn."Rhea ID"), '') AS INTEGER))
                    -- The above commented line would mean we only find updates if the reaction *key* (text+rhea) was the same previously.
                    -- The current simpler join on text means we find Rhea ID changes even if the previous Rhea ID was different.
            """

        query_rxn_field_update = f"""
            SELECT
                prod_r.reaction_id AS entity_key, 
                'reaction' AS entity_type,
                'rhea_id' AS field_name,
                CAST(prod_r.rhea_id AS TEXT) AS old_value_in_production,
                TRIM(se."Rhea ID") AS new_value_from_current_sl, 
                {select_prev_sl_rhea_val_sql} AS previous_value_from_snapshot_sl,
                {manual_override_rhea_sql} AS is_manual_override_in_production,
                {change_in_sl_rhea_sql} AS is_change_in_sl_since_snapshot
            FROM {staging_enzymes_table} se
            JOIN lipograph.reactions prod_r ON se."Reaction text" = prod_r.reaction_text
            {join_prev_rxn_sql}
            WHERE {where_rhea_differs}; 
        """
        df_rxn_updates = self._query_to_df(query_rxn_field_update)
        if not df_rxn_updates.empty:
            updated_entries_list.append(df_rxn_updates)

        if prev_rxn_table_created_and_used:
            try:
                self.cursor.execute(f"DROP TABLE IF EXISTS {temp_prev_rxn_table_name};")
                self.logger.info(f"Dropped temporary reaction snapshot table: {temp_prev_rxn_table_name}")
            except psycopg.Error as e:
                self.logger.warning(f"Could not drop {temp_prev_rxn_table_name}: {e}")

        # --- Final Concatenation ---
        if not updated_entries_list:
            # Return empty DataFrame with all expected columns
            return pd.DataFrame(columns=['entity_key', 'entity_type', 'field_name',
                                         'old_value_in_production', 'new_value_from_current_sl',
                                         'previous_value_from_snapshot_sl',
                                         'is_manual_override_in_production',
                                         'is_change_in_sl_since_snapshot'])

        final_updates_df = pd.concat(updated_entries_list, ignore_index=True)
        self.logger.info(f"Total {len(final_updates_df)} potential field updates found (3-way).")
        return final_updates_df

    def _generate_single_delta_file(self, df: pd.DataFrame, filename_suffix: str,
                                    output_directory: Path):  # Changed param
        """Saves a DataFrame to a CSV file in the specified output directory."""
        if df is None or df.empty:
            self.logger.info(f"No data to write for '{filename_suffix}'. Skipping file generation.")
            return

        output_directory.mkdir(parents=True, exist_ok=True)  # Ensure dir exists
        output_path = output_directory / filename_suffix  # Use passed directory
        try:
            df.to_csv(output_path, index=False, sep='\t')
            self.logger.info(f"Generated delta file: {output_path}")
        except Exception as e:
            self.logger.error(f"Failed to write delta file '{output_path}': {e}")

    def _generate_summary_report(self, summary_data: Dict[str, Any], filename_suffix: str, output_directory: Path): # Changed param
        """Generates a text summary report in the specified output directory."""
        output_directory.mkdir(parents=True, exist_ok=True)  # Ensure dir exists
        output_path = output_directory / filename_suffix  # Use passed directory

        report_content = f"LipidCRED Database Update Summary\n"
        report_content += f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        report_content += f"Data Timestamp (source files): {summary_data.get('data_timestamp', 'N/A')}\n\n"  # You might pass this in

        report_content += f"New Organisms Found: {summary_data.get('new_organisms', 0)}\n"
        report_content += f"New Molecules Found: {summary_data.get('new_molecules', 0)}\n"
        report_content += f"New Reaction Definitions Found: {summary_data.get('new_reaction_definitions', 0)}\n"
        report_content += f"New Enzyme Definitions Found: {summary_data.get('new_enzyme_definitions', 0)}\n"
        report_content += f"New Reaction-Enzyme links Found: {summary_data.get('new_reaction_links', 0)}\n"
        report_content += f"Potential Field Updates Found: {summary_data.get('updated_entries', 0)}\n\n"

        report_content += "Review the generated CSV files in the same directory.\n"
        report_content += "After review, rename relevant files to 'approved_...' and use the 'apply' mode.\n"

        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(report_content)
            self.logger.info(f"Generated summary report: {output_path}")
        except Exception as e:
            self.logger.error(f"Failed to write summary report '{output_path}': {e}")

    def _generate_new_additions_report(self, run_output_dir: Path):
        report_path = run_output_dir / "report_newly_inserted_entities.txt"
        self.logger.info(f"Generating report of newly inserted entities to: {report_path}")

        try:
            with open(report_path, "w", encoding='utf-8') as f:
                f.write("Summary of New Entities Inserted or Linked During 'apply' Phase\n")
                f.write("=" * 70 + "\n\n")

                if not self.new_additions_summary:
                    f.write("No new entities were inserted into the database during this run.\n")
                    return

                for entity_type, items in self.new_additions_summary.items():
                    f.write(f"--- Newly Inserted {entity_type.capitalize()} ({len(items)}) ---\n")
                    if items:
                        for item in sorted(items):  # Sort for consistent output
                            f.write(f"- {item}\n")
                    else:
                        f.write("None\n")
                    f.write("\n")
        except Exception as e:
            self.logger.error(f"Failed to write new additions report: {e}")

    def _generate_abbreviation_warning_report(self, run_output_dir: Path):
        """Generates a report of products in newly created pairs that are missing abbreviations."""
        if not self.abbreviation_warnings:
            self.logger.info("No products with missing abbreviations were found in new reaction pairs.")
            return

        report_path = run_output_dir / "report_reaction_products_missing_abbreviations.tsv"
        self.logger.warning(
            f"Found {len(self.abbreviation_warnings)} products in new reaction pairs missing critical abbreviations. "
            f"Generating report at: {report_path}")

        # Use pandas to easily handle list of dicts and avoid duplicates
        df = pd.DataFrame(self.abbreviation_warnings)
        df.drop_duplicates(inplace=True)

        df.to_csv(report_path, sep='\t', index=False)


    # Placeholder for the main orchestrator method for the 'generate' mode
    def process_files_for_review(self, lipids_file_path: str, enzymes_file_path: str,
                                 lipids2uniprot_file_path: Optional[str] = None,
                                 current_data_run_timestamp: Optional[str] = None):  # Renamed data_timestamp_str
        """
        Top-level method to ingest data into staging and then generate delta files for review.
        This is the main entry point for the 'generate' mode.
        """
        run_ts = current_data_run_timestamp if current_data_run_timestamp else datetime.now().strftime(
            "%Y%m%d_%H%M%S_auto")
        self.logger.info(f"Processing source files for review. Run timestamp: {run_ts}")

        # Step 1: Ingest data into staging tables
        # This uses the already updated ingest_swisslipids_data
        self.ingest_swisslipids_data(lipids_file_path, enzymes_file_path, lipids2uniprot_file_path)

        # Step 2: Generate delta files based on comparison.
        # For the first run or if no snapshots, pass None for snapshot paths.
        snapshot_dir = self.output_dir / "snapshots"  # Define where snapshots are stored/expected
        prev_lipids_snapshot_path = str(snapshot_dir / "previous_lipids_snapshot.tsv")
        prev_enzymes_snapshot_path = str(snapshot_dir / "previous_enzymes_snapshot.tsv")

        # Check if snapshot files actually exist to pass to generate_delta_files_for_review
        actual_prev_lipids_path = prev_lipids_snapshot_path if Path(prev_lipids_snapshot_path).exists() else None
        actual_prev_enzymes_path = prev_enzymes_snapshot_path if Path(prev_enzymes_snapshot_path).exists() else None

        summary = self.generate_delta_files_for_review(
            current_data_run_timestamp=run_ts,  # This is for naming the output sub-directory for this run
            previous_snapshot_lipids_file_path=actual_prev_lipids_path,
            previous_snapshot_enzymes_file_path=actual_prev_enzymes_path
        )

        self.logger.info(f"Review file generation process completed. Summary: {summary}")
        return summary

    def _insert_entity(self, table_name: str, conflict_columns: List[str], data_dict: Dict[str, Any],
                       return_id_col: Optional[str] = None) -> Tuple[Optional[int], bool]:
        """
        Helper to insert a single entity. If it already exists based on conflict_columns,
        it does nothing (ON CONFLICT DO NOTHING). Returns the ID and a boolean indicating if a new row was created.
        Assumes `self.cursor` is active.

        This method handles both dict and tuple results to support different cursor configurations.
        """
        if not data_dict:
            self.logger.warning(f"Data dictionary is empty for inserting into {table_name}. Skipping.")
            return None, False
        if not conflict_columns and return_id_col:
            raise ValueError("Conflict columns list cannot be empty for _insert_entity with ON CONFLICT.")

        # Prepare for INSERT
        cols_to_insert_quoted = [f'"{k}"' for k in data_dict.keys()]
        vals_placeholders = ["%s"] * len(data_dict)
        actual_values = list(data_dict.values())

        conflict_target_sql = ", ".join([f'"{c}"' for c in conflict_columns])
        returning_sql = f"RETURNING \"{return_id_col}\"" if return_id_col else ""

        sql = f"""
            INSERT INTO {table_name} ({', '.join(cols_to_insert_quoted)})
            VALUES ({', '.join(vals_placeholders)})
            ON CONFLICT ({conflict_target_sql}) DO NOTHING
            {returning_sql};
        """

        try:
            self.cursor.execute(sql, tuple(actual_values))
            inserted_id: Optional[int] = None
            was_inserted: bool = False

            if return_id_col:
                result = self.cursor.fetchone()
                if result:
                    # New row was inserted
                    inserted_id = self._extract_value_from_result(result, return_id_col, 0)
                    was_inserted = True
                else:
                    # Conflict occurred, fetch existing row
                    self.logger.debug(f"Insert into {table_name} resulted in DO NOTHING. "
                                      f"Fetching existing ID using conflict columns: {conflict_columns}.")

                    # Build WHERE clause for SELECT
                    select_where_clauses = []
                    select_where_values = []
                    for col_name in conflict_columns:
                        if col_name not in data_dict:
                            raise ValueError(f"Conflict column '{col_name}' missing from data_dict for {table_name}")

                        value = data_dict[col_name]
                        if value is None:
                            select_where_clauses.append(f'"{col_name}" IS NULL')
                        else:
                            select_where_clauses.append(f'"{col_name}" = %s')
                            select_where_values.append(value)

                    if select_where_clauses:
                        select_sql = f'SELECT "{return_id_col}" FROM {table_name} WHERE {" AND ".join(select_where_clauses)} LIMIT 1;'
                        self.cursor.execute(select_sql, tuple(select_where_values))
                        existing_result = self.cursor.fetchone()

                        if existing_result:
                            inserted_id = self._extract_value_from_result(existing_result, return_id_col, 0)
                            self.logger.debug(f"Fetched existing ID {inserted_id} for conflict keys in {table_name}.")
                        else:
                            conflict_values_for_log = {c: data_dict.get(c, "MISSING") for c in conflict_columns}
                            self.logger.warning(
                                f"Could not retrieve ID for conflict values {conflict_values_for_log} in {table_name}")

            return inserted_id, was_inserted

        except psycopg.Error as e:
            self.logger.error(
                f"Database error in _insert_entity for table {table_name}: "
                f"{e.diag.message_primary if e.diag else e}",
                exc_info=True)
            raise
        except Exception as e:
            self.logger.error(
                f"Unexpected error in _insert_entity for table {table_name}: {e}",
                exc_info=True)
            raise

    def _extract_value_from_result(self, result, column_name: str, index: int = 0):
        """
        Extract a value from a query result that could be either a dict or tuple.

        Args:
            result: The result from fetchone()
            column_name: The column name to extract (for dict results)
            index: The column index to extract (for tuple results)

        Returns:
            The extracted value
        """
        if result is None:
            return None

        if isinstance(result, dict):
            return result.get(column_name)
        else:
            # Tuple or list result
            return result[index] if index < len(result) else None

    def _apply_approved_new_organisms(self, approved_new_organisms_file: Path):
        if not approved_new_organisms_file.exists():
            self.logger.info(f"No approved new organisms file found: {approved_new_organisms_file}")
            return 0

        df = pd.read_csv(approved_new_organisms_file, sep='\t', dtype=str, keep_default_na=False, na_values=[''])
        count = 0
        for _, row in df.iterrows():
            # organism_id from CSV is text, production table is INTEGER
            org_id_val = row.get('organism_id')
            if not org_id_val or not org_id_val.strip():
                self.logger.warning(f"Skipping organism with missing ID: {row}")
                continue
            try:
                organism_data = {
                    "organism_id": int(org_id_val),
                    "taxon_scientific_name": row.get('taxon_scientific_name')
                }
                self._insert_entity("lipograph.organism", ["organism_id"], organism_data)
                count += 1
            except Exception as e:
                self.logger.error(f"Failed to insert new organism {row.to_dict()}: {e}")
        self.logger.info(f"Applied {count}/{len(df)} new organisms.")
        return count

    def _apply_approved_new_molecules(self, approved_new_molecules_file: Path):
        if not approved_new_molecules_file.exists():
            self.logger.info(f"No approved new molecules file found: {approved_new_molecules_file}")
            return 0

        df = pd.read_csv(approved_new_molecules_file, sep='\t', dtype=str, keep_default_na=False, na_values=[''])
        count = 0
        for _, row in df.iterrows():
            # Clean abbreviation before insert
            cleaned_abbr = row.get('abbreviation', '')  # Use original abbreviation for cleaning
            if pd.isna(cleaned_abbr) or not cleaned_abbr:  # if abbreviation is empty, cleaned is also empty/None
                final_cleaned_abbr = None
            else:
                final_cleaned_abbr = self.lipid_translator.clean_lipid_abbreviation(cleaned_abbr)
                if pd.isna(final_cleaned_abbr): final_cleaned_abbr = None

            molecule_data = {
                "molecule_name": row.get('molecule_name'),
                "swisslipids_id": row.get('swisslipids_id'),  # Conflict target
                "abbreviation": row.get('abbreviation') if pd.notna(row.get('abbreviation')) else None,
                "cleaned_abbreviation": final_cleaned_abbr,
                "molecule_synonyms": row.get('molecule_synonyms') if pd.notna(row.get('molecule_synonyms')) else None,
                "sl_lipid_class": row.get('sl_lipid_class') if pd.notna(row.get('sl_lipid_class')) else None
            }
            if not molecule_data["swisslipids_id"]:
                self.logger.warning(f"Skipping molecule with missing swisslipids_id: {row.to_dict()}")
                continue
            try:
                self._insert_entity("lipograph.molecules", ["swisslipids_id"], molecule_data,
                                    return_id_col="molecule_id")
                count += 1
            except Exception as e:
                self.logger.error(f"Failed to insert new molecule {row.to_dict()}: {e}")
        self.logger.info(f"Applied {count}/{len(df)} new molecules.")
        return count

    def _apply_approved_new_reaction_definitions(self, approved_new_reactions_file: Path):
        if not approved_new_reactions_file.exists():
            self.logger.info(f"No approved new reactions file found: {approved_new_reactions_file}")
            return 0

        df = pd.read_csv(approved_new_reactions_file, sep='\t', dtype=str, keep_default_na=False, na_values=[''])
        count = 0
        for _, row in df.iterrows():
            rhea_id_val = row.get('rhea_id')
            reaction_data = {
                "reaction_text": row.get('reaction_text'),  # Part of conflict target
                "rhea_id": int(float(rhea_id_val)) if pd.notna(rhea_id_val) and rhea_id_val.strip() != '' else None,
                # Part of conflict target
                "doi": row.get('doi') if pd.notna(row.get('doi')) else None
            }
            if not reaction_data["reaction_text"]:
                self.logger.warning(f"Skipping reaction with missing reaction_text: {row.to_dict()}")
                continue
            try:
                self._insert_entity("lipograph.reactions", ["reaction_text", "rhea_id"], reaction_data,
                                    return_id_col="reaction_id")
                count += 1
            except Exception as e:
                self.logger.error(f"Failed to insert new reaction definition {row.to_dict()}: {e}")
        self.logger.info(f"Applied {count}/{len(df)} new reaction definitions.")
        return count

    def _apply_approved_new_enzyme_definitions(self, approved_new_enzymes_file: Path):
        if not approved_new_enzymes_file.exists():
            self.logger.info(f"No approved new enzymes file found: {approved_new_enzymes_file}")
            return 0

        df = pd.read_csv(approved_new_enzymes_file, sep='\t', dtype=str, keep_default_na=False, na_values=[''])
        count = 0
        for _, row in df.iterrows():
            org_id_val = row.get('organism_id')
            enzyme_data = {
                "enzyme_name": row.get('enzyme_name') if pd.notna(row.get('enzyme_name')) else None,
                "uniprot_id": row.get('uniprot_id'),  # Conflict target
                "organism_id": int(float(org_id_val)) if pd.notna(org_id_val) and org_id_val.strip() != '' else None,
                "swisslipids_p_id": row.get('swisslipids_p_id') if pd.notna(row.get('swisslipids_p_id')) else None
            }
            if not enzyme_data["uniprot_id"]:
                self.logger.warning(f"Skipping enzyme with missing uniprot_id: {row.to_dict()}")
                continue
            try:
                self._insert_entity("lipograph.enzymes", ["uniprot_id"], enzyme_data, return_id_col="enzyme_id")
                count += 1
            except Exception as e:
                self.logger.error(f"Failed to insert new enzyme definition {row.to_dict()}: {e}")
        self.logger.info(f"Applied {count}/{len(df)} new enzyme definitions.")
        return count

    def _apply_approved_reaction_links(self, approved_links_file: Path):
        """
        Reads the approved links file and processes only those rows.
        """
        if not approved_links_file.exists():
            self.logger.info("No approved new reaction links file found. Skipping.")
            return

        links_to_process_df = pd.read_csv(approved_links_file, sep='\t', dtype=str)

        if links_to_process_df.empty:
            self.logger.info("Approved new links file is empty.")
            return

        self.logger.info(f"[Linker] Starting to process {len(links_to_process_df)} approved new reaction links...")
        for idx, row in links_to_process_df.iterrows():
            self._process_and_link_reaction_entry(row.to_dict())

    def _resolve_molecule_by_abbr(self, cleaned_abbreviation: str) -> Optional[Dict[str, Any]]:
        """
        Resolves a molecule's full details from its cleaned_abbreviation.
        1. Checks the curated family_molecule_id_map first.
        2. If not found, queries the database by cleaned_abbreviation.
        3. Handles ambiguity by sorting and selecting the first result.
        """
        if not cleaned_abbreviation or not cleaned_abbreviation.strip():
            self.logger.warning("Attempted to resolve an empty or None cleaned_abbreviation.")
            return None

        self.logger.debug(f"Resolving molecule by cleaned_abbreviation: '{cleaned_abbreviation}'")

        # Step 1: Check curated family map
        if cleaned_abbreviation in self.family_molecule_id_map:
            molecule_id = self.family_molecule_id_map[cleaned_abbreviation]
            self.logger.info(f"Resolved '{cleaned_abbreviation}' to molecule_id {molecule_id} via curated family map.")
            # ## CORRECTED: Use the new helper method ##
            return self._get_molecule_details(molecule_id)

        # Step 2: Query the database by cleaned_abbreviation
        query_db = "SELECT * FROM lipograph.molecules WHERE cleaned_abbreviation = %s;"
        df_result = self._query_to_df(query_db, (cleaned_abbreviation,))

        if df_result.empty:
            self.logger.debug(f"Molecule with cleaned_abbreviation '{cleaned_abbreviation}' not found in DB.")
            return None

        if len(df_result) > 1:
            df_result_sorted = df_result.sort_values(by='molecule_id')
            chosen_molecule_data = df_result_sorted.iloc[0].to_dict()
            self.logger.warning(
                f"Ambiguity: Found {len(df_result)} molecules with cleaned_abbreviation '{cleaned_abbreviation}'. "
                f"Sorted by ID, chose ID: {chosen_molecule_data['molecule_id']}. "
                f"All found IDs: {df_result['molecule_id'].tolist()}"
            )
            return chosen_molecule_data

        # Exactly one result
        resolved_data = df_result.iloc[0].to_dict()
        self.logger.debug(
            f"Resolved '{cleaned_abbreviation}' to molecule_id {resolved_data['molecule_id']} via database query.")
        return resolved_data

    @lru_cache(maxsize=8192)
    def _resolve_molecule_name_from_sl_text(self, sl_molecule_name: str) -> Optional[Dict[str, Any]]:
        """
        Resolves a molecule name from SwissLipids reaction text to its full production record.
        Returns a dictionary of the molecule's data or None if not found.
        """
        if not sl_molecule_name or not sl_molecule_name.strip():
            self.logger.warning("Attempted to resolve an empty SL molecule name.")
            return None


        original_sl_name_stripped = sl_molecule_name.strip()
        self.logger.debug(f"Resolving SL reaction text component: '{original_sl_name_stripped}'")

        # --- Attempt 1: Exact match on molecule_name ---
        query_exact_name = "SELECT * FROM lipograph.molecules WHERE molecule_name = %s;"
        df_exact_name = self._query_to_df(query_exact_name, (original_sl_name_stripped,))
        if not df_exact_name.empty:
            # Handle potential ambiguity by picking the first result after sorting by ID
            chosen_row = df_exact_name.sort_values(by='molecule_id').iloc[0]
            self.logger.debug(
                f"Resolved '{original_sl_name_stripped}' by exact name match to ID {chosen_row['molecule_id']}.")
            return chosen_row.to_dict()

        # --- Attempt 2: Match within molecule_synonyms ---
        synonym_pattern = f"%|{original_sl_name_stripped}|%"
        query_synonym = "SELECT * FROM lipograph.molecules WHERE molecule_synonyms LIKE %s;"
        df_synonym = self._query_to_df(query_synonym, (synonym_pattern,))
        if not df_synonym.empty:
            chosen_row = df_synonym.sort_values(by='molecule_id').iloc[0]
            self.logger.debug(
                f"Resolved '{original_sl_name_stripped}' by synonym match to ID {chosen_row['molecule_id']}.")
            return chosen_row.to_dict()

        # --- Attempt 3: Check curated family_molecule_id_map ---
        if original_sl_name_stripped in self.family_molecule_id_map:
            molecule_id = self.family_molecule_id_map[original_sl_name_stripped]
            self.logger.info(
                f"Resolved '{original_sl_name_stripped}' to molecule_id {molecule_id} via curated family map.")
            # ## CORRECTED: Use the internal ID helper ##
            return self._get_molecule_details(molecule_id)

        # --- Attempt 4: Clean the SL name and resolve by abbreviation ---
        cleaned_abbr_from_sl_name = ""
        try:
            if self.lipid_translator:
                cleaned_abbr_from_sl_name = self.lipid_translator.clean_lipid_abbreviation(
                    original_sl_name_stripped)
        except Exception as e_trans:
            self.logger.error(f"Failed to clean SL name '{original_sl_name_stripped}': {e_trans}")

        if cleaned_abbr_from_sl_name:
            # This call now returns a dict or None
            mol_info_dict = self._resolve_molecule_by_abbr(cleaned_abbr_from_sl_name)
            if mol_info_dict:
                self.logger.debug(
                    f"Resolved '{original_sl_name_stripped}' (cleaned to '{cleaned_abbr_from_sl_name}') to ID {mol_info_dict['molecule_id']}.")
                return mol_info_dict

        # If all local attempts fail, return None. The calling function will handle logging this as unresolved.
        return None

        # # --- If all attempts fail ---
        # self.logger.warning(f"Could not resolve SL reaction component '{original_sl_name_stripped}' to any molecule_id "
        #                     f"(cleaned attempt: '{cleaned_abbr_from_sl_name}').")
        # return None

    def _normalize_stereochemistry(self, sl_molecule_name: str) -> Optional[str]:
        """
        Attempts to normalize a molecule name by handling stereochemistry notation
        like "(11R,12R)-..." -> "11,12-...". It now includes a check to ensure
        'R' or 'S' is present to avoid false positives.
        """
        # Step 1: Match any name starting with "(...)-"
        match = re.match(r'\(([^)]+)\)-(.*)', sl_molecule_name)

        if match:
            stereo_part = match.group(1)  # e.g., "11R,12R" or "9Z,12Z"
            rest_of_name = match.group(2)

            # Step 2: **Safety Check** - Only proceed if it's actually stereochemistry
            if 'R' in stereo_part or 'S' in stereo_part:
                # It's a match, now perform the transformation
                numbers_only = re.sub(r'[RS]', '', stereo_part)
                normalized_name = f"{numbers_only}-{rest_of_name}"

                self.logger.info(f"Normalized stereochemistry: '{sl_molecule_name}' -> '{normalized_name}'")
                return normalized_name
            # If 'R' or 'S' are not in the parentheses, it's not the pattern we want to fix.
            # Example: "(9Z,12Z)-..." will be ignored.

        return None

    def _strip_html(self, text: str) -> str:
        if not text: return ""
        from bs4 import BeautifulSoup  # Add 'pip install beautifulsoup4' to requirements
        soup = BeautifulSoup(text, "html.parser")
        return soup.get_text()

    @lru_cache(maxsize=4036)  # Cache up to 1024 unique API lookups
    def _resolve_lipid_via_swisslipids_api(self,
                                           unresolved_lipid_name: str,
                                           reaction_context_slp_id: str,
                                           reaction_context_rhea_id_str: Optional[str]) -> Optional[
        str]:  # Returns SLM ID
        import requests
        from thefuzz import fuzz

        if not hasattr(self, 'api_slp_cache'):
            self.api_slp_cache = {}

        api_base_url = "https://www.swisslipids.org/api/index.php/entity/"

        try:
            if reaction_context_slp_id in self.api_slp_cache:
                slp_data = self.api_slp_cache[reaction_context_slp_id]
            else:
                self.logger.debug(f"API Call: Fetching {api_base_url + reaction_context_slp_id}")
                response = requests.get(api_base_url + reaction_context_slp_id, timeout=10)
                response.raise_for_status()
                slp_data = response.json()
                self.api_slp_cache[reaction_context_slp_id] = slp_data
        except requests.exceptions.RequestException as e:
            self.logger.warning(f"API request failed for SLP ID {reaction_context_slp_id}: {e}")
            return None
        except json.JSONDecodeError as e:
            self.logger.warning(f"Failed to decode JSON from API for SLP ID {reaction_context_slp_id}: {e}")
            return None

        if not slp_data or 'reactions' not in slp_data or not slp_data['reactions']:
            self.logger.debug(f"No reactions found in API response for SLP ID {reaction_context_slp_id}")
            return None

        # --- Simplified Matching Logic ---
        MIN_ACCEPTABLE_SCORE = 90  # Tune this! If names are very close, this can be high.
        # For "acyl-CoA" vs "acyl-CoAs", fuzz.ratio might be ~95
        # fuzz.token_set_ratio might be 100. Let's try fuzz.ratio for stricter length.

        target_rhea_id_int = None
        if reaction_context_rhea_id_str and reaction_context_rhea_id_str.strip().isdigit():
            target_rhea_id_int = int(reaction_context_rhea_id_str.strip())

        candidate_molecules_in_reaction = []  # List of (slm_id, cleaned_api_name)

        found_target_reaction = False
        if target_rhea_id_int is not None:
            for slcr_id, reaction_data_wrapper in slp_data['reactions'].items():
                reaction_details = reaction_data_wrapper.get('reaction', {})
                rhea_info = reaction_details.get('rhea', {})
                api_rhea_id_str = rhea_info.get('rhea_id')

                if api_rhea_id_str and api_rhea_id_str.isdigit() and int(api_rhea_id_str) == target_rhea_id_int:
                    found_target_reaction = True
                    for participant_type in ['rhea_reactants', 'rhea_products']:
                        participants = rhea_info.get(participant_type, {})
                        for slm_id, molecule_info in participants.items():
                            api_molecule_name_raw = molecule_info.get('name', '')
                            api_molecule_name_clean = self._strip_html(api_molecule_name_raw).strip()
                            if api_molecule_name_clean:
                                candidate_molecules_in_reaction.append((slm_id, api_molecule_name_clean))
                    break
            if not found_target_reaction:
                self.logger.debug(
                    f"Target Rhea ID {target_rhea_id_int} not found in API response for SLP {reaction_context_slp_id}.")
                return None
        else:
            self.logger.debug(
                f"No Rhea ID context for SLP {reaction_context_slp_id}, API fallback requires Rhea ID for precision. Skipping API lookup for '{unresolved_lipid_name}'.")
            return None  # Require Rhea ID to ensure we're looking at the right reaction context.

        if not candidate_molecules_in_reaction:
            self.logger.debug(
                f"No candidate molecules extracted from API for SLP {reaction_context_slp_id}, Rhea {target_rhea_id_int} to match '{unresolved_lipid_name}'.")
            return None

        # Score candidates from the specific reaction
        best_match_slm_id = None
        highest_score = -1
        ambiguous_matches = []  # Store (score, slm_id, name) for logging if ambiguous

        for slm_id, api_name in candidate_molecules_in_reaction:
            # Using fuzz.ratio is stricter about sequence and length than token_set_ratio
            # This should help differentiate "L-serine" from "xxx-L-serine" better.
            score = fuzz.ratio(unresolved_lipid_name.lower(), api_name.lower())

            self.logger.debug(
                f"API Fuzzy: '{unresolved_lipid_name}' vs API_NAME '{api_name}' (SLM: {slm_id}) -> score {score}")

            if score >= MIN_ACCEPTABLE_SCORE:
                if score > highest_score:
                    highest_score = score
                    best_match_slm_id = slm_id
                    ambiguous_matches = [(score, slm_id, api_name)]  # New best, reset ambiguity list
                elif score == highest_score:
                    # Same score as current best, this is an ambiguity
                    ambiguous_matches.append((score, slm_id, api_name))

        if not best_match_slm_id:
            self.logger.info(
                f"API Fallback: No matches >= {MIN_ACCEPTABLE_SCORE} for '{unresolved_lipid_name}' in SLP {reaction_context_slp_id}, Rhea {target_rhea_id_int}.")
            return None

        if len(ambiguous_matches) > 1:
            # This means multiple different SLM IDs (or same SLM appearing multiple times with same name) achieved the same highest score
            # Or multiple different names achieved the same high score for different SLMs.
            # We need to check if all SLM IDs in ambiguous_matches are actually the same.
            unique_slm_ids_in_ambiguity = set(item[1] for item in ambiguous_matches)
            if len(unique_slm_ids_in_ambiguity) > 1:
                self.logger.warning(
                    f"API Fallback AMBIGUITY for '{unresolved_lipid_name}' (SLP {reaction_context_slp_id}, Rhea {target_rhea_id_int}): "
                    f"Multiple distinct SLM IDs achieved score {highest_score}. Matches: {ambiguous_matches}. Skipping."
                )
                return None
            else:
                # All ambiguous matches point to the SAME SLM_ID, just different textual representations perhaps, or same name.
                # This is acceptable. The best_match_slm_id is already set to one of these.
                self.logger.info(
                    f"API Fallback: Multiple representations for SLM ID '{best_match_slm_id}' achieved score {highest_score} for '{unresolved_lipid_name}'. "
                    f"Selected SLM: {best_match_slm_id} (e.g., from API name '{ambiguous_matches[0][2]}')."
                )
                return best_match_slm_id

        # Only one clear best match (or multiple identical matches for the same SLM ID)
        self.logger.info(
            f"API Fallback: Best match for '{unresolved_lipid_name}' is API_NAME '{ambiguous_matches[0][2]}' (SLM: {best_match_slm_id}) "
            f"with score {highest_score} (SLP {reaction_context_slp_id}, Rhea {target_rhea_id_int})."
        )
        return best_match_slm_id

    def _ensure_organism(self, organism_id_str: str, organism_name: Optional[str] = None) -> Optional[int]:
        """Ensure organism exists, insert if not, return its integer ID."""
        if not organism_id_str or not organism_id_str.strip():
            self.logger.error("Organism ID is missing for ensure_organism.")
            return None
        try:
            organism_id = int(organism_id_str.strip())
        except ValueError:
            self.logger.error(f"Invalid organism ID format: {organism_id_str}")
            return None

        # ## MODIFIED: Capture `was_inserted` and log/summarize ##
        prod_organism_id, was_inserted = self._insert_entity(
            "lipograph.organism", ["organism_id"],
            {"organism_id": organism_id, "taxon_scientific_name": organism_name},
            return_id_col="organism_id"
        )

        if was_inserted:
            self.logger.info(f"[Linker] INSERTED new organism: ID {organism_id}, Name: {organism_name}")
            self.new_additions_summary['organisms'].append(f"ID={organism_id}, Name={organism_name}")
        else:
            self.logger.debug(f"[Linker] Found existing organism: ID {organism_id}")

        return prod_organism_id

    def _ensure_enzyme(self, uniprot_id: str, gene_name: Optional[str], organism_id: int,
                       swisslipids_p_id: Optional[str] = None) -> Optional[int]:
        """Ensure enzyme exists by UniProt ID + organism ID, insert if not, return its enzyme_id."""
        if not uniprot_id or not uniprot_id.strip():
            self.logger.error("[Linker] UniProt ID is missing for ensure_enzyme.")
            return None
        if organism_id is None:
            self.logger.error(f"[Linker] Organism ID is missing for enzyme {uniprot_id}.")
            return None

        enzyme_data = {
            "uniprot_id": uniprot_id.strip(),
            "enzyme_name": gene_name if pd.notna(gene_name) else None,
            "organism_id": organism_id,
            "swisslipids_p_id": swisslipids_p_id if pd.notna(swisslipids_p_id) else None
        }

        # ## MODIFIED: Capture `was_inserted` and log/summarize ##
        prod_enzyme_id, was_inserted = self._insert_entity(
            "lipograph.enzymes", ["uniprot_id"],  # Corrected conflict key
            enzyme_data,
            return_id_col="enzyme_id"
        )

        if was_inserted:
            log_msg = f"[Linker] INSERTED new enzyme: ID {prod_enzyme_id}, UniProt: {uniprot_id}, Organism: {organism_id}, Name: {gene_name or 'N/A'}"
            self.logger.info(log_msg)
            self.new_additions_summary['enzymes'].append(log_msg)
        else:
            self.logger.debug(
                f"[Linker] Found existing enzyme: ID {prod_enzyme_id} for UniProt {uniprot_id}, Organism {organism_id}")
            # Optional: Add logic here to UPDATE swisslipids_p_id if needed

        return prod_enzyme_id

    def _ensure_reaction(self, reaction_text: str, rhea_id_str: Optional[str], doi: Optional[str]) -> Optional[int]:
        """Ensure reaction definition exists, insert if not, return its reaction_id."""
        if not reaction_text or not reaction_text.strip():
            self.logger.error("[Linker] Reaction text is missing for ensure_reaction.")
            return None

        rhea_id = None
        if rhea_id_str and pd.notna(rhea_id_str) and str(rhea_id_str).strip():
            try:
                rhea_id = int(float(str(rhea_id_str).strip()))
            except ValueError:
                self.logger.warning(f"[Linker] Invalid Rhea ID format '{rhea_id_str}', storing as NULL.")

        reaction_data = {
            "reaction_text": reaction_text.strip(),
            "rhea_id": rhea_id,
            "doi": doi if doi and pd.notna(doi) and doi.strip() else None
        }

        # ## MODIFIED: Capture `was_inserted` and log/summarize ##
        prod_reaction_id, was_inserted = self._insert_entity(
            "lipograph.reactions", ["reaction_text", "rhea_id"],
            reaction_data,
            return_id_col="reaction_id"
        )

        if was_inserted:
            log_msg = f"[Linker] INSERTED new reaction: ID {prod_reaction_id}, Text: '{reaction_text[:80]}...', Rhea: {rhea_id or 'N/A'}"
            self.logger.info(log_msg)
            self.new_additions_summary['reactions'].append(log_msg)
        else:
            self.logger.debug(f"[Linker] Found existing reaction: ID {prod_reaction_id}")

        return prod_reaction_id

    def _ensure_reaction_enzyme(self, reaction_id: int, enzyme_id: int) -> Optional[int]:
        """Ensure reaction_enzyme link exists, insert if not, return its reaction_enzyme_id."""
        if reaction_id is None or enzyme_id is None:
            self.logger.error(
                f"[Linker] Missing reaction_id ({reaction_id}) or enzyme_id ({enzyme_id}) for reaction_enzyme link.")
            return None

        # ## MODIFIED: Capture `was_inserted` and log/summarize ##
        prod_reaction_enzyme_id, was_inserted = self._insert_entity(
            "lipograph.reaction_enzyme", ["reaction_id", "enzyme_id"],
            {"reaction_id": reaction_id, "enzyme_id": enzyme_id},
            return_id_col="reaction_enzyme_id"
        )

        if was_inserted:
            log_msg = f"[Linker] INSERTED new reaction_enzyme link: RE_ID {prod_reaction_enzyme_id} (Reaction ID: {reaction_id}, Enzyme ID: {enzyme_id})"
            self.logger.info(log_msg)
            self.new_additions_summary['reaction_enzyme_links'].append(log_msg)
        else:
            self.logger.debug(f"[Linker] Found existing reaction_enzyme link: RE_ID {prod_reaction_enzyme_id}")

        return prod_reaction_enzyme_id

    def _ensure_reaction_pairs(self,
                               reaction_enzyme_id: int,
                               reactant_molecule_ids: List[int],
                               product_molecule_ids: List[int],
                               is_sl_reversible: bool) -> int:
        """
        Ensures reaction_pairs exist for all generated combinations of reactant and product IDs.
        Returns the count of pairs successfully created or found.
        """
        if reaction_enzyme_id is None:
            self.logger.error("[Linker] Missing reaction_enzyme_id for creating reaction_pairs.")
            return 0

        if not reactant_molecule_ids or not product_molecule_ids:
            self.logger.warning(
                f"[Linker] Empty reactant or product ID list for RE_ID:{reaction_enzyme_id}. No pairs to create.")
            return 0

        id_combinations = generate_reaction_combinations(
            reactant_molecule_ids,
            product_molecule_ids,
            is_sl_reversible
        )

        created_or_found_count = 0
        # No longer need ignored_ids here, as filtering happens before this method is called.

        for r_id, p_id in id_combinations:
            # ## MODIFIED: Capture `was_inserted` and log/summarize for pairs ##
            pair_id, was_inserted = self._insert_entity(
                table_name="lipograph.reaction_pairs",
                conflict_columns=["reaction_enzyme_id", "reactant_molecule_id", "product_molecule_id"],
                data_dict={
                    "reaction_enzyme_id": reaction_enzyme_id,
                    "reactant_molecule_id": r_id,
                    "product_molecule_id": p_id
                },
                return_id_col="pair_id"
            )

            if pair_id is not None:
                created_or_found_count += 1
                if was_inserted:
                    log_msg = f"[Linker] INSERTED new reaction pair: Pair ID {pair_id} (RE_ID: {reaction_enzyme_id}, R_mol: {r_id}, P_mol: {p_id})"
                    self.logger.info(log_msg)
                    # Adding to summary might be too verbose, but possible:
                    self.new_additions_summary['reaction_pairs'].append(log_msg)
                else:
                    self.logger.debug(
                        f"[Linker] Found existing reaction pair: Pair ID {pair_id} (RE_ID: {reaction_enzyme_id}, R_mol: {r_id}, P_mol: {p_id})")

        # ## REFINEMENT: Log summary for this specific reaction-enzyme link ##
        self.logger.info(
            f"[Linker] Ensured {created_or_found_count} of {len(id_combinations)} possible reaction pairs for RE_ID:{reaction_enzyme_id}.")
        return created_or_found_count

    def _apply_manual_reactions(self, manual_reactions_file: Path) -> int:
        """
        Processes the approved manual reactions file by re-using the main linking logic.
        """
        if not manual_reactions_file.exists():
            self.logger.info(f"No manual reactions file found: {manual_reactions_file}. Skipping.")
            return 0

        self.logger.info(f"Processing manual reactions from {manual_reactions_file}...")
        manual_df = pd.read_csv(manual_reactions_file, sep='\t', dtype=str).fillna('')
        if manual_df.empty: return 0

        successful_additions = 0
        for idx, row in manual_df.iterrows():
            try:
                reaction_text_abbr = row.get('InputReactionText_CleanedAbbr', '').strip()
                if not reaction_text_abbr: continue

                # 1. Parse the abbreviation-based text
                parsed_components, is_reversible = split_reaction_text(reaction_text_abbr)
                if not parsed_components:
                    self.logger.warning(f"Could not parse manual reaction row {idx + 1}: {reaction_text_abbr}")
                    continue

                # 2. Resolve abbreviations to get full molecule names
                reactant_details = [self._resolve_molecule_by_abbr(abbr) for abbr in parsed_components[0]]
                product_details = [self._resolve_molecule_by_abbr(abbr) for abbr in parsed_components[1]]

                if None in reactant_details or None in product_details:
                    self.logger.error(
                        f"Could not resolve all abbreviations for manual reaction row {idx + 1}. Skipping.")
                    continue

                # 3. Reconstruct the full reaction text using production names
                reactant_names = " + ".join([info['molecule_name'] for info in reactant_details])
                product_names = " + ".join([info['molecule_name'] for info in product_details])
                arrow = " <=> " if is_reversible else " => "
                full_reaction_text = f"{reactant_names}{arrow}{product_names}"

                # 4. Construct a "virtual" entry dictionary that matches the staging table format
                virtual_entry = {
                    "SwissLipids ID": None,  # Manual entries don't have an SLP ID
                    "UniProtKB AC(s)": row.get('Enzyme_UniProtKB_AC'),
                    "Reaction text": full_reaction_text,  # Use the newly constructed text
                    "Rhea ID": row.get('Rhea_ID'),
                    "Protein taxon": row.get('Organism_ID'),
                    "Taxon scientific name": row.get('Organism_ScientificName'),
                    "Gene name": row.get('Enzyme_GeneName'),
                    "DOI": row.get('DOI')  # Add DOI for the ensure_reaction call
                }

                # 5. Call the unified helper method with the virtual entry
                self.logger.info(f"Processing manual reaction row {idx + 1} as virtual entry...")
                self._process_and_link_reaction_entry(virtual_entry)

                successful_additions += 1

            except Exception as e:
                self.logger.error(f"Critical error on manual reaction row {idx + 1}: {e}", exc_info=True)
                raise

        self.logger.info(f"Finished processing {successful_additions}/{len(manual_df)} manual reactions.")
        return successful_additions

    def _archive_current_staging_data(self, run_timestamp_for_snapshot_naming: str):
        """Saves current staging data to snapshot TSV files."""
        self.logger.info(
            f"Archiving current staging data as snapshot for timestamp: {run_timestamp_for_snapshot_naming}...")
        snapshot_dir = self.output_dir / "snapshots"
        snapshot_dir.mkdir(parents=True, exist_ok=True)

        for table_name_config_key, snapshot_file_name_base in [
            ('get_staging_lipids_table', "lipids_snapshot"),
            ('get_staging_enzymes_table', "enzymes_snapshot")
        ]:
            try:
                table_name = getattr(self.config, table_name_config_key)()
                df = self._query_to_df(f"SELECT * FROM {table_name};")  # Use semicolon for safety
                if not df.empty:
                    # Use a fixed name for "previous" snapshot to simplify lookup next time
                    # Or use a timestamped archive for history (more complex lookup)
                    # For now, fixed name:
                    output_path = snapshot_dir / f"previous_{snapshot_file_name_base}.tsv"
                    # output_path = snapshot_dir / f"{snapshot_file_name_base}_{run_timestamp_for_snapshot_naming}.tsv" # Alternative: dated snapshots

                    df.to_csv(output_path, sep='\t', index=False, encoding='utf-8')
                    self.logger.info(f"Archived '{table_name}' data to '{output_path}'")
                else:
                    self.logger.info(f"Staging table '{table_name}' was empty. No snapshot archived for it.")
            except Exception as e:
                self.logger.error(f"Failed to archive staging table data for config key '{table_name_config_key}': {e}",
                                  exc_info=True)

    def apply_all_approved_updates(self, approved_delta_dir_path: Path,
                                   run_timestamp_of_deltas: str,
                                   commit_changes: bool = True):  # Added commit_changes for dry-run
        """
        Applies all approved new entities, field updates, and manual additions.
        Manages a single transaction for all database modifications.
        Archives staging data upon successful commit.

        Args:
            approved_delta_dir_path (Path): Directory containing 'approved_delta_*.csv'
                                            and 'approved_manual_add_*.tsv' files.
            run_timestamp_of_deltas (str): The timestamp string associated with the delta files
                                           (used to find correct approved files).
            commit_changes (bool): If True, commits changes. If False, rolls back (dry run).
        """
        self.logger.info(f"Starting to apply all approved updates. Commit mode: {commit_changes}")
        self.logger.info(
            f"Looking for approved files in: {approved_delta_dir_path} with timestamp affix: {run_timestamp_of_deltas}")

        # Define paths to approved files
        # Assuming reviewers rename/prepare files like: approved_delta_new_organisms_YYYYMMDD_HHMMSS_run.csv
        # The run_timestamp_of_deltas helps pick the correct set.
        approved_new_org_file = approved_delta_dir_path / f"approved_delta_new_organisms.tsv"  # Simpler: assume fixed name after approval
        approved_new_mol_file = approved_delta_dir_path / f"approved_delta_new_molecules.tsv"
        approved_new_rxn_file = approved_delta_dir_path / f"approved_delta_new_reaction_definitions.tsv"
        approved_new_enz_file = approved_delta_dir_path / f"approved_delta_new_enzyme_definitions.tsv"
        approved_new_links_file = approved_delta_dir_path / "approved_delta_new_reaction_links.tsv"
        approved_updates_file = approved_delta_dir_path / f"approved_delta_updated_entries.tsv"
        approved_manual_rxn_file = approved_delta_dir_path / f"approved_template_manual_add_reactions.tsv"  # Reviewer renames template

        if self.conn is None or self.conn.closed:
            raise ConnectionError("Database connection is not active for applying updates.")

        try:
            # Start transaction (psycopg3 connections start in no-transaction state,
            # first DML/DDL starts one, or explicitly with `with self.conn.transaction():`)
            # self.conn.autocommit = False # Ensure we are not in autocommit mode. Default is False.

            self.new_additions_summary.clear()
            self.abbreviation_warnings.clear()  # ## NEW: Clear warnings at start of apply ##

            with self.conn.transaction():  # psycopg3 preferred way for transactions
                self.logger.info("Transaction started for applying updates.")

                # 1. Apply approved new entities from SwissLipids diffs
                self._apply_approved_new_organisms(approved_new_org_file)
                self._apply_approved_new_molecules(approved_new_mol_file)
                self._apply_approved_new_reaction_definitions(approved_new_rxn_file)
                self._apply_approved_new_enzyme_definitions(approved_new_enz_file)
                self._apply_approved_reaction_links(approved_new_links_file)

                # 2. Apply approved updates to existing entities
                self._apply_approved_field_updates(approved_updates_file)

                # 3. Apply approved manual additions
                self._apply_manual_reactions(approved_manual_rxn_file)
                # Add calls for other manual additions (enzymes, molecules) if templates exist

                # 4. Create/update reaction_enzyme and reaction_pairs links
                #    This method needs to be robust. It processes current staging data
                #    to link entities that are now confirmed to be in production.
                # self._create_reaction_enzyme_and_pairs()  # Needs careful implementation

                if not commit_changes:
                    raise DryRunRollbackException("Dry run requested, changes will be rolled back.")

                self.logger.info("All updates applied successfully within transaction.")

            # This block executes ONLY if the transaction was not rolled back
            if commit_changes:
                self.logger.info("Transaction committed successfully.")
                # Archive staging data only after a successful commit
                self._archive_current_staging_data(run_timestamp_of_deltas)

        except DryRunRollbackException as dre:
            self.logger.info(f"Dry run completed successfully: {dre}")
        except Exception as e:
            self.logger.error(f"Failed to apply updates, transaction was rolled back: {e}", exc_info=True)
            raise
        finally:
            # Generate reports regardless of commit/rollback, as they reflect the work done
            self.logger.info("Generating final reports for the 'apply' run...")
            self._generate_new_additions_report(approved_delta_dir_path)
            self._generate_abbreviation_warning_report(approved_delta_dir_path)

    @lru_cache(maxsize=8192)  # Increased cache size slightly
    def _get_molecule_details_from_slm_id(self, slm_id: str) -> Optional[Dict[str, Any]]:
        """
        ## RENAMED & REFACTORED ##
        Helper to get the full molecule record from a SwissLipids ID.
        """
        if not slm_id: return None
        query = "SELECT * FROM lipograph.molecules WHERE swisslipids_id = %s LIMIT 1;"

        # Using dict_row factory ensures the result is a dict-like object
        temp_cursor = self.conn.cursor(row_factory=dict_row)
        temp_cursor.execute(query, (slm_id,))
        result = temp_cursor.fetchone()
        temp_cursor.close()

        if result:
            return dict(result)
        self.logger.warning(f"SLM ID '{slm_id}' from resolutions file not found in the production 'molecules' table.")
        return None

    @lru_cache(maxsize=8192)
    def _get_molecule_details(self, molecule_id: int) -> Optional[Dict[str, Any]]:
        """
        Fetches the full details for a single molecule from the production table
        using its internal primary key (molecule_id). The result is cached.

        Args:
            molecule_id (int): The primary key of the molecule to fetch.

        Returns:
            Optional[Dict[str, Any]]: A dictionary of the molecule's data or None if not found.
        """
        if molecule_id is None:
            return None

        self.logger.debug(f"Fetching details for molecule_id: {molecule_id}")
        query = """
            SELECT * 
            FROM lipograph.molecules 
            WHERE molecule_id = %s;
        """

        # Use a temporary cursor with dict_row factory for this self-contained lookup
        try:
            temp_cursor = self.conn.cursor(row_factory=dict_row)
            temp_cursor.execute(query, (molecule_id,))
            result = temp_cursor.fetchone()
            temp_cursor.close()

            if result:
                return dict(result)
            else:
                self.logger.warning(
                    f"Could not find details for molecule_id {molecule_id}, which was expected to exist.")
                return None
        except Exception as e:
            self.logger.error(f"Error fetching details for molecule_id {molecule_id}: {e}", exc_info=True)
            return None

    def _resolve_molecule_with_all_fallbacks(
            self,
            name: str,
            slp_id: Optional[str],
            rhea_id: Optional[str],
    ) -> Optional[Dict[str, Any]]:
        """
        Resolves a molecule name to a production molecule_id using the full resolution hierarchy.

        This is the central resolver function called for each potential lipid name from a
        SwissLipids reaction text. It follows a strict order of checks:

        1.  **Context-Specific YAML Cache:** Checks for a precise resolution rule you have
            manually created for this exact name in this specific reaction context
            (name + SLP ID + Rhea ID). This is the highest priority.

        2.  **Global Name YAML Cache:** Checks for a general resolution rule for this name
            that applies to all contexts. This is used if no context-specific rule exists.

        3.  **Local Database Resolution:** If not found in the YAML cache, it attempts to
            resolve the name using existing database information (exact name match,
            synonym match, or cleaned abbreviation match). This is the primary automated
            method for known lipids.

        4.  **Log for Manual Curation:** If all above methods fail, the name is considered
            unresolved. It is logged with its context, and its count is tracked in
            `self.unresolved_names_this_run`. This method does *not* make live API calls;
            that task is delegated to the separate interactive curation script.

        Args:
            name (str): The molecule name to resolve (from reaction text).
            slp_id (Optional[str]): The SwissLipids enzyme/process ID for context.
            rhea_id (Optional[str]): The Rhea ID for context.

        Returns:
            Optional[int]: The resolved `molecule_id` from the production database,
                           or `None` if the name could not be resolved or was explicitly
                           marked as "no match" in the YAML file.
        """
        if not name or not name.strip():
            self.logger.warning("[Resolver] Received an empty name for resolution.")
            return None

        # --- Step 1 & 2: Check Curated Resolutions from YAML file ---
        rhea_for_key = rhea_id if (rhea_id and str(rhea_id).strip()) else "N/A"
        context_key = f"{name}::{slp_id}::{rhea_for_key}"

        # Check for a highly specific, context-aware rule first.
        if context_key in self.context_resolutions:
            resolution = self.context_resolutions[context_key]
            self.logger.debug(
                f"[Resolver] YAML-CONTEXT: Found rule for '{name}' in context [{slp_id}/{rhea_for_key}] -> '{resolution}'")
            if resolution == NO_MATCH_MARKER:
                return None  # Explicitly marked as not a match in this context.

            # We have an SLM ID, now get the internal molecule_id.
            return self._get_molecule_details_from_slm_id(resolution)

        # If no context-specific rule, check for a global rule for this name.
        if name in self.global_name_resolutions:
            resolution = self.global_name_resolutions[name]
            self.logger.debug(f"[Resolver] YAML-GLOBAL: Found rule for '{name}' -> '{resolution}'")
            if resolution == NO_MATCH_MARKER:
                return None  # Explicitly marked as not a match globally.

            return self._get_molecule_details_from_slm_id(resolution)

        # Step 3: Local Database Resolution
        mol_details = self._resolve_molecule_name_from_sl_text(name)
        if mol_details is not None:
            return mol_details

        # --- Step 4: Fallback - Log for Manual Curation ---
        # If we've reached this point, the name is unresolved by all automated and cached means.
        # We will not attempt a live API call here. Instead, we log it clearly
        # so the interactive script can pick it up later.
        self.logger.warning(
            f"[Resolver] UNRESOLVED: Could not resolve '{name}' using YAML or local DB. "
            f"(Context: SLP {slp_id}, Rhea {rhea_for_key}). It needs to be curated."
        )
        # Track this failure for the final summary report.
        self.unresolved_names_this_run[name] += 1

        return None

    def _process_and_link_reaction_entry(self, reaction_entry: Dict[str, Any]):
        """
        A unified helper that takes a single reaction entry (either from staging or a manual
        file) and performs all the necessary ensuring, parsing, resolving, and linking.

        This is the core method that handles the entire linking process for a single reaction:
        1. Reads reaction-enzyme associations from the entry.
        2. Ensures the corresponding organism, reaction, and enzyme entities exist.
        3. Parses the reaction text, filters out common non-lipid metabolites.
        4. Resolves lipid names to production molecule_ids using fallback strategies.
        5. If all lipid participants are successfully resolved, creates reaction_enzyme
           links and corresponding reaction_pairs.
        """
        # --- Extract all necessary context from the entry ---
        slp_id = reaction_entry.get("SwissLipids ID")
        uniprot_acs_str = reaction_entry.get("UniProtKB AC(s)", "")
        reaction_text = reaction_entry.get("Reaction text", "")
        rhea_id_str = reaction_entry.get("Rhea ID")
        organism_id_str = reaction_entry.get("Protein taxon", "")
        organism_name = reaction_entry.get("Taxon scientific name")
        gene_name = reaction_entry.get("Gene name")
        doi = reaction_entry.get("DOI")  # For manual entries

        self.logger.info(f"[Linker] Processing entry: SLP {slp_id or 'N/A'}, Rhea {rhea_id_str or 'N/A'}")

        if not all([uniprot_acs_str, reaction_text, organism_id_str]):
            self.logger.warning(f"[Linker] Skipping entry due to missing essential data: {reaction_entry}")
            return

        try:
            # --- 1. Ensure Organism and Reaction entities exist in production ---
            organism_id_prod = self._ensure_organism(organism_id_str, organism_name)
            if organism_id_prod is None:
                self.logger.error(
                    f"[Linker] Failed to ensure organism for TaxonID '{organism_id_str}'. Skipping entry.")
                return

            prod_reaction_id = self._ensure_reaction(reaction_text, rhea_id_str, doi)
            if prod_reaction_id is None:
                self.logger.error(f"[Linker] Failed to ensure reaction for text '{reaction_text}'. Skipping entry.")
                return

            # --- 2. Parse reaction text and filter metabolites ---
            parsed_components_tuple = split_reaction_text(reaction_text)
            if parsed_components_tuple[0] is None:
                self.logger.warning(f"[Linker] Could not parse reaction text: '{reaction_text}'. Skipping.")
                return

            sl_reactant_names, sl_product_names = parsed_components_tuple[0]
            sl_is_reversible = parsed_components_tuple[1]

            ignored_metabolites_set = set(name.lower() for name in self.config.get_ignored_metabolites_for_pairs())
            reactants_to_resolve = [name for name in sl_reactant_names if
                                    name.strip().lower() not in ignored_metabolites_set]
            products_to_resolve = [name for name in sl_product_names if
                                   name.strip().lower() not in ignored_metabolites_set]

            if not reactants_to_resolve or not products_to_resolve:
                self.logger.info(
                    f"[Linker] Reaction '{reaction_text}' has no lipid reactants/products after filtering. Skipping.")
                return

            # --- 3. Resolve all lipid participants ---
            reactant_details_list = [
                self._resolve_molecule_with_all_fallbacks(name, slp_id, rhea_id_str) for name in reactants_to_resolve
            ]
            product_details_list = [
                self._resolve_molecule_with_all_fallbacks(name, slp_id, rhea_id_str) for name in products_to_resolve
            ]

            if None in reactant_details_list or None in product_details_list:
                self.logger.warning(f"[Linker] Skipping pair creation for '{reaction_text}' due to unresolved lipids.")
                return

            # --- 4. Validate abbreviations for ALL resolved lipids (Reactants and Products) ---
            all_resolved_lipids_details = reactant_details_list + product_details_list
            for lipid_details in all_resolved_lipids_details:
                abbr = lipid_details.get('abbreviation')
                cleaned_abbr = lipid_details.get('cleaned_abbreviation')

                is_missing_abbr = not abbr or not str(abbr).strip()
                is_missing_cleaned_abbr = not cleaned_abbr or not str(cleaned_abbr).strip()

                if is_missing_abbr or is_missing_cleaned_abbr:
                    warning_msg = (
                        f"ABBREVIATION WARNING (Reaction ID: {prod_reaction_id}): "
                        f"Lipid participant '{lipid_details['molecule_name']}' (ID: {lipid_details['molecule_id']}) is missing a critical abbreviation. "
                        f"[Abbr: {'MISSING' if is_missing_abbr else 'OK'}, "
                        f"Cleaned Abbr: {'MISSING' if is_missing_cleaned_abbr else 'OK'}]"
                    )
                    self.logger.warning(warning_msg)
                    self.abbreviation_warnings.append({
                        'reaction_id': prod_reaction_id,
                        'reaction_text': reaction_text,
                        'molecule_id': lipid_details['molecule_id'],
                        'molecule_name': lipid_details['molecule_name'],
                        'molecule_slm_id': lipid_details.get('swisslipids_id', 'N/A')
                    })

            # --- 5. If all lipids resolved, proceed to link enzymes and create pairs ---
            reactant_ids = [d['molecule_id'] for d in reactant_details_list]
            product_ids = [d['molecule_id'] for d in product_details_list]

            uniprot_acs = [ac.strip() for ac in uniprot_acs_str.split('|') if ac.strip()]
            for ac in uniprot_acs:
                # Ensure the enzyme exists for this specific organism
                prod_enzyme_id = self._ensure_enzyme(ac, gene_name, organism_id_prod, swisslipids_p_id=slp_id)
                if prod_enzyme_id is None:
                    self.logger.error(
                        f"[Linker] Failed to ensure enzyme for UniProt AC {ac}. Skipping this enzyme link.")
                    continue

                # Ensure the link between the reaction and the enzyme exists
                prod_reaction_enzyme_id = self._ensure_reaction_enzyme(prod_reaction_id, prod_enzyme_id)
                if prod_reaction_enzyme_id is None:
                    self.logger.error(
                        f"[Linker] Failed to ensure reaction_enzyme link for R_ID:{prod_reaction_id}, E_ID:{prod_enzyme_id}. Skipping pairs for this link.")
                    continue

                # Create the specific reactant->product pairs for this reaction-enzyme instance
                self._ensure_reaction_pairs(prod_reaction_enzyme_id, reactant_ids, product_ids, sl_is_reversible)

                self.logger.debug(
                    f"[Linker] Successfully processed links for Reaction ID {prod_reaction_id} and Enzyme ID {prod_enzyme_id} (RE_ID: {prod_reaction_enzyme_id}).")

        except Exception as e:
            self.logger.error(f"[Linker] A critical error occurred while processing entry: {reaction_entry}",
                              exc_info=True)
            raise  # Propagate the error to ensure the entire transaction is rolled back

class DryRunRollbackException(Exception):
    pass
