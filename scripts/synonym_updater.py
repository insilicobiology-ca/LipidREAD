import logging
from typing import List, Dict, Optional, Literal, Set
from src.data.data_access import DatabaseLipidDataAccess
from src.data.config_loader import ConfigLoader
from src.parsing.lipid_parser import LipidParser
from src.utils.lipid_translator import LipidTranslator
from src.utils.lipid_translator import LipidTranslator

logger = logging.getLogger(__name__)
lipid_parser = LipidParser()


class LipidSynonymUpdater:
    def __init__(
            self,
            data_access: DatabaseLipidDataAccess,
            lipid_parser: LipidParser,
            lipid_translator: LipidTranslator,
    ):
        self.data_access = data_access
        self.lipid_parser = lipid_parser
        self.lipid_translator = lipid_translator



    def _process_lipid_name(self, lipid_name: str, prefix: Literal['Glc', 'Gal']) -> Optional[str]:
        """
        Process a lipid name to generate its Hex variant
        """
        try:
            lipid_adduct = lipid_parser.parse(lipid_name)
            lipid_species_adduct = lipid_adduct.get_lipid_string(LipidLevel.SPECIES)

            print(f"Lipid name: {lipid_name} -> {LipidTranslator.clean_lipid_abbreviation(lipid_name)}")

            if prefix in lipid_species_adduct:
                return lipid_species_adduct.replace(prefix, 'Hex')
            return None

        except Exception as e:
            logger.error(f"Failed to process lipid {lipid_name}: {str(e)}")
            return None

    def add_synonym(self, lipid_id: str, new_synonym: str) -> bool:
        """
        Add a new synonym to the SYNONYMS field.
        Returns True if successful, False otherwise.
        """
        try:
            with self.conn:  # Uses transaction
                self.cursor.execute(
                    """
                    UPDATE "LipidTranslator".test_lipids 
                    SET "SYNONYMS" = CASE 
                        WHEN "SYNONYMS" IS NULL OR "SYNONYMS" = '' THEN %s
                        ELSE "SYNONYMS" || %s
                    END
                    WHERE "LM_ID" = %s
                    RETURNING "LM_ID"
                    """,
                    (f"|{new_synonym}|", new_synonym + "|", lipid_id)
                )
                return bool(self.cursor.fetchone())
        except psycopg2.Error as e:
            logger.error(f"Failed to add synonym for lipid {lipid_id}: {str(e)}")
            return False

    def process_ceramides(self, prefix: Literal['Glc', 'Gal'], batch_size: int = 200) -> Dict[str, int]:
        """
        Process ceramides with specific prefix (Glc or Gal) to generate Hex synonyms.
        Also track cases where processed name differs from ABBREVIATION column.
        """
        stats = {
            'processed': 0,
            'successful': 0,
            'failed': 0,
            'skipped': 0,
            'abbreviation_mismatches': 0
        }

        try:
            self.cursor.execute(
                """
                SELECT COUNT(*) FROM "LipidTranslator".test_lipids 
                WHERE "NAME" LIKE %s
                """,
                (f'{prefix}Cer%',)
            )
            total_rows = self.cursor.fetchone()[0]

            self.cursor.execute(
                """
                SELECT "LM_ID", "NAME", "SYNONYMS", "ABBREVIATION"
                FROM "LipidTranslator".test_lipids 
                WHERE "NAME" LIKE %s
                """,
                (f'{prefix}Cer%',)
            )

            while True:
                rows = self.cursor.fetchmany(batch_size)
                if not rows:
                    break

                for row in rows:
                    stats['processed'] += 1
                    try:
                        new_synonym = self._process_lipid_name(row['NAME'], prefix)

                        if not new_synonym:
                            stats['skipped'] += 1
                            continue

                        # Compare with ABBREVIATION column
                        if row['ABBREVIATION'] and new_synonym != row['ABBREVIATION']:
                            stats['abbreviation_mismatches'] += 1
                            self.abbreviation_discrepancies.add((
                                row['LM_ID'],
                                row['NAME'],
                                new_synonym,  # goslin's output
                                row['ABBREVIATION']  # current abbreviation in DB
                            ))
                            logger.warning(
                                f"Abbreviation mismatch for {row['LM_ID']}:\n"
                                f"  Goslin output: {new_synonym}\n"
                                f"  Current ABBREVIATION: {row['ABBREVIATION']}"
                            )

                        # Skip if new_synonym is already in SYNONYMS
                        if row['SYNONYMS'] and f"|{new_synonym}|" in row['SYNONYMS']:
                            stats['skipped'] += 1
                            continue

                        if self.add_synonym(row['LM_ID'], new_synonym):
                            stats['successful'] += 1
                        else:
                            stats['failed'] += 1

                    except Exception as e:
                        logger.error(f"Error processing lipid {row['LM_ID']}: {str(e)}")
                        stats['failed'] += 1

                logger.info(f"Progress: {stats['processed']}/{total_rows} rows")

        except Exception as e:
            logger.error(f"Failed to process {prefix}Cer lipids: {str(e)}")
            raise

        return stats

    def save_discrepancies(self, filename: str):
        """
        Save abbreviation discrepancies to a file
        """
        try:
            with open(filename, 'w') as f:
                f.write("LM_ID,NAME,GOSLIN_OUTPUT,CURRENT_ABBREVIATION\n")
                for entry in sorted(self.abbreviation_discrepancies):
                    f.write(','.join(str(x) for x in entry) + '\n')
            logger.info(f"Saved {len(self.abbreviation_discrepancies)} discrepancies to {filename}")
        except Exception as e:
            logger.error(f"Failed to save discrepancies to file: {str(e)}")

    def __del__(self):
        if self.cursor:
            self.cursor.close()
        if self.conn:
            self.conn.close()


def update_lipid_synonyms(config_path: str, discrepancies_file: str = "abbreviation_discrepancies.csv"):
    """
    Update lipid synonyms for both GlcCer and GalCer entries
    """
    config = ConfigLoader(config_path)
    updater = LipidSynonymUpdater(config)

    try:
        # Process GlcCer entries
        glc_stats = updater.process_ceramides('Glc')
        logger.info(f"GlcCer processing complete. Stats: {glc_stats}")

        # Process GalCer entries
        gal_stats = updater.process_ceramides('Gal')
        logger.info(f"GalCer processing complete. Stats: {gal_stats}")

        updater.save_discrepancies(discrepancies_file)

    except Exception as e:
        logger.error(f"Failed to update synonyms: {str(e)}")


update_lipid_synonyms("/config/lipid_config.yaml")

