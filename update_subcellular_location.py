import argparse
import logging
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Set, Optional, Dict

import pandas as pd
import requests

# --- Path Setup ---
try:
    script_dir = Path(__file__).resolve().parent
    src_path = script_dir / "src"
    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))
except Exception as e:
    print(f"Error setting up sys.path: {e}", file=sys.stderr)

from src.data.config_loader import ConfigLoader
from src.utils.db_updater import LipidDatabaseUpdater

# --- Global Configuration ---
logger = logging.getLogger(__name__)

BATCH_SIZE = 300  # UniProt ID Mapper can handle larger batches than SPARQL


class SubcellularLocationUpdater:
    """
    Updates the enzymes table with subcellular location data from UniProt ID Mapper.
    """

    def __init__(self, config_loader: ConfigLoader):
        self.config_loader = config_loader
        self.db_updater = LipidDatabaseUpdater(config_loader)
        logger.info("SubcellularLocationUpdater initialized.")

    def _get_uniprot_subcellular_locations(self, uniprot_ids: Set[str]) -> Optional[pd.DataFrame]:
        """
        Uses UniProt ID Mapper REST API to get subcellular location data.
        """
        if not uniprot_ids:
            return None

        uniprot_ids_list = list(uniprot_ids)

        # Step 1: Submit the job
        submit_data = {
            "from": "UniProtKB_AC-ID",
            "to": "UniProtKB-Swiss-Prot",
            "ids": ",".join(uniprot_ids_list)
        }

        logger.info(f"Submitting ID mapping job for {len(uniprot_ids)} UniProt IDs...")

        try:
            # Submit the job
            response = requests.post(
                "https://rest.uniprot.org/idmapping/run",
                data=submit_data
            )
            response.raise_for_status()
            job_id = response.json()["jobId"]
            logger.info(f"Job submitted with ID: {job_id}")

            # Step 2: Poll for completion
            max_attempts = 60  # Wait up to 5 minutes
            polling_interval = 5

            for attempt in range(max_attempts):
                status_response = requests.get(f"https://rest.uniprot.org/idmapping/status/{job_id}")
                status_response.raise_for_status()

                status_data = status_response.json()

                # Check if job is finished
                if "results" in status_data or "failedIds" in status_data:
                    logger.info(f"Job completed after {attempt + 1} attempts")
                    break
                elif "jobStatus" in status_data:
                    if status_data["jobStatus"] in ("NEW", "RUNNING"):
                        logger.debug(f"Job status: {status_data['jobStatus']}, retrying in {polling_interval}s")
                        time.sleep(polling_interval)
                    else:
                        logger.error(f"Job failed with status: {status_data['jobStatus']}")
                        return None
                else:
                    time.sleep(polling_interval)
            else:
                logger.error("UniProt ID mapping job timed out")
                return None

            # Step 3: Get the results URL
            details_response = requests.get(f"https://rest.uniprot.org/idmapping/details/{job_id}")
            details_response.raise_for_status()
            redirect_url = details_response.json()["redirectURL"]

            # Convert to stream URL and add our desired columns and format
            stream_url = redirect_url.replace("/results/", "/results/stream/")
            stream_url += "?format=tsv&fields=accession,cc_subcellular_location"

            # Step 4: Download results
            logger.info("Downloading subcellular location results...")
            results_response = requests.get(stream_url)
            results_response.raise_for_status()

            # Step 5: Parse into DataFrame
            from io import StringIO
            df = pd.read_csv(StringIO(results_response.text), sep='\t', dtype=str)
            logger.info(f"Retrieved subcellular location data for {len(df)} UniProt entries")

            return df

        except requests.exceptions.RequestException as e:
            logger.error(f"UniProt ID Mapper request failed: {e}")
            return None
        except KeyError as e:
            logger.error(f"Unexpected response structure from UniProt API: {e}")
            return None
        except Exception as e:
            logger.error(f"Error processing UniProt ID Mapper results: {e}")
            return None

    def add_subcellular_location_column(self, commit: bool = False) -> bool:
        """
        Adds subcellular_location column to enzymes table if it doesn't exist.
        """
        logger.info("Checking if subcellular_location column exists...")

        with self.db_updater as updater:
            # Set autocommit for direct execution
            updater.conn.autocommit = True

            try:
                # Check if column already exists
                column_check = updater._query_to_df(
                    """
                    SELECT column_name 
                    FROM information_schema.columns 
                    WHERE table_schema = 'lipograph' 
                    AND table_name = 'enzymes' 
                    AND column_name = 'subcellular_location';
                    """
                )

                if not column_check.empty:
                    logger.info("subcellular_location column already exists.")
                    return True

                logger.info("Adding subcellular_location column to enzymes table...")
                if not commit:
                    logger.warning("DRY RUN: Would add subcellular_location column.")
                    return True

                # Add the column directly with autocommit
                cursor = updater.conn.cursor()
                cursor.execute(
                    "ALTER TABLE lipograph.enzymes ADD COLUMN subcellular_location TEXT;"
                )
                logger.info("Successfully added subcellular_location column.")
                return True

            except Exception as e:
                logger.error(f"Failed to add subcellular_location column: {e}")
                return False

    def update_subcellular_locations(self, commit: bool = False) -> None:
        """
        Updates all enzymes with subcellular location data from UniProt.
        """
        logger.info("--- Updating Subcellular Locations ---")

        if not commit:
            logger.warning("DRY RUN MODE: No changes will be committed to the database.")

        # First, ensure the column exists
        if not self.add_subcellular_location_column(commit):
            logger.error("Could not add subcellular_location column. Aborting.")
            return

        with self.db_updater as updater:
            # Set autocommit for direct execution
            updater.conn.autocommit = True

            try:
                # Get all UniProt IDs from enzymes table
                logger.info("Fetching all UniProt IDs from enzymes table...")
                enzymes_df = updater._query_to_df(
                    "SELECT enzyme_id, uniprot_id FROM lipograph.enzymes WHERE uniprot_id IS NOT NULL and subcellular_location is NULL;"
                )

                if enzymes_df.empty:
                    logger.warning("No enzymes with UniProt IDs found.")
                    return

                logger.info(f"Found {len(enzymes_df)} enzymes with UniProt IDs")

                # Get unique UniProt IDs
                unique_uniprot_ids = set(enzymes_df['uniprot_id'].unique())
                logger.info(f"Processing {len(unique_uniprot_ids)} unique UniProt IDs")

                # Process in batches
                uniprot_ids_list = list(unique_uniprot_ids)
                location_data = {}

                for i in range(0, len(uniprot_ids_list), BATCH_SIZE):
                    batch_ids = set(uniprot_ids_list[i:i + BATCH_SIZE])
                    batch_num = i // BATCH_SIZE + 1
                    total_batches = (len(uniprot_ids_list) + BATCH_SIZE - 1) // BATCH_SIZE

                    logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch_ids)} IDs)...")

                    batch_df = self._get_uniprot_subcellular_locations(batch_ids)
                    if batch_df is not None:
                        for _, row in batch_df.iterrows():
                            uniprot_id = row.get('Entry', '')
                            location = row.get('Subcellular location [CC]', '')

                            if uniprot_id and pd.notna(location) and location.strip():
                                location_data[uniprot_id] = location.strip()

                    # Small delay between batches to be nice to UniProt
                    if i + BATCH_SIZE < len(uniprot_ids_list):
                        time.sleep(1)

                logger.info(f"Retrieved subcellular location data for {len(location_data)} proteins")

                if not location_data:
                    logger.warning("No subcellular location data retrieved. Nothing to update.")
                    return

                # Update the database
                if not commit:
                    logger.info("DRY RUN: Would update the following enzymes:")
                    for uniprot_id, location in list(location_data.items())[:10]:  # Show first 10
                        logger.info(f"  {uniprot_id}: {location[:100]}{'...' if len(location) > 100 else ''}")
                    if len(location_data) > 10:
                        logger.info(f"  ... and {len(location_data) - 10} more")
                    return

                # Actual database update
                logger.info("Updating enzymes table with subcellular location data...")
                update_count = 0

                cursor = updater.conn.cursor()
                for uniprot_id, location in location_data.items():
                    cursor.execute(
                        """
                        UPDATE lipograph.enzymes 
                        SET subcellular_location = %s 
                        WHERE uniprot_id = %s AND (subcellular_location IS NULL OR subcellular_location = '');
                        """,
                        (location, uniprot_id)
                    )
                    if cursor.rowcount > 0:
                        update_count += 1

                logger.info(f"Successfully updated {update_count} enzymes with subcellular location data.")

            except Exception as e:
                logger.error(f"Error during subcellular location update: {e}", exc_info=True)

    def show_statistics(self) -> None:
        """
        Shows statistics about subcellular location data coverage.
        """
        logger.info("--- Subcellular Location Statistics ---")

        with self.db_updater as updater:
            try:
                # Check if column exists
                column_check = updater._query_to_df(
                    """
                    SELECT column_name 
                    FROM information_schema.columns 
                    WHERE table_schema = 'lipograph' 
                    AND table_name = 'enzymes' 
                    AND column_name = 'subcellular_location';
                    """
                )

                if column_check.empty:
                    logger.info("subcellular_location column does not exist yet.")
                    return

                # Get statistics
                stats_df = updater._query_to_df(
                    """
                    SELECT 
                        COUNT(*) as total_enzymes,
                        COUNT(subcellular_location) as enzymes_with_location,
                        COUNT(*) - COUNT(subcellular_location) as enzymes_without_location,
                        ROUND(100.0 * COUNT(subcellular_location) / COUNT(*), 1) as coverage_percentage
                    FROM lipograph.enzymes 
                    WHERE uniprot_id IS NOT NULL;
                    """
                )

                if not stats_df.empty:
                    row = stats_df.iloc[0]
                    logger.info(f"Total enzymes with UniProt IDs: {row['total_enzymes']}")
                    logger.info(f"Enzymes with subcellular location: {row['enzymes_with_location']}")
                    logger.info(f"Enzymes without subcellular location: {row['enzymes_without_location']}")
                    logger.info(f"Coverage: {row['coverage_percentage']}%")

                # Show some example locations
                examples_df = updater._query_to_df(
                    """
                    SELECT uniprot_id, enzyme_name, subcellular_location
                    FROM lipograph.enzymes 
                    WHERE subcellular_location IS NOT NULL
                    LIMIT 5;
                    """
                )

                if not examples_df.empty:
                    logger.info("\nExample subcellular locations:")
                    for _, row in examples_df.iterrows():
                        location = row['subcellular_location']
                        display_location = location[:100] + "..." if len(location) > 100 else location
                        logger.info(f"  {row['uniprot_id']} ({row['enzyme_name']}): {display_location}")

            except Exception as e:
                logger.error(f"Error getting statistics: {e}")


def main():
    """Main function to run the CLI."""
    parser = argparse.ArgumentParser(
        description="Update enzymes table with subcellular location data from UniProt.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--config-file",
        type=str,
        default="config/lipid_config.yaml",
        help="Path to the configuration file.",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument(
        "--action",
        choices=["update", "stats", "add-column"],
        required=True,
        help="Action to perform: update subcellular locations, show stats, or just add column.",
    )
    parser.add_argument(
        "--commit",
        action="store_true",
        help="Commit changes to the database. Default is a dry run.",
    )

    args = parser.parse_args()

    # Setup logging
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    log_file_path = log_dir / f"subcellular_location_{args.action}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    logging.basicConfig(
        level=args.log_level.upper(),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file_path, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )

    try:
        config = ConfigLoader(args.config_file)
        updater = SubcellularLocationUpdater(config)

        if args.action == "add-column":
            updater.add_subcellular_location_column(args.commit)
        elif args.action == "update":
            updater.update_subcellular_locations(args.commit)
        elif args.action == "stats":
            updater.show_statistics()

    except FileNotFoundError as e:
        logger.critical(f"A required file was not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.critical(f"An unexpected error occurred: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()