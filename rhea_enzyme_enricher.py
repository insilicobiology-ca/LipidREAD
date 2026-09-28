import argparse
import json
import logging
import sys
import time
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Set

import pandas as pd
import requests

from src.utils.text_utils import split_reaction_text

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

# Constants for filenames
DISCOVERED_DATA_FILE = "discovered_enzymes.json"
MAP_TEMPLATE_FILE = "uniprot_to_gene_map_template.tsv"
UNIPROT_DETAILS_FILE = "uniprot_details.tsv"
BATCH_SIZE = 40
SPARQL_TIMEOUT = 300  # 5-minute timeout for SPARQL queries


class RheaEnzymeEnricher:
    """
    Enriches the LipidCRED database by discovering and ingesting new enzymes
    for reactions with existing Rhea IDs, using UniProt ID Mapper and
    a two-phase interactive workflow.
    """

    def __init__(self, config_loader: ConfigLoader):
        self.config_loader = config_loader
        self.db_updater = LipidDatabaseUpdater(config_loader)
        self.uniprot_sparql_endpoint = "https://sparql.uniprot.org/sparql"
        logger.info("RheaEnzymeEnricher initialized.")

    def _execute_batched_sparql_query(
            self, rhea_ids_batch: List[int]
    ) -> Optional[List[Dict]]:
        """
        Executes a SPARQL query for a BATCH of MASTER Rhea IDs to get UniProt protein IDs.
        We'll use the ID Mapper for detailed information afterwards.
        """
        if not rhea_ids_batch:
            return []

        values_clause = " ".join(
            [f"<http://rdf.rhea-db.org/{rhea_id}>" for rhea_id in rhea_ids_batch]
        )

        query = f"""
        PREFIX up: <http://purl.uniprot.org/core/>

        SELECT DISTINCT ?reaction ?protein
        WHERE {{
          VALUES ?reaction {{ {values_clause} }}

          SERVICE <{self.uniprot_sparql_endpoint}> {{
            ?protein a up:Protein ;
                     up:reviewed true ;
                     up:annotation ?annotation .

            ?annotation a up:Catalytic_Activity_Annotation ;
                       up:catalyticActivity/up:catalyzedReaction ?reaction .
          }}
        }}
        """
        try:
            logger.debug(
                f"Executing SPARQL query for batch of {len(rhea_ids_batch)} master Rhea IDs."
            )
            response = requests.post(
                self.uniprot_sparql_endpoint,
                data={"query": query},
                headers={"Accept": "application/sparql-results+json"},
                timeout=SPARQL_TIMEOUT,
            )
            response.raise_for_status()
            results = response.json()["results"]["bindings"]
            logger.debug(f"SPARQL query for batch returned {len(results)} results.")
            return results
        except requests.exceptions.RequestException as e:
            logger.error(f"SPARQL request failed for batch: {e}")
        except json.JSONDecodeError as e:
            logger.error(f"Failed to decode JSON from SPARQL response for batch: {e}")
        return None

    def _get_uniprot_details_via_id_mapper(
            self, uniprot_ids: Set[str], output_path: Path
    ) -> Optional[pd.DataFrame]:
        """
        Uses UniProt ID Mapper REST API to get detailed information about proteins,
        including their catalytic activities and specific Rhea IDs.
        """
        if not uniprot_ids:
            return None

        uniprot_ids_list = list(uniprot_ids)

        # Step 1: Submit the job using the REST API
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
            stream_url += "?format=tsv&fields=accession,gene_primary,organism_id,rhea,cc_catalytic_activity"

            # Step 4: Download results
            logger.info("Downloading mapping results...")
            results_response = requests.get(stream_url)
            results_response.raise_for_status()

            # Step 5: Save to file and load into DataFrame
            details_file = output_path / UNIPROT_DETAILS_FILE
            with open(details_file, 'w', encoding='utf-8') as f:
                f.write(results_response.text)

            logger.info(f"UniProt details saved to: {details_file}")

            # Load into DataFrame
            df = pd.read_csv(details_file, sep='\t', dtype=str)
            logger.info(f"Loaded {len(df)} UniProt entries with detailed information")

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

    def _parse_catalytic_activity(self, catalytic_activity: str, target_rhea_ids: Set[int]) -> List[Tuple[int, str]]:
        """
        Parses the catalytic activity string to extract specific directional Rhea IDs
        that match our targets, along with their physiological directions.

        Returns list of tuples: (specific_rhea_id, direction)
        where direction is something like "left-to-right", "right-to-left", etc.
        """
        if pd.isna(catalytic_activity) or not catalytic_activity:
            return []

        results = []

        # Split by "CATALYTIC ACTIVITY:" to get individual reaction segments
        segments = catalytic_activity.split("CATALYTIC ACTIVITY:")

        for segment in segments:
            if not segment.strip():
                continue

            # Look for master Rhea ID in this segment
            master_match = re.search(r'Xref=Rhea:RHEA:(\d+)', segment)
            if not master_match:
                continue

            master_rhea_id = int(master_match.group(1))

            # Look for physiological direction and corresponding specific Rhea ID
            # Pattern: PhysiologicalDirection=direction; Xref=Rhea:RHEA:specific_id;
            direction_pattern = r'PhysiologicalDirection=([^;]+);\s*Xref=Rhea:RHEA:(\d+);'
            direction_match = re.search(direction_pattern, segment)

            if direction_match:
                direction = direction_match.group(1).strip()
                specific_rhea_id = int(direction_match.group(2))

                # Check if this specific Rhea ID is one we're interested in
                if specific_rhea_id in target_rhea_ids:
                    results.append((specific_rhea_id, direction))
                    logger.debug(
                        f"Found matching catalytic activity: Master={master_rhea_id}, Specific={specific_rhea_id}, Direction={direction}")

        return results

    def discover_enzymes(self, rhea_directions_path: str, output_dir: str) -> None:
        """
        Phase 1: Discovers new potential enzymes using master Rhea IDs and UniProt ID Mapper.
        """
        logger.info("--- Phase 1: Discovering New Enzymes (ID Mapper Approach) ---")
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        try:
            directions_df = pd.read_csv(rhea_directions_path, sep="\t", dtype=str)
            for col in directions_df.columns:
                directions_df[col] = directions_df[col].str.strip()
            logger.info(
                f"Loaded and cleaned {len(directions_df)} Rhea direction quartets."
            )
        except Exception as e:
            logger.critical(
                f"Failed to load or process Rhea directions file '{rhea_directions_path}': {e}"
            )
            return

        discovered_data = {}
        valid_reactions_to_query = []
        target_rhea_ids = set()  # All specific Rhea IDs we're interested in

        with self.db_updater as updater:
            local_reactions_df = updater._query_to_df(
                "SELECT reaction_id, reaction_text, rhea_id FROM lipograph.reactions WHERE rhea_id IS NOT NULL;"
            )
            logger.info(
                f"Found {len(local_reactions_df)} reactions with Rhea IDs in local DB for validation."
            )

            # Step 1: Validate reactions and collect master IDs and target specific IDs
            for _, row in local_reactions_df.iterrows():
                reaction_id, rhea_id = row["reaction_id"], row["rhea_id"]
                reaction_text = row.get("reaction_text", "")

                components, is_reversible = split_reaction_text(reaction_text)
                if not components:
                    logger.warning(
                        f"Could not parse reaction text for reaction_id {reaction_id}. Skipping."
                    )
                    continue

                if is_reversible:
                    direction_inferred = "bi"
                elif ' <= ' in reaction_text:
                    direction_inferred = "rl"  # right-to-left
                else:
                    direction_inferred = "lr"  # left-to-right (=> or default)
                rhea_id_str = str(rhea_id).strip()

                quartet_rows = directions_df[
                    directions_df.isin([rhea_id_str]).any(axis=1)
                ]
                if quartet_rows.empty:
                    logger.warning(
                        f"Rhea ID '{rhea_id_str}' not in directions file. Skipping reaction_id {reaction_id}."
                    )
                    continue

                quartet_row = quartet_rows.iloc[0]

                is_valid_direction = (
                                             direction_inferred == "lr" and quartet_row["RHEA_ID_LR"] == rhea_id_str
                                     ) or (
                                             direction_inferred == "rl" and quartet_row["RHEA_ID_RL"] == rhea_id_str
                                     # Add this line
                                     ) or (
                                             direction_inferred == "bi" and quartet_row["RHEA_ID_BI"] == rhea_id_str
                                     )

                if is_valid_direction:
                    master_rhea_id = int(quartet_row["RHEA_ID_MASTER"])
                    valid_reactions_to_query.append((reaction_id, master_rhea_id, rhea_id))
                    target_rhea_ids.add(rhea_id)  # Add the specific Rhea ID we care about
                    logger.debug(
                        f"Reaction {reaction_id} (Rhea:{rhea_id_str}) validated. Will query with master ID: {master_rhea_id}."
                    )

            logger.info(
                f"Found {len(valid_reactions_to_query)} reactions with valid directionality to query."
            )
            if not valid_reactions_to_query:
                return

            # Step 2: Get UniProt IDs using master Rhea IDs
            unique_master_ids = sorted(
                list(set(r[1] for r in valid_reactions_to_query))
            )

            all_uniprot_ids = set()
            master_id_to_reaction_ids_map = {}

            for reaction_id, master_id, specific_rhea_id in valid_reactions_to_query:
                if master_id not in master_id_to_reaction_ids_map:
                    master_id_to_reaction_ids_map[master_id] = []
                master_id_to_reaction_ids_map[master_id].append((reaction_id, specific_rhea_id))

            # Query SPARQL to get UniProt IDs
            num_master_ids = len(unique_master_ids)
            for i in range(0, num_master_ids, BATCH_SIZE):
                batch_master_ids = unique_master_ids[i: i + BATCH_SIZE]
                logger.info(
                    f"Processing SPARQL batch {i // BATCH_SIZE + 1}/{(num_master_ids + BATCH_SIZE - 1) // BATCH_SIZE}..."
                )

                sparql_results = self._execute_batched_sparql_query(batch_master_ids)
                if not sparql_results:
                    continue

                for result in sparql_results:
                    uniprot_ac = result.get("protein", {}).get("value", "").split("/")[-1]
                    if uniprot_ac:
                        all_uniprot_ids.add(uniprot_ac)

                time.sleep(1)

            logger.info(f"Found {len(all_uniprot_ids)} unique UniProt IDs from SPARQL queries")

            # Step 3: Use ID Mapper to get detailed information
            uniprot_details_df = self._get_uniprot_details_via_id_mapper(all_uniprot_ids, output_path)
            if uniprot_details_df is None:
                logger.error("Failed to get UniProt details via ID Mapper")
                return

            # Step 4: Process results and filter by specific Rhea IDs
            all_found_uniprot_acs = {}

            for _, row in uniprot_details_df.iterrows():
                uniprot_ac = row.get('Entry', '')
                gene_name = row.get('Gene Names (primary)', '')
                organism_id = row.get('Organism (ID)', '')
                rhea_ids = row.get('Rhea ID', '')
                catalytic_activity = row.get('Catalytic activity', '')

                if not uniprot_ac:
                    continue

                # Parse catalytic activity to find matching Rhea IDs with directions
                matching_activities = self._parse_catalytic_activity(catalytic_activity, target_rhea_ids)

                if matching_activities:
                    # Map back to our reaction IDs
                    for specific_rhea_id, direction in matching_activities:
                        # Find which of our reactions correspond to this specific Rhea ID
                        for reaction_id, master_id, target_specific_rhea_id in valid_reactions_to_query:
                            if target_specific_rhea_id == specific_rhea_id:
                                enzyme_info = {
                                    "uniprot_ac": uniprot_ac,
                                    "taxon_id": organism_id,
                                    "organism_name": None,  # Could be enhanced by organism lookup
                                    "gene_name": gene_name,
                                    "rhea_master_ids": rhea_ids,  # From "Rhea ID" column
                                    "catalytic_activities": [(specific_rhea_id, direction)],
                                    "physiological_direction": direction
                                }

                                if reaction_id not in discovered_data:
                                    discovered_data[reaction_id] = []
                                discovered_data[reaction_id].append(enzyme_info)

                                if uniprot_ac not in all_found_uniprot_acs:
                                    all_found_uniprot_acs[uniprot_ac] = gene_name

                                logger.debug(
                                    f"Associated enzyme {uniprot_ac} with reaction_id {reaction_id} via specific Rhea {specific_rhea_id} (direction: {direction})."
                                )

        # Step 5: Save outputs
        discovered_data_path = output_path / DISCOVERED_DATA_FILE
        discovered_data_str_keys = {str(k): v for k, v in discovered_data.items()}
        with open(discovered_data_path, "w", encoding="utf-8") as f:
            json.dump(discovered_data_str_keys, f, indent=2)
        logger.info(
            f"Saved {len(discovered_data)} reactions with potential new enzymes to: {discovered_data_path}"
        )

        if all_found_uniprot_acs:
            map_template_path = output_path / MAP_TEMPLATE_FILE
            map_df = pd.DataFrame(
                all_found_uniprot_acs.items(), columns=["uniprot_ac", "gene_name"]
            )
            map_df.sort_values("uniprot_ac", inplace=True)
            map_df.to_csv(map_template_path, sep="\t", index=False)
            logger.info(
                f"Generated gene mapping template with {len(map_df)} UniProt ACs at: {map_template_path}"
            )

        logger.info("--- Phase 1: Discovery Complete ---")

    def debug_rhea_id(self, rhea_directions_path: str, rhea_id: int, output_dir: str = None) -> None:
        """
        Debug phase: Traces the complete pipeline for a specific Rhea ID to identify issues.
        """
        logger.info(f"--- Debug Mode: Tracing Rhea ID {rhea_id} ---")

        try:
            directions_df = pd.read_csv(rhea_directions_path, sep="\t", dtype=str)
            for col in directions_df.columns:
                directions_df[col] = directions_df[col].str.strip()
            logger.info(f"Loaded {len(directions_df)} Rhea direction quartets.")
        except Exception as e:
            logger.critical(f"Failed to load Rhea directions file: {e}")
            return

        with self.db_updater as updater:
            # Step 1: Check if Rhea ID exists in local database
            logger.info(f"Step 1: Checking local database for Rhea ID {rhea_id}")
            local_reactions_df = updater._query_to_df(
                "SELECT reaction_id, reaction_text, rhea_id FROM lipograph.reactions WHERE rhea_id = %s;",
                (rhea_id,)
            )

            if local_reactions_df.empty:
                logger.error(f"❌ Rhea ID {rhea_id} not found in local database!")
                return

            reaction_row = local_reactions_df.iloc[0]
            reaction_id = reaction_row["reaction_id"]
            reaction_text = reaction_row["reaction_text"]

            logger.info(f"✅ Found in database:")
            logger.info(f"   Reaction ID: {reaction_id}")
            logger.info(f"   Reaction Text: {reaction_text}")

            # Step 2: Parse reaction direction
            logger.info(f"Step 2: Parsing reaction directionality")
            components, is_reversible = split_reaction_text(reaction_text)
            if not components:
                logger.error(f"❌ Could not parse reaction text: {reaction_text}")
                return

            if is_reversible:
                direction_inferred = "bi"
            elif ' <= ' in reaction_text:
                direction_inferred = "rl"  # right-to-left
            else:
                direction_inferred = "lr"  # left-to-right (=> or default)
            logger.info(
                f"✅ Parsed direction: {direction_inferred} ({'bidirectional' if is_reversible else 'left-to-right'})")

            # Step 3: Check directions file
            logger.info(f"Step 3: Looking up in directions file")
            rhea_id_str = str(rhea_id).strip()

            quartet_rows = directions_df[directions_df.isin([rhea_id_str]).any(axis=1)]
            if quartet_rows.empty:
                logger.error(f"❌ Rhea ID {rhea_id_str} not found in directions file!")
                logger.info("   This could be why it's not being processed.")
                return

            quartet_row = quartet_rows.iloc[0]
            logger.info(f"✅ Found quartet:")
            logger.info(f"   Master: {quartet_row['RHEA_ID_MASTER']}")
            logger.info(f"   LR: {quartet_row['RHEA_ID_LR']}")
            logger.info(f"   RL: {quartet_row['RHEA_ID_RL']}")
            logger.info(f"   BI: {quartet_row['RHEA_ID_BI']}")

            # Step 4: Check directionality validation
            logger.info(f"Step 4: Validating directionality")
            is_valid_direction = (
                                         direction_inferred == "lr" and quartet_row["RHEA_ID_LR"] == rhea_id_str
                                 ) or (
                                         direction_inferred == "rl" and quartet_row["RHEA_ID_RL"] == rhea_id_str
                                 # Add this line
                                 ) or (
                                         direction_inferred == "bi" and quartet_row["RHEA_ID_BI"] == rhea_id_str
                                 )

            if not is_valid_direction:
                logger.error(f"❌ Directionality validation failed!")
                logger.error(f"   Inferred direction: {direction_inferred}")
                logger.error(f"   Expected for LR: {quartet_row['RHEA_ID_LR']}")
                logger.error(f"   Expected for BI: {quartet_row['RHEA_ID_BI']}")
                logger.error(f"   Actual Rhea ID: {rhea_id_str}")
                logger.info("   This is likely why the Rhea ID is being excluded!")

                # Check if it might be RL (right-to-left)
                if quartet_row["RHEA_ID_RL"] == rhea_id_str:
                    logger.info(f"   💡 This appears to be a RIGHT-TO-LEFT reaction (RL)!")
                    logger.info(f"   Your reaction text suggests right-to-left: {reaction_text}")
                    logger.info(f"   But your parser inferred: {direction_inferred}")
                    logger.info(f"   Consider updating the reaction parsing logic to detect RL direction.")

                return

            logger.info(f"✅ Directionality validation passed")
            master_rhea_id = int(quartet_row["RHEA_ID_MASTER"])
            logger.info(f"   Will query with master ID: {master_rhea_id}")

            # Step 5: Test SPARQL query
            logger.info(f"Step 5: Testing SPARQL query with master ID {master_rhea_id}")
            sparql_results = self._execute_batched_sparql_query([master_rhea_id])

            if not sparql_results:
                logger.error(f"❌ SPARQL query returned no results for master ID {master_rhea_id}")
                logger.info("   This means no UniProt proteins are associated with this master Rhea ID.")
                return

            uniprot_ids = set()
            for result in sparql_results:
                uniprot_ac = result.get("protein", {}).get("value", "").split("/")[-1]
                if uniprot_ac:
                    uniprot_ids.add(uniprot_ac)

            logger.info(f"✅ SPARQL found {len(uniprot_ids)} UniProt IDs:")
            for uniprot_id in sorted(list(uniprot_ids)[:10]):  # Show first 10
                logger.info(f"   {uniprot_id}")
            if len(uniprot_ids) > 10:
                logger.info(f"   ... and {len(uniprot_ids) - 10} more")

            # Step 6: Test ID Mapper (optional, only if output_dir provided)
            if output_dir and uniprot_ids:
                logger.info(f"Step 6: Testing ID Mapper with {len(uniprot_ids)} UniProt IDs")
                output_path = Path(output_dir)
                output_path.mkdir(parents=True, exist_ok=True)

                # Test with subset to avoid long waits
                test_ids = set(list(uniprot_ids)[:20])  # Test with first 20 IDs
                uniprot_details_df = self._get_uniprot_details_via_id_mapper(test_ids, output_path)

                if uniprot_details_df is None:
                    logger.error(f"❌ ID Mapper failed")
                    return

                logger.info(f"✅ ID Mapper returned {len(uniprot_details_df)} entries")

                # Step 7: Test catalytic activity parsing
                logger.info(f"Step 7: Testing catalytic activity parsing for Rhea ID {rhea_id}")
                target_rhea_ids = {rhea_id}
                matching_count = 0

                for _, row in uniprot_details_df.iterrows():
                    uniprot_ac = row.get('Entry', '')
                    catalytic_activity = row.get('Catalytic activity', '')

                    if pd.isna(catalytic_activity) or not catalytic_activity:
                        continue

                    matching_activities = self._parse_catalytic_activity(catalytic_activity, target_rhea_ids)

                    if matching_activities:
                        matching_count += 1
                        logger.info(f"   ✅ {uniprot_ac}: Found {len(matching_activities)} matching activities")
                        for specific_rhea_id, direction in matching_activities:
                            logger.info(f"      Rhea {specific_rhea_id}: {direction}")

                        # Show a sample of the catalytic activity text
                        if catalytic_activity:
                            sample = catalytic_activity[:200] + "..." if len(
                                catalytic_activity) > 200 else catalytic_activity
                            logger.info(f"      Sample activity: {sample}")

                if matching_count == 0:
                    logger.error(f"❌ No proteins had catalytic activities matching Rhea ID {rhea_id}")
                    logger.info("   This could be because:")
                    logger.info(
                        "   1. The proteins don't have the specific directional Rhea ID in their catalytic activities")
                    logger.info("   2. The catalytic activity parsing is not finding the pattern")
                    logger.info("   3. UniProt data doesn't include this specific directional ID")
                else:
                    logger.info(f"✅ Found {matching_count} proteins with matching catalytic activities")

            # Step 8: Check discovered_enzymes.json if it exists
            if output_dir:
                discovered_file = Path(output_dir) / DISCOVERED_DATA_FILE
                if discovered_file.exists():
                    logger.info(f"Step 8: Checking discovered_enzymes.json")
                    try:
                        with open(discovered_file, 'r', encoding='utf-8') as f:
                            discovered_data = json.load(f)

                        reaction_id_str = str(reaction_id)
                        if reaction_id_str in discovered_data:
                            enzymes = discovered_data[reaction_id_str]
                            logger.info(
                                f"✅ Reaction {reaction_id} found in discovered_enzymes.json with {len(enzymes)} enzymes")
                            for i, enzyme in enumerate(enzymes[:5]):  # Show first 5
                                logger.info(
                                    f"   {i + 1}. {enzyme.get('uniprot_ac', 'N/A')} - {enzyme.get('gene_name', 'N/A')}")
                        else:
                            logger.error(f"❌ Reaction {reaction_id} NOT found in discovered_enzymes.json")
                            logger.info("   Available reaction IDs:")
                            for rid in sorted(discovered_data.keys())[:10]:
                                logger.info(f"   {rid}")
                    except Exception as e:
                        logger.error(f"❌ Could not read discovered_enzymes.json: {e}")

        logger.info(f"--- Debug Complete for Rhea ID {rhea_id} ---")

    def ingest_mapped_enzymes(self, mapping_file_path: str, discovered_data_path: str, commit: bool) -> None:
        """
        Phase 2: Ingests mapped enzymes and reuses existing reaction pair definitions
        to create new connections in the database.
        """
        logger.info("--- Phase 2: Ingesting Mapped Enzymes (Optimized) ---")
        if not commit:
            logger.warning("DRY RUN MODE: No changes will be committed to the database.")

        try:
            map_df = pd.read_csv(mapping_file_path, sep="\t", dtype=str)
            gene_map = map_df.dropna(subset=["gene_name"]).set_index("uniprot_ac")["gene_name"].to_dict()
            logger.info(f"Loaded {len(gene_map)} UniProt AC to gene name mappings.")

            with open(discovered_data_path, "r", encoding="utf-8") as f:
                discovered_data = json.load(f)
            logger.info(f"Loaded discovered data for {len(discovered_data)} reactions.")

        except FileNotFoundError as e:
            logger.critical(f"Required input file not found: {e}. Aborting.");
            return
        except Exception as e:
            logger.critical(f"Failed to load input files: {e}");
            return

        with self.db_updater as updater:
            try:
                with updater.conn.transaction():
                    logger.info("Database transaction started.")

                    # Step 1: Pre-fetch all existing reaction pairs for the reactions we're about to modify.
                    # This is much more efficient than querying inside the loop.
                    all_reaction_ids = [int(rid) for rid in discovered_data.keys()]
                    if not all_reaction_ids:
                        logger.info("No discovered data to ingest.");
                        return

                    logger.info(f"Pre-fetching existing reaction pairs for {len(all_reaction_ids)} reactions...")
                    existing_pairs_df = updater._query_to_df(
                        """
                        SELECT re.reaction_id, rp.reactant_molecule_id, rp.product_molecule_id
                        FROM lipograph.reaction_pairs rp
                        JOIN lipograph.reaction_enzyme re ON rp.reaction_enzyme_id = re.reaction_enzyme_id
                        WHERE re.reaction_id = ANY(%s)
                        GROUP BY re.reaction_id, rp.reactant_molecule_id, rp.product_molecule_id;
                        """,
                        (all_reaction_ids,)
                    )

                    # Group the pairs by reaction_id for easy lookup
                    pairs_by_reaction = existing_pairs_df.groupby('reaction_id')[
                        ['reactant_molecule_id', 'product_molecule_id']].apply(
                        lambda x: [tuple(row) for row in x.to_numpy()]
                    ).to_dict()
                    logger.info(f"Found existing pairs for {len(pairs_by_reaction)} reactions.")

                    # Step 2: Iterate through discovered data and create new links
                    for reaction_id_str, enzymes in discovered_data.items():
                        reaction_id = int(reaction_id_str)

                        # Get the predefined molecule pairs for this reaction
                        molecule_pairs = pairs_by_reaction.get(reaction_id)
                        if not molecule_pairs:
                            logger.warning(
                                f"Reaction ID {reaction_id} has no existing pairs in the database. Cannot create new enzyme links for it. Please verify the initial import.")
                            continue

                        for enzyme_info in enzymes:
                            uniprot_ac = enzyme_info["uniprot_ac"]
                            if uniprot_ac not in gene_map:
                                continue

                            gene_name = gene_map[uniprot_ac]
                            logger.info(
                                f"Processing mapped enzyme {uniprot_ac} ('{gene_name}') for reaction_id {reaction_id}")

                            # A. Ensure organism and enzyme exist
                            organism_id = updater._ensure_organism(enzyme_info["taxon_id"],
                                                                   enzyme_info["organism_name"])
                            if not organism_id: continue

                            enzyme_id = updater._ensure_enzyme(uniprot_ac, gene_name, organism_id)
                            if not enzyme_id: continue

                            # B. Ensure the link between reaction and enzyme exists
                            re_id, was_inserted = updater._insert_entity(
                                "lipograph.reaction_enzyme", ["reaction_id", "enzyme_id"],
                                {"reaction_id": reaction_id, "enzyme_id": enzyme_id},
                                return_id_col="reaction_enzyme_id"
                            )

                            # C. If the reaction_enzyme link is NEW, create all its pairs.
                            #    If it already existed, we assume its pairs are also correct.
                            if was_inserted and re_id:
                                logger.debug(
                                    f"New reaction_enzyme link created (ID: {re_id}). Creating its reaction pairs...")
                                for reactant_id, product_id in molecule_pairs:
                                    updater._ensure_reaction_pairs(re_id, [reactant_id], [product_id],
                                                                   is_sl_reversible=False)
                                logger.info(
                                    f"Created {len(molecule_pairs)} reaction pairs for new link (re_id: {re_id}).")
                            elif not was_inserted:
                                logger.debug(
                                    f"Reaction_enzyme link for reaction {reaction_id} and enzyme {enzyme_id} already exists. Skipping pair creation.")

                    if not commit:
                        raise RuntimeError("DRY RUN: Rolling back transaction.")

                logger.info("Transaction committed successfully.")

            except Exception as e:
                logger.error(f"An error occurred during ingestion. Transaction will be rolled back. Error: {e}",
                             exc_info=True)

        logger.info("--- Phase 2: Ingestion Complete ---")

    def preview_ingest(self, mapping_file_path: str, discovered_data_path: str) -> None:
        """
        Phase 2 Preview: Provides detailed statistics on new and existing connections
        that would be made by the 'ingest' phase, without changing the database.
        """
        logger.info("--- Phase 2 Preview: Analyzing for New Connections ---")

        try:
            map_df = pd.read_csv(mapping_file_path, sep="\t", dtype=str)
            gene_map = (
                map_df.dropna(subset=["gene_name"])
                .set_index("uniprot_ac")["gene_name"]
                .to_dict()
            )
            logger.info(f"Loaded {len(gene_map)} UniProt AC to gene name mappings.")

            with open(discovered_data_path, "r", encoding="utf-8") as f:
                discovered_data = json.load(f)
            logger.info(f"Loaded discovered data for {len(discovered_data)} reactions.")

        except FileNotFoundError as e:
            logger.critical(f"Required input file not found: {e}. Aborting.")
            return
        except Exception as e:
            logger.critical(f"Failed to load input files: {e}")
            return

        # Data collection lists for detailed stats
        new_links_details = []

        with self.db_updater as updater:
            logger.info("Fetching existing data from the database for comparison...")
            # Fetch existing reaction-enzyme connections
            existing_links_df = updater._query_to_df(
                """
                SELECT re.reaction_id, e.uniprot_id
                FROM lipograph.reaction_enzyme re
                JOIN lipograph.enzymes e ON re.enzyme_id = e.enzyme_id;
                """
            )
            existing_connections_set = set(
                tuple(row) for row in existing_links_df.itertuples(index=False)
            )

            # Fetch existing enzymes
            existing_enzymes_df = updater._query_to_df("SELECT uniprot_id FROM lipograph.enzymes;")
            existing_uniprot_ids = set(existing_enzymes_df['uniprot_id'])

            # Fetch existing organism IDs
            existing_organisms_df = updater._query_to_df("SELECT organism_id FROM lipograph.organism;")
            existing_organism_ids = set(existing_organisms_df['organism_id'])

            # Fetch reaction texts for reporting
            all_reaction_ids_in_file = [int(rid) for rid in discovered_data.keys()]
            reaction_text_df = updater._query_to_df(
                "SELECT reaction_id, reaction_text FROM lipograph.reactions WHERE reaction_id = ANY(%s)",
                (all_reaction_ids_in_file,)
            )
            reaction_id_to_text = reaction_text_df.set_index('reaction_id')['reaction_text'].to_dict()

            logger.info("Comparing discovered data against current database state...")
            # Simulate the ingestion loop
            for reaction_id_str, enzymes in discovered_data.items():
                reaction_id = int(reaction_id_str)

                for enzyme_info in enzymes:
                    uniprot_ac = enzyme_info["uniprot_ac"]
                    if uniprot_ac not in gene_map:
                        continue

                    if (reaction_id, uniprot_ac) not in existing_connections_set:
                        new_links_details.append({
                            "reaction_id": reaction_id,
                            "reaction_text": reaction_id_to_text.get(reaction_id, "N/A"),
                            "uniprot_id": uniprot_ac,
                            "is_new_enzyme": uniprot_ac not in existing_uniprot_ids,
                            "organism_id": int(enzyme_info["taxon_id"]) if enzyme_info.get("taxon_id") else None
                        })

        # --- Generate Detailed Statistics ---
        if not new_links_details:
            logger.info("--- Ingestion Preview Summary ---")
            logger.info("✅ No new connections would be added. The database is up-to-date with these files.")
            return

        stats_df = pd.DataFrame(new_links_details)

        # New entity counts
        total_new_links = len(stats_df)
        newly_affected_reactions = stats_df['reaction_id'].nunique()
        new_uniprot_ids_to_add = set(stats_df[stats_df['is_new_enzyme']]['uniprot_id'])
        total_new_enzymes = len(new_uniprot_ids_to_add)
        all_organism_ids_in_new_links = set(stats_df.dropna(subset=['organism_id'])['organism_id'].astype(int))
        new_organism_ids_to_add = all_organism_ids_in_new_links - existing_organism_ids
        total_new_organisms = len(new_organism_ids_to_add)

        # *** NEW: Calculate how many existing enzymes get new connections ***
        # Filter for new links that involve enzymes already in the DB
        existing_enzymes_getting_new_links_df = stats_df[stats_df['is_new_enzyme'] == False]
        # Count how many unique enzymes this corresponds to
        count_existing_enzymes_with_new_links = existing_enzymes_getting_new_links_df['uniprot_id'].nunique()

        # Existing entity counts
        existing_enzyme_count = len(existing_uniprot_ids)
        existing_organism_count = len(existing_organism_ids)
        existing_connection_count = len(existing_connections_set)

        # Top 10 reactions
        top_10_reactions = (
            stats_df.groupby(['reaction_id', 'reaction_text'])
            .size()
            .reset_index(name='new_enzyme_count')
            .sort_values(by='new_enzyme_count', ascending=False)
            .head(10)
        )

        # --- Print Final Report ---
        print("\n" + "=" * 65)
        print("              Ingestion Preview: Statistical Report")
        print("=" * 65)

        print("\n--- Database State Summary ---")
        print(f"                                   {'Existing':<12} | {'To Be Added':<12} | {'New Total':<12}")
        print(f"----------------------------------- {'-' * 12} + {'-' * 12} + {'-' * 12}")
        print(
            f"Organisms                          {existing_organism_count:<12} | {total_new_organisms:<12} | {existing_organism_count + total_new_organisms:<12}")
        print(
            f"Enzymes (Unique UniProt IDs)       {existing_enzyme_count:<12} | {total_new_enzymes:<12} | {existing_enzyme_count + total_new_enzymes:<12}")
        print(
            f"Reaction-Enzyme Connections        {existing_connection_count:<12} | {total_new_links:<12} | {existing_connection_count + total_new_links:<12}")

        # *** NEW: Print the new detailed breakdown ***
        print("\n--- Breakdown of New Connections ---")
        print(f"- The {total_new_links} new connections will affect {newly_affected_reactions} unique reactions.")
        print(
            f"- Of these, {count_existing_enzymes_with_new_links} EXISTING enzymes will receive new reaction connections.")

        print("\n--- Top 10 Reactions Receiving the Most New Enzyme Connections ---")
        print(top_10_reactions.to_string(index=False))

        print("\n" + "=" * 65)
        logger.info("This was a read-only check. No changes were made to the database.")
        logger.info("Run the 'ingest --commit' command to apply these changes.")


def main():
    """Main function to run the CLI."""
    parser = argparse.ArgumentParser(
        description="Enrich LipidCRED database with enzymes from Rhea/UniProt.",
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
        "--phase",
        choices=["discover", "ingest", "debug", "preview-ingest"],
        required=True,
        help="The operational phase to run.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=f"db_updates/enrichment_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        help="Directory for output files from 'discover' phase and input for 'ingest' phase.",
    )
    parser.add_argument(
        "--rhea-directions",
        type=str,
        default="data/rhea-directions.tsv",
        help="[discover/debug phase] Path to the rhea-directions.tsv file.",  # Updated help
    )
    parser.add_argument(
        "--map-file",
        type=str,
        help="[ingest phase] Path to the user-completed uniprot_to_gene_map.tsv file.",
    )
    parser.add_argument(
        "--commit",
        action="store_true",
        help="[ingest phase] Commit changes to the database. Default is a dry run.",
    )
    parser.add_argument(
        "--rhea-id",
        type=int,
        help="[debug phase] Specific Rhea ID to debug.",
    )  # Added new argument

    args = parser.parse_args()

    log_dir = Path(args.output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file_path = log_dir / f"enricher_run_{args.phase}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    logging.basicConfig(
        level=args.log_level.upper(),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        # This now configures handlers for both file and console
        handlers=[
            logging.FileHandler(log_file_path, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )

    try:
        config = ConfigLoader(args.config_file)
        enricher = RheaEnzymeEnricher(config)

        if args.phase == "discover":
            if not args.rhea_directions:
                parser.error("--rhea-directions is required for the 'discover' phase.")
            enricher.discover_enzymes(args.rhea_directions, args.output_dir)

        elif args.phase == "ingest":
            default_map_file = Path(args.output_dir) / MAP_TEMPLATE_FILE.replace(
                "_template", ""
            )
            map_file = args.map_file or str(default_map_file)
            discovered_file = Path(args.output_dir) / DISCOVERED_DATA_FILE

            if not Path(map_file).exists() or not discovered_file.exists():
                parser.error(
                    f"--map-file ('{map_file}') and the discovered data file ('{discovered_file}') "
                    f"must exist for the 'ingest' phase. Ensure you've run 'discover' first."
                )

            enricher.ingest_mapped_enzymes(map_file, str(discovered_file), args.commit)

        elif args.phase == "preview-ingest":
            default_map_file = Path(args.output_dir) / MAP_TEMPLATE_FILE.replace(
                "_template", ""
            )
            map_file = args.map_file or str(default_map_file)
            discovered_file = Path(args.output_dir) / DISCOVERED_DATA_FILE

            if not Path(map_file).exists() or not discovered_file.exists():
                parser.error(
                    f"--map-file ('{map_file}') and the discovered data file ('{discovered_file}') "
                    f"must exist for the 'preview-ingest' phase."
                )

            enricher.preview_ingest(map_file, str(discovered_file))

        elif args.phase == "debug":  # Added debug phase handler
            if not args.rhea_id:
                parser.error("--rhea-id is required for the 'debug' phase.")
            if not args.rhea_directions:
                parser.error("--rhea-directions is required for the 'debug' phase.")

            enricher.debug_rhea_id(args.rhea_directions, args.rhea_id, args.output_dir)

    except FileNotFoundError as e:
        logger.critical(f"A required file was not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.critical(f"An unexpected error occurred: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
