import sys
import collections
import pickle
import json
from pathlib import Path
import pandas as pd
import psycopg
from psycopg.rows import dict_row
from typing import Set, Tuple, Optional, Dict


# --- Path Setup ---
def setup_project_path():
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir
    src_path = project_root / "src"
    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))
    return project_root


PROJECT_ROOT = setup_project_path()
from src.data.config_loader import ConfigLoader

# Import the LipidParser directly
from src.parsing.lipid_parser import LipidParser, LipidComponents, LipidParserFactory, LipidComponentCache

from src.utils.db_updater import LipidDatabaseUpdater
from src.utils.text_utils import split_reaction_text, generate_reaction_combinations


class GeneralizationLinker:
    """
    A cached script to link generalized reactions to their specific counterparts
    by comparing the classes of their lipid participants.
    """

    def __init__(self, config_loader: ConfigLoader):
        self.config_loader = config_loader
        self.conn = None
        self.cursor = None
        self.lipid_parser: Optional[LipidParser] = None
        self.db_updater: Optional[LipidDatabaseUpdater] = None

        # Cache file paths
        self.cache_dir = PROJECT_ROOT / "cache"
        self.cache_dir.mkdir(exist_ok=True)
        self.reaction_cache_file = self.cache_dir / "reaction_classes_cache.pickle"
        self.reaction_list_cache_file = self.cache_dir / "reaction_lists_cache.json"

    def connect_and_setup(self):
        """Establishes database connection and initializes the LipidParser."""
        try:
            # Create connection
            self.conn = psycopg.connect(
                dbname=self.config_loader.get_db_name(),
                user=self.config_loader.get_db_user(),
                password=self.config_loader.get_db_password(),
                host=self.config_loader.get_db_host(),
                port=self.config_loader.get_db_port(),
                options=f"-c search_path={self.config_loader.get_db_schema()}",
            )
            self.conn.row_factory = dict_row
            self.cursor = self.conn.cursor()

            # Handle any existing transaction state
            self.handle_transaction_state()

            # Initialize LipidParser
            lipid_component_cache = LipidComponentCache()
            lipid_parser_factory = LipidParserFactory(self.config_loader)
            self.lipid_parser = LipidParser(
                factory=lipid_parser_factory,
                components_cache=lipid_component_cache
            )
            self.db_updater = LipidDatabaseUpdater(self.config_loader)
            self.db_updater.conn = self.conn
            self.db_updater.cursor = self.cursor

            print("Successfully connected to the database and initialized LipidParser.")
        except Exception as e:
            print(f"Error during setup: {e}", file=sys.stderr)
            sys.exit(1)

    def handle_transaction_state(self):
        """Handle any existing transaction state on the connection."""
        try:
            print(f"Connection transaction status: {self.conn.info.transaction_status}")

            # If we're in a transaction, commit or rollback to clean state
            if self.conn.info.transaction_status == 2:  # INTRANS
                print("Connection is in transaction state. Committing to clean state...")
                self.conn.commit()
                print("Transaction committed. Connection is now clean.")
            elif self.conn.info.transaction_status == 3:  # INERROR
                print("Connection is in error state. Rolling back...")
                self.conn.rollback()
                print("Transaction rolled back. Connection is now clean.")

        except Exception as e:
            print(f"Warning: Could not clean transaction state: {e}")
            # If we can't clean it, try rollback as last resort
            try:
                self.conn.rollback()
                print("Rollback successful.")
            except:
                print("Could not rollback. Connection may need to be recreated.")

    def close(self):
        """Closes the database connection."""
        if self.cursor:
            self.cursor.close()
        if self.conn:
            self.conn.close()
        print("Database connection closed.")

    def load_reaction_lists_cache(self) -> Optional[Dict]:
        """Load cached reaction IDs if available and recent."""
        if not self.reaction_list_cache_file.exists():
            return None

        try:
            with open(self.reaction_list_cache_file, 'r') as f:
                cache_data = json.load(f)

            print(f"Loaded reaction lists from cache: {len(cache_data.get('generalized', []))} generalized, {len(cache_data.get('specific', []))} specific")
            return cache_data
        except Exception as e:
            print(f"Could not load reaction lists cache: {e}")
            return None

    def save_reaction_lists_cache(self, generalized_ids: list, specific_ids: list):
        """Save reaction ID lists to cache."""
        try:
            cache_data = {
                'generalized': generalized_ids,
                'specific': specific_ids,
                'timestamp': pd.Timestamp.now().isoformat()
            }
            with open(self.reaction_list_cache_file, 'w') as f:
                json.dump(cache_data, f)
            print(f"Saved reaction lists to cache: {len(generalized_ids)} generalized, {len(specific_ids)} specific")
        except Exception as e:
            print(f"Could not save reaction lists cache: {e}")

    def load_reaction_classes_cache(self) -> Dict:
        """Load cached reaction class data if available."""
        if not self.reaction_cache_file.exists():
            return {}

        try:
            with open(self.reaction_cache_file, 'rb') as f:
                cache_data = pickle.load(f)
            print(f"Loaded {len(cache_data)} cached reaction class analyses")
            return cache_data
        except Exception as e:
            print(f"Could not load cache: {e}")
            return {}

    def save_reaction_classes_cache(self, cache_data: Dict):
        """Save reaction class data to cache."""
        try:
            with open(self.reaction_cache_file, 'wb') as f:
                pickle.dump(cache_data, f)
            print(f"Saved {len(cache_data)} reaction analyses to cache")
        except Exception as e:
            print(f"Could not save cache: {e}")

    def get_reaction_participants_classes(
            self, reaction_id: int
    ) -> Optional[Tuple[Set[str], Set[str]]]:
        """For a reaction_id, fetches its molecule participants and returns their headgroup classes."""
        try:
            # Fetch reactants' cleaned abbreviations
            self.cursor.execute(
                """
                SELECT DISTINCT m.cleaned_abbreviation FROM lipograph.reaction_pairs rp
                JOIN lipograph.molecules m ON rp.reactant_molecule_id = m.molecule_id
                WHERE rp.reaction_enzyme_id IN (SELECT reaction_enzyme_id FROM lipograph.reaction_enzyme WHERE reaction_id = %s);
            """,
                (reaction_id,),
            )
            reactant_abbrs = {
                row["cleaned_abbreviation"]
                for row in self.cursor.fetchall()
                if row["cleaned_abbreviation"]
            }

            # Fetch products' cleaned abbreviations
            self.cursor.execute(
                """
                SELECT DISTINCT m.cleaned_abbreviation FROM lipograph.reaction_pairs rp
                JOIN lipograph.molecules m ON rp.product_molecule_id = m.molecule_id
                WHERE rp.reaction_enzyme_id IN (SELECT reaction_enzyme_id FROM lipograph.reaction_enzyme WHERE reaction_id = %s);
            """,
                (reaction_id,),
            )
            product_abbrs = {
                row["cleaned_abbreviation"]
                for row in self.cursor.fetchall()
                if row["cleaned_abbreviation"]
            }

            if not reactant_abbrs and not product_abbrs:
                return None

            reactant_classes = set()
            for abbr in reactant_abbrs:
                try:
                    if abbr.startswith('FA('):
                        continue
                    else:
                        parsed = self.lipid_parser.parse_lipid(abbr)
                        reactant_classes.add(parsed.headgroup)
                except Exception:
                    # Skip other non-lipid molecules
                    pass

            product_classes = set()
            for abbr in product_abbrs:
                try:
                    if abbr.startswith('FA('):
                        continue
                    else:
                        parsed = self.lipid_parser.parse_lipid(abbr)
                        product_classes.add(parsed.headgroup)
                except Exception:
                    # Skip other non-lipid molecules
                    pass

            # Only return if we found at least some lipids
            if reactant_classes or product_classes:
                return reactant_classes, product_classes
            else:
                return None

        except Exception as e:
            print(
                f"  - Warning: Could not parse participants for reaction_id {reaction_id}. Error: {e}",
                file=sys.stderr,
            )
            return None

    def get_reaction_texts(self, reaction_ids: list) -> Dict[int, str]:
        """Fetches reaction texts for a list of IDs."""
        if not reaction_ids:
            return {}
        self.cursor.execute(
            "SELECT reaction_id, reaction_text FROM lipograph.reactions WHERE reaction_id = ANY(%s)",
            (reaction_ids,),
        )
        return {
            row["reaction_id"]: row["reaction_text"] for row in self.cursor.fetchall()
        }

    def fetch_reaction_ids(self) -> Tuple[list, list]:
        """Fetch generalized and specific reaction IDs, using cache if available."""
        # Try to load from cache first
        cached_lists = self.load_reaction_lists_cache()
        if cached_lists:
            use_cache = input("Use cached reaction lists? (yes/no): ").lower()
            if use_cache == 'yes':
                return cached_lists['generalized'], cached_lists['specific']

        print("Fetching reaction IDs from database...")

        print("Fetching candidate generalized reactions (rhea_id IS NULL)...")
        self.cursor.execute(
            "SELECT reaction_id FROM lipograph.reactions WHERE rhea_id IS NULL;"
        )
        generalized_reactions_ids = [
            row["reaction_id"] for row in self.cursor.fetchall()
        ]

        print("Fetching all specific reactions (rhea_id IS NOT NULL)...")
        self.cursor.execute(
            "SELECT reaction_id FROM lipograph.reactions WHERE rhea_id IS NOT NULL;"
        )
        specific_reactions_ids = [row["reaction_id"] for row in self.cursor.fetchall()]

        # Save to cache
        self.save_reaction_lists_cache(generalized_reactions_ids, specific_reactions_ids)

        return generalized_reactions_ids, specific_reactions_ids

    def build_reaction_classes_cache(self, all_reaction_ids: list) -> Dict:
        """Build cache of reaction classes, using existing cache where possible."""
        # Load existing cache
        reaction_class_cache = self.load_reaction_classes_cache()

        # Find which reactions we still need to process
        missing_ids = [rid for rid in all_reaction_ids if rid not in reaction_class_cache]

        if missing_ids:
            print(f"Need to analyze {len(missing_ids)} new reactions (have {len(reaction_class_cache)} cached)")

            for i, reaction_id in enumerate(missing_ids):
                if i > 0 and i % 100 == 0:
                    print(f"  Processed {i}/{len(missing_ids)} new reactions...")
                    # Save cache periodically
                    self.save_reaction_classes_cache(reaction_class_cache)

                participants = self.get_reaction_participants_classes(reaction_id)
                if participants:
                    reaction_class_cache[reaction_id] = participants

            # Save final cache
            self.save_reaction_classes_cache(reaction_class_cache)
        else:
            print("All reactions already cached!")

        return reaction_class_cache

    def debug_specific_reaction(self, reaction_id: int):
        """Debug why a specific reaction isn't being matched."""
        print(f"\n=== DEBUGGING REACTION {reaction_id} ===")

        # Check if it exists and get basic info
        self.cursor.execute(
            "SELECT reaction_id, reaction_text, rhea_id FROM lipograph.reactions WHERE reaction_id = %s",
            (reaction_id,)
        )
        reaction_info = self.cursor.fetchone()
        if not reaction_info:
            print(f"❌ Reaction {reaction_id} not found in database!")
            return

        print(f"✅ Reaction exists: {reaction_info['reaction_text']}")
        print(f"   RHEA ID: {reaction_info['rhea_id']}")

        # Check if it would be included in our specific reactions list
        is_specific = reaction_info['rhea_id'] is not None
        print(f"   Is specific reaction (has RHEA ID): {is_specific}")

        # Get participant classes for this reaction
        participants = self.get_reaction_participants_classes(reaction_id)
        if participants:
            reactant_classes, product_classes = participants
            print(f"   Reactant classes: {reactant_classes}")
            print(f"   Product classes: {product_classes}")

            # Find generalized reactions with same participant classes
            print(f"\n🔍 Looking for generalized reactions with matching classes...")
            self.cursor.execute(
                "SELECT reaction_id FROM lipograph.reactions WHERE rhea_id IS NULL AND (generalized_from_reaction_ids IS NULL OR generalized_from_reaction_ids = '')"
            )
            generalized_ids = [row["reaction_id"] for row in self.cursor.fetchall()]

            matches = []
            for gen_id in generalized_ids:
                gen_participants = self.get_reaction_participants_classes(gen_id)
                if gen_participants:
                    gen_reactant_classes, gen_product_classes = gen_participants
                    if (reactant_classes == gen_reactant_classes and
                            product_classes == gen_product_classes):
                        matches.append(gen_id)

            if matches:
                print(f"   ✅ Found {len(matches)} matching generalized reactions: {matches}")

                # Check if any of these already have this reaction linked
                for gen_id in matches:
                    self.cursor.execute(
                        "SELECT reaction_id, generalized_from_reaction_ids FROM lipograph.reactions WHERE reaction_id = %s",
                        (gen_id,)
                    )
                    gen_info = self.cursor.fetchone()
                    current_links = gen_info['generalized_from_reaction_ids'] or ""
                    if f"|{reaction_id}|" in current_links:
                        print(f"   ✅ Reaction {reaction_id} is already linked to generalized reaction {gen_id}")
                    else:
                        print(f"   ❌ Reaction {reaction_id} is NOT linked to generalized reaction {gen_id}")
                        print(f"      Current links for {gen_id}: {current_links}")
            else:
                print(f"   ❌ No matching generalized reactions found!")
        else:
            print(f"   ❌ Could not parse participants for reaction {reaction_id}")

    def check_existing_links(self, specific_reaction_id: int):
        """Check if a specific reaction is already linked to any generalized reaction."""
        print(f"\n=== CHECKING EXISTING LINKS FOR REACTION {specific_reaction_id} ===")

        self.cursor.execute(
            """
            SELECT reaction_id, reaction_text, generalized_from_reaction_ids 
            FROM lipograph.reactions 
            WHERE generalized_from_reaction_ids LIKE %s
            """,
            (f"%|{specific_reaction_id}|%",)
        )

        linked_reactions = self.cursor.fetchall()
        if linked_reactions:
            print(f"✅ Reaction {specific_reaction_id} is linked to {len(linked_reactions)} generalized reactions:")
            for row in linked_reactions:
                print(f"   - Gen ID {row['reaction_id']}: {row['reaction_text']}")
                print(f"     All links: {row['generalized_from_reaction_ids']}")
        else:
            print(f"❌ Reaction {specific_reaction_id} is not linked to any generalized reactions")

    def compare_reactions(self, reaction_id1: int, reaction_id2: int):
        """Compare two reactions to see why they might not be matching."""
        print(f"\n=== COMPARING REACTIONS {reaction_id1} vs {reaction_id2} ===")

        for rid in [reaction_id1, reaction_id2]:
            self.cursor.execute(
                "SELECT reaction_id, reaction_text, rhea_id FROM lipograph.reactions WHERE reaction_id = %s",
                (rid,)
            )
            info = self.cursor.fetchone()
            if info:
                print(f"Reaction {rid}: {info['reaction_text']}")
                print(f"  RHEA ID: {info['rhea_id']}")

                participants = self.get_reaction_participants_classes(rid)
                if participants:
                    reactant_classes, product_classes = participants
                    print(f"  Reactant classes: {reactant_classes}")
                    print(f"  Product classes: {product_classes}")
                else:
                    print(f"  ❌ Could not parse participants")
            else:
                print(f"❌ Reaction {rid} not found!")

        # Compare if both exist
        p1 = self.get_reaction_participants_classes(reaction_id1)
        p2 = self.get_reaction_participants_classes(reaction_id2)

        if p1 and p2:
            r1_reactants, r1_products = p1
            r2_reactants, r2_products = p2

            print(f"\n🔍 COMPARISON RESULTS:")
            print(f"  Reactants match: {r1_reactants == r2_reactants}")
            print(f"  Products match: {r1_products == r2_products}")
            print(f"  Overall match: {r1_reactants == r2_reactants and r1_products == r2_products}")

            if r1_reactants != r2_reactants:
                print(f"  Reactant differences:")
                print(f"    Only in {reaction_id1}: {r1_reactants - r2_reactants}")
                print(f"    Only in {reaction_id2}: {r2_reactants - r1_reactants}")

            if r1_products != r2_products:
                print(f"  Product differences:")
                print(f"    Only in {reaction_id1}: {r1_products - r2_products}")
                print(f"    Only in {reaction_id2}: {r2_products - r1_products}")

    def perform_updates_with_fresh_connection(self, updates_to_perform):
        """Perform updates with a completely fresh, dedicated connection."""
        print("\nCreating fresh connection specifically for updates...")

        # Create a brand new connection just for updates
        update_conn = psycopg.connect(
            dbname=self.config_loader.get_db_name(),
            user=self.config_loader.get_db_user(),
            password=self.config_loader.get_db_password(),
            host=self.config_loader.get_db_host(),
            port=self.config_loader.get_db_port(),
            options=f"-c search_path={self.config_loader.get_db_schema()}"
        )
        update_conn.row_factory = dict_row

        try:
            with update_conn.cursor() as cursor:
                update_query = "UPDATE lipograph.reactions SET generalized_from_reaction_ids = %s WHERE reaction_id = %s;"

                # Use explicit transaction block
                with update_conn.transaction():
                    for gen_id, spec_ids in updates_to_perform.items():
                        # Get existing links
                        cursor.execute(
                            "SELECT generalized_from_reaction_ids FROM lipograph.reactions WHERE reaction_id = %s",
                            (gen_id,)
                        )
                        existing = cursor.fetchone()['generalized_from_reaction_ids'] or ""

                        # Parse existing IDs
                        existing_ids = set()
                        if existing:
                            existing_ids = {int(x) for x in existing.strip('|').split('|') if x}

                        # Add new IDs
                        all_ids = existing_ids.union(set(spec_ids))
                        id_string = f"|{'|'.join(map(str, sorted(all_ids)))}|"

                        cursor.execute(update_query, (id_string, gen_id))
                        print(f"  Updated reaction_id {gen_id}: added {spec_ids} to existing {existing_ids}")

                print("Transaction block completed - changes should be committed.")

                # Verify immediately after transaction
                cursor.execute(
                    "SELECT COUNT(*) as count FROM lipograph.reactions WHERE generalized_from_reaction_ids IS NOT NULL AND generalized_from_reaction_ids != ''"
                )
                count_result = cursor.fetchone()
                print(f"Post-commit verification: {count_result['count']} updated records")

                # Show some actual updated records
                cursor.execute(
                    "SELECT reaction_id, generalized_from_reaction_ids FROM lipograph.reactions WHERE generalized_from_reaction_ids IS NOT NULL AND generalized_from_reaction_ids != '' LIMIT 5"
                )
                sample_results = cursor.fetchall()
                print("Sample updated records:")
                for row in sample_results:
                    print(f"  - ID {row['reaction_id']}: {row['generalized_from_reaction_ids']}")

            print(f"\nSUCCESS: Updated {len(updates_to_perform)} reactions.")
            print("You can now check in PyCharm SQL console - changes should be visible.")

        except Exception as e:
            print(f"Fresh connection update failed: {e}")
            raise
        finally:
            update_conn.close()
            print("Fresh update connection closed.")

    def propagate_enzymes_and_pairs(self, gen_reaction_id: int):
        """
        For a single generalized reaction, finds all enzymes from its specific
        children and ensures they are linked to the parent. It creates the
        necessary reaction_pairs for any new links by first checking for
        existing pairs, and falling back to parsing the reaction text if none exist.
        """
        print(f"\n  Propagating enzymes for Generalized Reaction ID: {gen_reaction_id}")

        # --- Part 1: Determine which new enzyme links to create ---
        dict_cursor = self.cursor  # Use the main dict cursor for all read operations
        dict_cursor.execute(
            "SELECT generalized_from_reaction_ids FROM lipograph.reactions WHERE reaction_id = %s",
            (gen_reaction_id,)
        )
        specific_ids_str = (dict_cursor.fetchone() or {}).get('generalized_from_reaction_ids')
        if not specific_ids_str:
            print(f"    - No specific children listed for reaction {gen_reaction_id}. Nothing to propagate.")
            return

        specific_child_ids = [int(x) for x in specific_ids_str.strip('|').split('|') if x]
        if not specific_child_ids:
            return

        dict_cursor.execute(
            "SELECT DISTINCT enzyme_id FROM lipograph.reaction_enzyme WHERE reaction_id = ANY(%s);",
            (specific_child_ids,)
        )
        child_enzyme_ids = {row['enzyme_id'] for row in dict_cursor.fetchall()}

        if not child_enzyme_ids:
            print("    - No enzymes found from children. Nothing to propagate.")
            return

        dict_cursor.execute(
            "SELECT enzyme_id FROM lipograph.reaction_enzyme WHERE reaction_id = %s;",
            (gen_reaction_id,)
        )
        parent_enzyme_ids = {row['enzyme_id'] for row in dict_cursor.fetchall()}

        new_enzyme_ids_to_link = child_enzyme_ids - parent_enzyme_ids

        if not new_enzyme_ids_to_link:
            print("    - Parent is already up-to-date with all child enzymes.")
            return

        print(f"    - Found {len(new_enzyme_ids_to_link)} new enzyme connections to propagate.")

        # --- Part 2: Define the molecular participants for the generalized reaction ---
        molecule_pairs = []
        is_reversible = False  # Default unless text parsing says otherwise

        # A. First, try to get pairs from existing links for this generalized reaction
        dict_cursor.execute(
            """
            SELECT rp.reactant_molecule_id, rp.product_molecule_id
            FROM lipograph.reaction_pairs rp
            JOIN lipograph.reaction_enzyme re ON rp.reaction_enzyme_id = re.reaction_enzyme_id
            WHERE re.reaction_id = %s
            GROUP BY rp.reactant_molecule_id, rp.product_molecule_id
            """,
            (gen_reaction_id,)
        )
        existing_pairs = dict_cursor.fetchall()

        if existing_pairs:
            molecule_pairs = existing_pairs
        else:
            # B. If no pairs exist (it's a new generalized reaction), parse its text
            print(
                f"    - No existing pairs found for reaction {gen_reaction_id}. Parsing reaction text to define participants.")
            dict_cursor.execute("SELECT reaction_text FROM lipograph.reactions WHERE reaction_id = %s",
                                (gen_reaction_id,))
            reaction_text_row = dict_cursor.fetchone()

            if not reaction_text_row or not reaction_text_row['reaction_text']:
                print(f"    - CRITICAL ERROR: Reaction {gen_reaction_id} has no text and no pairs. Cannot proceed.")
                return

            reaction_text = reaction_text_row['reaction_text']

            # Use the imported utility function
            components, is_reversible_from_text = split_reaction_text(reaction_text)
            if not components:
                print(f"    - CRITICAL ERROR: Could not parse reaction text '{reaction_text}'. Cannot create pairs.")
                return

            is_reversible = is_reversible_from_text
            reactant_names, product_names = components
            reactant_details = [self.db_updater._resolve_molecule_name_from_sl_text(name) for name in reactant_names]
            product_details = [self.db_updater._resolve_molecule_name_from_sl_text(name) for name in product_names]

            if any(m is None for m in reactant_details + product_details):
                print(
                    f"    - CRITICAL ERROR: Could not resolve all molecules in '{reaction_text}'. Cannot create pairs.")
                return

            reactant_ids = [m['molecule_id'] for m in reactant_details]
            product_ids = [m['molecule_id'] for m in product_details]

            # Use the imported utility function to get all pairs correctly
            combinations = generate_reaction_combinations(reactant_ids, product_ids, is_reversible)
            # Convert to the list-of-dicts format the rest of the function expects
            molecule_pairs = [{'reactant_molecule_id': r, 'product_molecule_id': p} for r, p in combinations]

        # --- Part 3: Create the new links and their corresponding pairs ---
        for enzyme_id in new_enzyme_ids_to_link:
            # Use the robust helper method from db_updater
            re_id = self.db_updater._ensure_reaction_enzyme(gen_reaction_id, enzyme_id)

            if re_id:
                # _ensure_reaction_enzyme already logs if the link was new or existing.
                # We only need to create pairs for the link (new or not).
                # The _ensure_reaction_pairs method handles duplicates gracefully.
                if molecule_pairs:
                    for pair in molecule_pairs:
                        reactant_id = pair['reactant_molecule_id']
                        product_id = pair['product_molecule_id']
                        # is_sl_reversible=False because we are creating specific pairs one-by-one, not generating new combinations
                        self.db_updater._ensure_reaction_pairs(re_id, [reactant_id], [product_id],
                                                               is_sl_reversible=False)
                    print(f"      - Ensured {len(molecule_pairs)} molecular pairs exist for link re_id: {re_id}.")
                else:
                    # This case is now handled by the check above, but we keep the log for safety.
                    print(
                        f"    - CRITICAL WARNING: No molecular pairs found or derived for reaction {gen_reaction_id}. Cannot create pairs for new link re_id: {re_id}.")

    def run_refactor_workflow(self):
        """Runs an interactive workflow to correct an improper generalization link."""
        print("\n--- Interactive Reaction Refactoring Workflow ---")
        try:
            # Get user input
            child_id_str = input("Enter the ID of the SPECIFIC reaction to move (e.g., 393): ")
            parent_id_str = input("Enter the ID of the INCORRECT GENERALIZED parent to unlink from (e.g., 2796): ")
            new_parent_text = input(
                "Enter the text for the NEW, correct generalized reaction (e.g., 'Ceramide => Ceramide'): ")

            child_id = int(child_id_str)
            parent_id = int(parent_id_str)

            with self.conn.transaction():
                # Step 1: Unlink from incorrect parent
                print(f"\n1. Unlinking child {child_id} from parent {parent_id}...")
                self.cursor.execute(
                    "SELECT generalized_from_reaction_ids FROM lipograph.reactions WHERE reaction_id = %s",
                    (parent_id,))
                result = self.cursor.fetchone()
                if not result or not result['generalized_from_reaction_ids']:
                    raise ValueError(f"Parent {parent_id} not found or has no links.")

                links_str = result['generalized_from_reaction_ids']
                child_ids = {int(x) for x in links_str.strip('|').split('|') if x}
                if child_id not in child_ids:
                    raise ValueError(f"Child {child_id} is not linked to parent {parent_id}.")

                child_ids.remove(child_id)
                new_links_str = f"|{'|'.join(map(str, sorted(list(child_ids))))}|" if child_ids else ""
                self.cursor.execute(
                    "UPDATE lipograph.reactions SET generalized_from_reaction_ids = %s WHERE reaction_id = %s",
                    (new_links_str, parent_id))
                print(f"   ...Unlinked successfully.")

                # Step 2: Create the new, correct generalized reaction
                print(f"\n2. Creating new generalized reaction '{new_parent_text}'...")
                initial_links = f"|{child_id}|"
                self.cursor.execute(
                    "INSERT INTO lipograph.reactions (reaction_text, generalized_from_reaction_ids) VALUES (%s, %s) RETURNING reaction_id;",
                    (new_parent_text, initial_links)
                )
                new_parent_id = self.cursor.fetchone()['reaction_id']
                print(f"   ...New parent created with ID: {new_parent_id}")

                # Step 3: Propagate enzymes and pairs to the new parent
                print(f"\n3. Propagating enzymes and pairs to new parent {new_parent_id}...")
                self.propagate_enzymes_and_pairs(new_parent_id)

                # Step 4: Ask to permanently exclude the incorrect link
                if input(
                        f"\n4. Permanently exclude the incorrect link ({parent_id}:{child_id}) from future runs? (yes/no): ").lower() == 'yes':
                    exclusions = self.load_exclusions()
                    exclusions.add((parent_id, child_id))
                    self.save_exclusions(exclusions)

                confirm = input("\nReview the changes above. Commit them to the database? (yes/no): ").lower()
                if confirm != 'yes':
                    raise RuntimeError("User aborted. Rolling back transaction.")

            print("\n✅ Transaction successfully committed.")

        except (ValueError, RuntimeError, psycopg.Error) as e:
            print(f"\n❌ An error occurred: {e}. Transaction has been rolled back.")

    def load_exclusions(self) -> Set[Tuple[int, int]]:
        """Loads a set of (gen_id, spec_id) tuples from the exclusion file."""
        exclusion_file = Path("../scripts/generalization_exclusions.json")
        if not exclusion_file.exists():
            return set()
        try:
            with open(exclusion_file, 'r') as f:
                # The file stores lists, but we convert them to tuples to be hashable for a set
                excluded_pairs = [tuple(pair) for pair in json.load(f)]
            print(f"Loaded {len(excluded_pairs)} exclusions from {exclusion_file.name}")
            return set(excluded_pairs)
        except (json.JSONDecodeError, TypeError) as e:
            print(f"Warning: Could not read or parse exclusion file. Starting with none. Error: {e}")
            return set()

    def save_exclusions(self, excluded_pairs: Set[Tuple[int, int]]):
        """Saves the set of excluded pairs back to the JSON file."""
        exclusion_file = Path("../scripts/generalization_exclusions.json")
        try:
            # Convert set of tuples back to a list of lists for JSON serialization
            with open(exclusion_file, 'w') as f:
                json.dump(sorted(list(excluded_pairs)), f, indent=2)  # Sort for consistent file output
            print(f"Saved {len(excluded_pairs)} exclusions to {exclusion_file.name}")
        except Exception as e:
            print(f"Error: Could not save exclusion file. Error: {e}")

    def run_propagation_workflow(self):
        """Runs enzyme and pair propagation for ALL existing generalized reactions."""
        print("\n--- Full Enzyme Propagation Workflow ---")

        self.cursor.execute(
            "SELECT reaction_id FROM lipograph.reactions WHERE generalized_from_reaction_ids IS NOT NULL AND generalized_from_reaction_ids != ''"
        )
        all_linked_gen_ids = [row['reaction_id'] for row in self.cursor.fetchall()]

        if not all_linked_gen_ids:
            print("No existing generalized reactions found to propagate.")
            return

        print(f"Found {len(all_linked_gen_ids)} generalized reactions to process.")
        if input("Proceed with propagation? (yes/no): ").lower() != 'yes':
            print("Aborting.")
            return

        try:
            with self.conn.transaction():
                for gen_id in all_linked_gen_ids:
                    self.propagate_enzymes_and_pairs(gen_id)
            print("\n✅ All enzyme propagation changes have been successfully committed.")
        except Exception as e:
            print(f"\n❌ An error occurred during propagation. Transaction was rolled back. Error: {e}", file=sys.stderr)

    def run_linking_workflow(self):
        """Main execution method with persistent exclusions and a two-phase update process."""
        # --- Optional Debug Mode ---
        if input("Enter debug mode? (yes/no): ").lower() == 'yes':
            while True:
                action = input("Enter 'debug [ID]', 'compare [ID1] [ID2]', or 'quit': ").strip()
                if action.lower() == 'quit':
                    break
                parts = action.split()
                if len(parts) == 2 and parts[0].lower() == 'debug':
                    try:
                        self.debug_specific_reaction(int(parts[1]))
                        self.check_existing_links(int(parts[1]))
                    except ValueError:
                        print("Please enter a valid reaction ID number")
                elif len(parts) == 3 and parts[0].lower() == 'compare':
                    try:
                        self.compare_reactions(int(parts[1]), int(parts[2]))
                    except ValueError:
                        print("Please enter valid reaction ID numbers")
                else:
                    print("Usage: 'debug [ID]', 'compare [ID1] [ID2]', or 'quit'")
            if input("Continue with main process? (yes/no): ").lower() != 'yes':
                print("Exiting after debug session.")
                return

        # --- Phase 1: Find and Apply Generalization Links ---
        print("\n--- Phase 1: Finding and Applying Generalization Links ---")

        # Load existing exclusions first
        exclusions = self.load_exclusions()

        # Find all potential new links
        generalized_reactions_ids, specific_reactions_ids = self.fetch_reaction_ids()
        all_reaction_ids = generalized_reactions_ids + specific_reactions_ids
        reaction_class_cache = self.build_reaction_classes_cache(all_reaction_ids)

        updates_to_perform = collections.defaultdict(list)
        for gen_id in generalized_reactions_ids:
            if gen_id in reaction_class_cache:
                self.cursor.execute(
                    "SELECT generalized_from_reaction_ids FROM lipograph.reactions WHERE reaction_id = %s", (gen_id,))
                existing_links_str = (self.cursor.fetchone() or {}).get('generalized_from_reaction_ids', "") or ""
                existing_ids = {int(x) for x in existing_links_str.strip('|').split('|') if x}

                gen_reactant_classes, gen_product_classes = reaction_class_cache[gen_id]
                for spec_id in specific_reactions_ids:
                    if spec_id in existing_ids or (gen_id, spec_id) in exclusions:
                        continue

                    if spec_id in reaction_class_cache:
                        spec_reactant_classes, spec_product_classes = reaction_class_cache[spec_id]
                        if gen_reactant_classes == spec_reactant_classes and gen_product_classes == spec_product_classes:
                            updates_to_perform[gen_id].append(spec_id)

        if not updates_to_perform:
            print("\nNo new, un-excluded generalization links found.")
        else:
            print(f"\nFound {sum(len(v) for v in updates_to_perform.values())} new potential links to create.")

            # Generate the summary report for review
            report_data = []
            all_ids_to_fetch_text = list(updates_to_perform.keys())
            for spec_ids in updates_to_perform.values():
                all_ids_to_fetch_text.extend(spec_ids)
            reaction_texts = self.get_reaction_texts(all_ids_to_fetch_text)
            for gen_id, spec_ids in updates_to_perform.items():
                for spec_id in sorted(spec_ids):
                    report_data.append(
                        {"generalized_reaction_id": gen_id, "generalized_reaction_text": reaction_texts.get(gen_id),
                         "matched_specific_reaction_id": spec_id,
                         "specific_reaction_text": reaction_texts.get(spec_id)})

            report_df = pd.DataFrame(report_data)
            report_path = PROJECT_ROOT / "proposed_generalization_links.tsv"
            report_df.to_csv(report_path, sep="\t", index=False)
            print(f"Summary report of NEW links saved to: {report_path.resolve()}")

            user_input_links = input("\nApply these new generalization links? (yes/no/edit): ").lower()

            if user_input_links == "edit":
                while True:
                    remove_input = input("Enter pair to remove ('gen_id:spec_id' or 'gen_id') or 'done': ").strip()
                    if remove_input.lower() == 'done': break
                    try:
                        if ':' in remove_input:
                            gen_str, spec_str = remove_input.split(':')
                            gen_id, spec_id = int(gen_str), int(spec_str)
                            if gen_id in updates_to_perform and spec_id in updates_to_perform[gen_id]:
                                updates_to_perform[gen_id].remove(spec_id)
                                if not updates_to_perform[gen_id]: del updates_to_perform[gen_id]
                                print(f"Removed link: {gen_id} -> {spec_id}")
                                if input(f"  Permanently exclude this link? (yes/no): ").lower() == 'yes':
                                    exclusions.add((gen_id, spec_id))
                        else:
                            gen_id = int(remove_input)
                            if gen_id in updates_to_perform:
                                removed_specs = updates_to_perform.pop(gen_id)
                                print(f"Removed all new links for generalized reaction {gen_id}")
                                if input(
                                        f"  Permanently exclude all {len(removed_specs)} of these links? (yes/no): ").lower() == 'yes':
                                    for spec_id in removed_specs:
                                        exclusions.add((gen_id, spec_id))
                    except (ValueError, KeyError):
                        print("Invalid format or pair not found.")

                self.save_exclusions(exclusions)

            if user_input_links in ("yes", "edit"):
                if not updates_to_perform:
                    print("No new links remaining after edit. Skipping update.")
                else:
                    try:
                        with self.conn.transaction():
                            print("\nApplying new generalization links...")
                            update_query = "UPDATE lipograph.reactions SET generalized_from_reaction_ids = %s WHERE reaction_id = %s;"
                            for gen_id, spec_ids_to_add in updates_to_perform.items():
                                if not spec_ids_to_add: continue
                                self.cursor.execute(
                                    "SELECT generalized_from_reaction_ids FROM lipograph.reactions WHERE reaction_id = %s",
                                    (gen_id,))
                                existing = (self.cursor.fetchone() or {}).get('generalized_from_reaction_ids', "") or ""
                                existing_ids = {int(x) for x in existing.strip('|').split('|') if x}
                                all_ids = existing_ids.union(set(spec_ids_to_add))
                                id_string = f"|{'|'.join(map(str, sorted(all_ids)))}|"
                                self.cursor.execute(update_query, (id_string, gen_id))
                        print("✅ Successfully committed new generalization links.")
                    except Exception as e:
                        print(f"❌ Error applying links. Transaction rolled back. Error: {e}");
                        return
            else:
                print("Skipping generalization linking.")

        # --- Phase 2: Propagate Enzyme Connections ---
        print("\n" + "=" * 50)
        print("--- Phase 2: Propagating Enzyme Connections ---")
        propagate_choice = input("Choose propagation scope ('updated', 'all', 'no'): ").lower()

        if propagate_choice not in ("updated", "all"):
            print("Skipping enzyme propagation. Workflow complete.");
            return

        ids_to_propagate = []
        if propagate_choice == "updated":
            ids_to_propagate = list(updates_to_perform.keys())
            if not ids_to_propagate:
                print("No reactions were updated, so there's nothing to propagate for 'updated'.");
                return
            print(f"Will propagate enzymes for {len(ids_to_propagate)} newly updated generalized reactions.")
        elif propagate_choice == "all":
            self.cursor.execute(
                "SELECT reaction_id FROM lipograph.reactions WHERE generalized_from_reaction_ids IS NOT NULL AND generalized_from_reaction_ids != ''")
            ids_to_propagate = [row['reaction_id'] for row in self.cursor.fetchall()]
            if not ids_to_propagate:
                print("No generalized reactions found in the database to propagate.");
                return
            print(f"Will propagate enzymes for all {len(ids_to_propagate)} generalized reactions in the database.")

        try:
            with self.conn.transaction():
                for gen_id in ids_to_propagate:
                    self.propagate_enzymes_and_pairs(gen_id)
            print("\n✅ All enzyme propagation changes have been successfully committed.")
        except Exception as e:
            print(f"\n❌ An error occurred during propagation. Transaction was rolled back. Error: {e}", file=sys.stderr)

    def run(self):
        """Main execution method to choose a workflow."""
        print("\n--- Generalization Linker and Refactoring Tool ---")
        mode = input("Choose a mode: [l]ink new, [p]ropagate all, [r]efactor, or [q]uit: ").lower()

        if mode == 'l':
            self.run_linking_workflow()
        elif mode == 'p':
            self.run_propagation_workflow()
        elif mode == 'r':
            self.run_refactor_workflow()
        elif mode == 'q':
            print("Exiting.")
        else:
            print("Invalid mode selected.")


def main():
    config_path = Path("../config/lipid_config.yaml")
    if not config_path.exists():
        print(f"Error: Configuration file not found at {config_path}", file=sys.stderr)
        sys.exit(1)

    config = ConfigLoader(str(config_path))
    linker = GeneralizationLinker(config)

    try:
        linker.connect_and_setup()
        linker.run()
    finally:
        linker.close()


if __name__ == "__main__":
    main()