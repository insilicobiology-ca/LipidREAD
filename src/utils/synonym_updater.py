#!/usr/bin/env python3
"""
Interactive Synonym Updater CLI Tool

An interactive command-line interface for managing molecule synonyms in the LipidCRED database.
Provides menu-driven interface with confirmation prompts before making changes.
Uses psycopg3.

Usage:
    python interactive_synonym_cli.py
"""

import logging
import sys
from pathlib import Path
from typing import List, Dict, Optional, Any
import json
import csv

# Add the project root to the path so we can import our modules
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from src.data.data_access import DatabaseLipidDataAccess
from src.data.config_loader import ConfigLoader

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class InteractiveSynonymUpdater:
    """
    Interactive synonym updater with psycopg3 support.
    """

    def __init__(self, data_access):
        """Initialize with DatabaseLipidDataAccess instance."""
        self.data_access = data_access
        self.data_access.conn.autocommit = True
        self.conn = data_access.conn
        self.cursor = data_access.cursor

    def add_synonyms_to_molecule(
        self, molecule_id: int, new_synonyms: List[str]
    ) -> bool:
        """Add new synonyms to a specific molecule."""
        if not new_synonyms:
            logger.warning(f"No synonyms provided for molecule_id {molecule_id}")
            return False

        try:
            # Get current synonyms
            current_synonyms = self._get_current_synonyms(molecule_id)

            # Merge with new synonyms (avoiding duplicates)
            updated_synonyms = self._merge_synonyms(current_synonyms, new_synonyms)

            # Start transaction and update
            with self.conn.transaction():
                self._update_molecule_synonyms(molecule_id, updated_synonyms)
                # Transaction auto-commits here

            logger.info(f"Successfully updated synonyms for molecule_id {molecule_id}")
            return True

        except Exception as e:
            logger.error(
                f"Failed to update synonyms for molecule_id {molecule_id}: {e}"
            )
            return False

    def _execute_with_data_access_transaction(
        self, operation: str, molecule_id: int, synonyms: List[str]
    ) -> bool:
        """Fallback method using data_access methods directly."""
        try:
            # Get current synonyms
            current_synonyms = self._get_current_synonyms(molecule_id)

            if operation == "add":
                updated_synonyms = self._merge_synonyms(current_synonyms, synonyms)
            elif operation == "replace":
                updated_synonyms = synonyms
            elif operation == "remove":
                updated_synonyms = [
                    syn for syn in current_synonyms if syn not in synonyms
                ]
            else:
                raise ValueError(f"Unknown operation: {operation}")

            # Update using data_access
            synonyms_string = self._format_synonyms_string(updated_synonyms)

            # Try to use data_access methods for update
            if hasattr(self.data_access, "execute_query"):
                query = "UPDATE lipograph.molecules SET molecule_synonyms = %s WHERE molecule_id = %s"
                self.data_access.execute_query(query, (synonyms_string, molecule_id))
            else:
                # Fallback: use cursor directly
                cursor = self.cursor or self.data_access.cursor
                cursor.execute(
                    "UPDATE lipograph.molecules SET molecule_synonyms = %s WHERE molecule_id = %s",
                    (synonyms_string, molecule_id),
                )

                # Try to commit if connection available
                if hasattr(self.data_access, "commit"):
                    self.data_access.commit()
                elif hasattr(cursor, "connection"):
                    cursor.connection.commit()

            logger.info(f"Successfully updated synonyms for molecule_id {molecule_id}")
            return True

        except Exception as e:
            logger.error(
                f"Failed to update synonyms for molecule_id {molecule_id}: {e}"
            )
            return False

    def batch_update_synonyms(
        self, synonym_updates: Dict[int, List[str]]
    ) -> Dict[str, int]:
        """Batch update synonyms for multiple molecules."""
        results = {"success": 0, "failed": 0, "skipped": 0}

        try:
            # Single transaction for the entire batch (psycopg3 style)
            with self.connection.transaction():
                for molecule_id, new_synonyms in synonym_updates.items():
                    try:
                        if not new_synonyms:
                            results["skipped"] += 1
                            continue

                        # Get current synonyms
                        current_synonyms = self._get_current_synonyms(molecule_id)

                        # Merge with new synonyms
                        updated_synonyms = self._merge_synonyms(
                            current_synonyms, new_synonyms
                        )

                        # Only update if there are actual changes
                        if set(updated_synonyms) != set(current_synonyms):
                            self._update_molecule_synonyms(
                                molecule_id, updated_synonyms
                            )
                            results["success"] += 1
                            logger.debug(
                                f"Updated synonyms for molecule_id {molecule_id}"
                            )
                        else:
                            results["skipped"] += 1
                            logger.debug(
                                f"No changes needed for molecule_id {molecule_id}"
                            )

                    except Exception as e:
                        logger.error(f"Failed to update molecule_id {molecule_id}: {e}")
                        results["failed"] += 1

        except Exception as e:
            logger.error(f"Batch update transaction failed: {e}")
            results["failed"] = len(synonym_updates)
            results["success"] = 0

        logger.info(f"Batch update completed: {results}")
        return results

    def replace_synonyms_for_molecule(
        self, molecule_id: int, new_synonyms: List[str]
    ) -> bool:
        """Replace all synonyms for a molecule."""
        try:
            with self.conn.transaction():
                self._update_molecule_synonyms(molecule_id, new_synonyms)
                # Transaction auto-commits here

            logger.info(f"Successfully replaced synonyms for molecule_id {molecule_id}")
            return True

        except Exception as e:
            logger.error(
                f"Failed to replace synonyms for molecule_id {molecule_id}: {e}"
            )
            return False

    def remove_synonyms_from_molecule(
        self, molecule_id: int, synonyms_to_remove: List[str]
    ) -> bool:
        """Remove specific synonyms from a molecule."""
        try:
            # Get current synonyms first (outside transaction)
            current_synonyms = self._get_current_synonyms(molecule_id)

            # Remove specified synonyms
            updated_synonyms = [
                syn for syn in current_synonyms if syn not in synonyms_to_remove
            ]

            # Update in transaction
            with self.conn.transaction():
                self._update_molecule_synonyms(molecule_id, updated_synonyms)
                # Transaction auto-commits here

            logger.info(f"Successfully removed synonyms from molecule_id {molecule_id}")
            return True

        except Exception as e:
            logger.error(
                f"Failed to remove synonyms from molecule_id {molecule_id}: {e}"
            )
            return False

    def get_molecules_with_synonyms(self, limit: Optional[int] = None) -> List[Dict]:
        """Get all molecules that have synonyms."""
        try:
            query = """
                SELECT molecule_id, molecule_name, molecule_synonyms, cleaned_abbreviation
                FROM lipograph.molecules 
                WHERE molecule_synonyms IS NOT NULL 
                AND molecule_synonyms != ''
                ORDER BY molecule_id
            """

            if limit:
                query += f" LIMIT {limit}"

            self.cursor.execute(query)
            results = self.cursor.fetchall()

            # Convert to dict and parse synonyms
            formatted_results = []
            for result in results:
                result_dict = dict(result)
                result_dict["parsed_synonyms"] = self._parse_synonyms_string(
                    result_dict["molecule_synonyms"]
                )
                formatted_results.append(result_dict)

            return formatted_results

        except Exception as e:
            logger.error(f"Failed to get molecules with synonyms: {e}")
            return []

    def search_molecules_by_synonym(self, synonym: str) -> List[Dict]:
        """Search for molecules that have a specific synonym."""
        try:
            # Use the pipe format to ensure exact matches only
            search_pattern = f"|{synonym}|"

            query = """
                SELECT molecule_id, molecule_name, cleaned_abbreviation, molecule_synonyms
                FROM lipograph.molecules 
                WHERE molecule_synonyms LIKE %s
                ORDER BY molecule_id
            """

            self.cursor.execute(query, (f"%{search_pattern}%",))
            results = self.cursor.fetchall()

            # Convert to dict and parse synonyms
            formatted_results = []
            for result in results:
                result_dict = dict(result)
                result_dict["parsed_synonyms"] = self._parse_synonyms_string(
                    result_dict["molecule_synonyms"]
                )
                formatted_results.append(result_dict)

            logger.info(
                f"Found {len(formatted_results)} molecules with synonym '{synonym}'"
            )
            return formatted_results

        except Exception as e:
            logger.error(f"Failed to search for synonym '{synonym}': {e}")
            return []

    def get_molecule_info(self, molecule_id: int) -> Optional[Dict]:
        """Get detailed information about a specific molecule."""
        try:
            query = """
                SELECT molecule_id, molecule_name, cleaned_abbreviation, 
                       molecule_synonyms, swisslipids_id, sl_lipid_class
                FROM lipograph.molecules 
                WHERE molecule_id = %s
            """

            self.cursor.execute(query, (molecule_id,))
            result = self.cursor.fetchone()

            if result:
                result_dict = dict(result)
                result_dict["parsed_synonyms"] = self._parse_synonyms_string(
                    result_dict["molecule_synonyms"]
                )
                return result_dict

            return None

        except Exception as e:
            logger.error(f"Failed to get molecule info for ID {molecule_id}: {e}")
            return None

    def validate_molecule_exists(self, molecule_id: int) -> bool:
        """Check if a molecule exists in the database."""
        query = "SELECT 1 FROM lipograph.molecules WHERE molecule_id = %s"
        self.cursor.execute(query, (molecule_id,))
        return self.cursor.fetchone() is not None

    def _get_current_synonyms(self, molecule_id: int) -> List[str]:
        """Get current synonyms for a molecule."""
        query = (
            "SELECT molecule_synonyms FROM lipograph.molecules WHERE molecule_id = %s"
        )
        self.cursor.execute(query, (molecule_id,))
        result = self.cursor.fetchone()

        if result and result["molecule_synonyms"]:
            return self._parse_synonyms_string(result["molecule_synonyms"])
        return []

    def _parse_synonyms_string(self, synonyms_string: Optional[str]) -> List[str]:
        """Parse the pipe-delimited synonyms string into a list."""
        if not synonyms_string:
            return []

        cleaned = synonyms_string.strip("|")
        if not cleaned:
            return []

        return [syn.strip() for syn in cleaned.split("|") if syn.strip()]

    def _format_synonyms_string(self, synonyms: List[str]) -> str:
        """Format synonyms list into pipe-delimited string."""
        if not synonyms:
            return ""

        clean_synonyms = [syn.strip() for syn in synonyms if syn.strip()]
        unique_synonyms = list(dict.fromkeys(clean_synonyms))

        return f"|{'|'.join(unique_synonyms)}|"

    def _merge_synonyms(
        self, current_synonyms: List[str], new_synonyms: List[str]
    ) -> List[str]:
        """Merge current and new synonyms, avoiding duplicates."""
        all_synonyms = current_synonyms + new_synonyms

        seen = set()
        merged = []
        for syn in all_synonyms:
            clean_syn = syn.strip()
            if clean_syn and clean_syn not in seen:
                seen.add(clean_syn)
                merged.append(clean_syn)

        return merged

    def _update_molecule_synonyms(self, molecule_id: int, synonyms: List[str]) -> None:
        """Update the molecule_synonyms field for a specific molecule."""
        synonyms_string = self._format_synonyms_string(synonyms)

        query = """
            UPDATE lipograph.molecules 
            SET molecule_synonyms = %s 
            WHERE molecule_id = %s
        """

        self.cursor.execute(query, (synonyms_string, molecule_id))

        if self.cursor.rowcount == 0:
            raise ValueError(f"No molecule found with ID {molecule_id}")

        logger.debug(f"Updated molecule_id {molecule_id} with {len(synonyms)} synonyms")


class InteractiveSynonymCLI:
    """Interactive command-line interface for synonym management."""

    def __init__(self, config_path: Optional[Path] = None):
        """Initialize with database connection."""
        if config_path is None:
            config_path = project_root / "config" / "lipid_config.yaml"

        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        print(f"🔗 Connecting to database using config: {config_path}")

        # Initialize components
        self.config_loader = ConfigLoader(config_path)
        self.data_access = DatabaseLipidDataAccess(self.config_loader)
        self.synonym_updater = InteractiveSynonymUpdater(self.data_access)

        print("✅ Connected to database successfully!")

    def run(self):
        """Main interactive loop."""
        print("\n" + "=" * 60)
        print("🧬 LipidCRED Synonym Manager")
        print("=" * 60)
        print("Interactive tool for managing molecule synonyms")
        print("All changes require confirmation before committing to database")
        print("=" * 60)

        while True:
            try:
                choice = self.show_main_menu()

                if choice == "1":
                    self.add_synonyms_mode()
                elif choice == "2":
                    self.replace_synonyms_mode()
                elif choice == "3":
                    self.remove_synonyms_mode()
                elif choice == "4":
                    self.search_mode()
                elif choice == "5":
                    self.view_molecule_mode()
                elif choice == "6":
                    self.list_molecules_mode()
                elif choice == "7":
                    self.batch_mode()
                elif choice == "8":
                    self.export_mode()
                elif choice == "0":
                    if self.confirm_action("Exit the program"):
                        print("👋 Goodbye!")
                        break
                else:
                    print("❌ Invalid choice. Please try again.")

            except KeyboardInterrupt:
                print("\n\n⚠️  Interrupted by user")
                if self.confirm_action("Exit the program"):
                    break
            except Exception as e:
                print(f"❌ Error: {e}")
                logger.exception("Unexpected error in main loop")

    def show_main_menu(self) -> str:
        """Display main menu and get user choice."""
        print("\n" + "-" * 50)
        print("📋 MAIN MENU")
        print("-" * 50)
        print("1. 🔧 Add synonyms to molecule")
        print("2. 🔄 Replace all synonyms for molecule")
        print("3. 🗑️  Remove synonyms from molecule")
        print("4. 🔍 Search molecules by synonym")
        print("5. 👁️  View molecule details")
        print("6. 📜 List molecules with synonyms")
        print("7. 📦 Batch operations")
        print("8. 💾 Export synonyms")
        print("0. 🚪 Exit")
        print("-" * 50)

        return input("Choose an option (0-8): ").strip()

    def add_synonyms_mode(self):
        """Interactive mode for adding synonyms."""
        print("\n🔧 ADD SYNONYMS MODE")
        print("-" * 30)

        # Get molecule ID
        molecule_id = self.get_molecule_id()
        if not molecule_id:
            return

        # Show current molecule info
        self.display_molecule_info(molecule_id)

        # Get new synonyms
        synonyms = self.get_synonyms_input("Enter synonyms to ADD (comma-separated)")
        if not synonyms:
            return

        # Preview changes
        current_synonyms = self.synonym_updater._get_current_synonyms(molecule_id)
        merged_synonyms = self.synonym_updater._merge_synonyms(
            current_synonyms, synonyms
        )

        print("\n📋 PREVIEW CHANGES:")
        print(f"Current synonyms: {current_synonyms}")
        print(f"Adding: {synonyms}")
        print(f"Result will be: {merged_synonyms}")

        # Confirm and execute
        if self.confirm_action("add these synonyms"):
            success = self.synonym_updater.add_synonyms_to_molecule(
                molecule_id, synonyms
            )
            if success:
                print("✅ Synonyms added successfully!")
            else:
                print("❌ Failed to add synonyms")

    def replace_synonyms_mode(self):
        """Interactive mode for replacing synonyms."""
        print("\n🔄 REPLACE SYNONYMS MODE")
        print("-" * 30)

        molecule_id = self.get_molecule_id()
        if not molecule_id:
            return

        self.display_molecule_info(molecule_id)

        synonyms = self.get_synonyms_input(
            "Enter NEW synonyms to REPLACE all existing (comma-separated)"
        )
        if synonyms is None:  # User cancelled
            return

        # Preview changes
        current_synonyms = self.synonym_updater._get_current_synonyms(molecule_id)

        print("\n📋 PREVIEW CHANGES:")
        print(f"Current synonyms: {current_synonyms}")
        print(f"Will be REPLACED with: {synonyms}")

        if self.confirm_action("replace ALL synonyms"):
            success = self.synonym_updater.replace_synonyms_for_molecule(
                molecule_id, synonyms
            )
            if success:
                print("✅ Synonyms replaced successfully!")
            else:
                print("❌ Failed to replace synonyms")

    def remove_synonyms_mode(self):
        """Interactive mode for removing synonyms."""
        print("\n🗑️ REMOVE SYNONYMS MODE")
        print("-" * 30)

        molecule_id = self.get_molecule_id()
        if not molecule_id:
            return

        self.display_molecule_info(molecule_id)

        current_synonyms = self.synonym_updater._get_current_synonyms(molecule_id)
        if not current_synonyms:
            print("ℹ️  This molecule has no synonyms to remove.")
            return

        print(f"\nCurrent synonyms: {current_synonyms}")
        synonyms = self.get_synonyms_input("Enter synonyms to REMOVE (comma-separated)")
        if not synonyms:
            return

        # Check if synonyms exist
        to_remove = [s for s in synonyms if s in current_synonyms]
        not_found = [s for s in synonyms if s not in current_synonyms]

        if not_found:
            print(f"⚠️  These synonyms were not found: {not_found}")

        if not to_remove:
            print("❌ No valid synonyms to remove")
            return

        # Preview changes
        remaining = [s for s in current_synonyms if s not in to_remove]

        print("\n📋 PREVIEW CHANGES:")
        print(f"Current synonyms: {current_synonyms}")
        print(f"Will remove: {to_remove}")
        print(f"Remaining: {remaining}")

        if self.confirm_action("remove these synonyms"):
            success = self.synonym_updater.remove_synonyms_from_molecule(
                molecule_id, to_remove
            )
            if success:
                print("✅ Synonyms removed successfully!")
            else:
                print("❌ Failed to remove synonyms")

    def search_mode(self):
        """Interactive mode for searching by synonym."""
        print("\n🔍 SEARCH MODE")
        print("-" * 30)

        synonym = input("Enter synonym to search for: ").strip()
        if not synonym:
            return

        print(f"🔍 Searching for molecules with synonym '{synonym}'...")
        results = self.synonym_updater.search_molecules_by_synonym(synonym)

        if not results:
            print("❌ No molecules found with that synonym")
            return

        print(f"\n✅ Found {len(results)} molecules:")
        print("=" * 80)

        for i, result in enumerate(results, 1):
            print(f"{i}. ID: {result['molecule_id']}")
            print(f"   Name: {result['molecule_name']}")
            print(f"   Abbreviation: {result['cleaned_abbreviation']}")
            print(f"   All synonyms: {result['parsed_synonyms']}")
            print("-" * 80)

    def view_molecule_mode(self):
        """Interactive mode for viewing molecule details."""
        print("\n👁️ VIEW MOLECULE MODE")
        print("-" * 30)

        molecule_id = self.get_molecule_id()
        if not molecule_id:
            return

        self.display_molecule_info(molecule_id, detailed=True)

    def list_molecules_mode(self):
        """Interactive mode for listing molecules with synonyms."""
        print("\n📜 LIST MOLECULES MODE")
        print("-" * 30)

        limit_input = input("Enter limit (press Enter for no limit): ").strip()
        limit = None
        if limit_input:
            try:
                limit = int(limit_input)
            except ValueError:
                print("❌ Invalid limit. Using no limit.")

        print("📜 Loading molecules with synonyms...")
        results = self.synonym_updater.get_molecules_with_synonyms(limit)

        if not results:
            print("❌ No molecules found with synonyms")
            return

        print(f"\n✅ Found {len(results)} molecules with synonyms:")
        print("=" * 100)

        for i, result in enumerate(results, 1):
            synonyms_count = len(result["parsed_synonyms"])
            print(
                f"{i:>3}. ID: {result['molecule_id']:>6} | "
                f"Name: {result['molecule_name'][:40]:<40} | "
                f"Synonyms ({synonyms_count}): {', '.join(result['parsed_synonyms'][:3])}"
                f"{'...' if synonyms_count > 3 else ''}"
            )

    def batch_mode(self):
        """Interactive mode for batch operations."""
        print("\n📦 BATCH OPERATIONS MODE")
        print("-" * 30)
        print("🚧 Coming soon! For now, use individual operations.")

    def export_mode(self):
        """Interactive mode for exporting synonyms."""
        print("\n💾 EXPORT MODE")
        print("-" * 30)

        filename = input("Enter filename (e.g., synonyms.csv): ").strip()
        if not filename:
            filename = "synonyms_export.csv"

        output_path = Path(filename)

        print("💾 Exporting synonyms...")
        results = self.synonym_updater.get_molecules_with_synonyms()

        try:
            with open(output_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(
                    ["molecule_id", "molecule_name", "cleaned_abbreviation", "synonyms"]
                )

                for result in results:
                    synonyms_str = "|".join(result["parsed_synonyms"])
                    writer.writerow(
                        [
                            result["molecule_id"],
                            result["molecule_name"],
                            result["cleaned_abbreviation"],
                            synonyms_str,
                        ]
                    )

            print(f"✅ Exported {len(results)} molecules to {output_path}")

        except Exception as e:
            print(f"❌ Export failed: {e}")

    def get_molecule_id(self) -> Optional[int]:
        """Get and validate molecule ID from user."""
        while True:
            try:
                molecule_input = input(
                    "Enter molecule ID (or 'back' to return): "
                ).strip()
                if molecule_input.lower() == "back":
                    return None

                molecule_id = int(molecule_input)

                if not self.synonym_updater.validate_molecule_exists(molecule_id):
                    print(f"❌ Molecule ID {molecule_id} does not exist")
                    continue

                return molecule_id

            except ValueError:
                print("❌ Please enter a valid number")
            except KeyboardInterrupt:
                return None

    def get_synonyms_input(self, prompt: str) -> Optional[List[str]]:
        """Get synonyms input from user."""
        try:
            synonyms_input = input(f"{prompt}: ").strip()
            if not synonyms_input:
                return []

            # Split by comma and clean
            synonyms = [s.strip() for s in synonyms_input.split(",") if s.strip()]
            return synonyms

        except KeyboardInterrupt:
            return None

    def display_molecule_info(self, molecule_id: int, detailed: bool = False):
        """Display molecule information."""
        info = self.synonym_updater.get_molecule_info(molecule_id)
        if not info:
            print(f"❌ Could not retrieve info for molecule {molecule_id}")
            return

        print("\n" + "=" * 60)
        print(f"🧬 MOLECULE INFO: {molecule_id}")
        print("=" * 60)
        print(f"Name: {info['molecule_name']}")
        print(f"Abbreviation: {info['cleaned_abbreviation']}")
        print(f"Current synonyms: {info['parsed_synonyms']}")

        if detailed:
            print(f"SwissLipids ID: {info['swisslipids_id']}")
            print(f"Lipid Class: {info['sl_lipid_class']}")

        print("=" * 60)

    def confirm_action(self, action: str) -> bool:
        """Get confirmation from user before proceeding."""
        while True:
            try:
                response = (
                    input(f"\n❓ Are you sure you want to {action}? (y/N): ")
                    .strip()
                    .lower()
                )
                if response in ["y", "yes"]:
                    return True
                elif response in ["n", "no", ""]:
                    return False
                else:
                    print("Please enter 'y' for yes or 'n' for no")
            except KeyboardInterrupt:
                return False

    def close(self):
        """Clean up database connections."""
        # DatabaseLipidDataAccess uses __del__ for cleanup, so no explicit close needed
        pass


def main():
    """Main entry point."""
    try:
        # Check for custom config path
        config_path = Path("../../config/lipid_config.yaml")
        if len(sys.argv) > 1:
            config_path = Path(sys.argv[1])

        # Initialize and run CLI
        cli = InteractiveSynonymCLI(config_path)
        cli.run()

    except FileNotFoundError as e:
        print(f"❌ Configuration error: {e}")
        return 1
    except Exception as e:
        print(f"❌ Failed to initialize: {e}")
        logger.exception("Initialization failed")
        return 1
    finally:
        if "cli" in locals():
            cli.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
