import yaml
from typing import List, Dict, Any, Tuple


class ConfigLoader:
    def __init__(self, config_file: str):
        with open(config_file, 'r') as file:
            self.config = yaml.safe_load(file)

    def get_lipid_prefixes(self) -> Dict[str, List[str]]:
        return self.config['lipid_prefixes']

    def get_lipid_categories(self) -> Dict[str, Dict[str, List[str]]]:
        """
        Returns the complete lipid_categories dictionary
        """
        return self.config.get('lipid_categories', {})

    def get_category_for_lipid_class(self, lipid_class: str) -> str:
        """
        Find the category for a given lipid class.

        Args:
            lipid_class: The lipid class to find the category for

        Returns:
            Category name as string, or "Other" if not found
        """
        categories = self.get_lipid_categories()

        for category, details in categories.items():
            if lipid_class in details['classes']:
                return category
        return 'Other'

    def get_matching_strategies(self) -> Dict[str, List[str]]:
        return self.config['matching_strategies']

    def get_super_matching_strategies(self) -> Dict[str, List[str]]:
        return self.config['super_matching_strategies']

    def get_db_name(self) -> str:
        return self.config['database']['name']

    def get_db_user(self) -> str:
        return self.config['database']['user']

    def get_db_password(self) -> str:
        return self.config['database']['password']

    def get_db_host(self) -> str:
        return self.config['database']['host']

    def get_db_port(self) -> int:
        return self.config['database']['port']

    def get_db_schema(self) -> str:
        return self.config['database'].get('schema', 'public')

    def get_super_reaction_ids(self) -> List[int]:
        return self.config['super_reaction_ids']

    def get_output_paths(self) -> Dict[str, str]:
        return self.config['output_paths']

    def get_known_headgroups(self) -> List[str]:
        return self.config['known_headgroups']

    def get_family_circle_params(self) -> Dict[str, Tuple[float, float, float, float]]:
        return self.config.get('family_circle_params', {})

    def get_enzyme_coordinates(self) -> Dict[str, Tuple[float, float]]:
        return self.config.get('enzyme_coordinates', {})

    def get_family_molecule_ids(self) -> Dict[str, int]:
        return self.config.get('family_molecule_ids', {})

    def get_updater_setting(self, key: str, default: Any = None) -> Any:
        """Generic getter for database_updater settings."""
        return self.config.get('database_updater', {}).get(key, default)

    def get_staging_lipids_table(self) -> str:
        return self.get_updater_setting('staging_lipids_table', 'lipograph.staging_swisslipids_lipids')

    def get_staging_enzymes_table(self) -> str:
        return self.get_updater_setting('staging_enzymes_table', 'lipograph.staging_swisslipids_enzymes')

    def get_staging_lipids2uniprot_table(self) -> str:
        return self.get_updater_setting('staging_lipids2uniprot_table', 'lipograph.staging_swisslipids_lipids2uniprot')

    def get_update_output_dir(self) -> str:
        return self.get_updater_setting('update_output_dir', 'db_updates')

    def get_backup_dir(self) -> str:
        return self.get_updater_setting('backup_dir', 'db_backups')

    def get_log_dir(self) -> str:
        return self.get_updater_setting('log_dir', 'logs')

    def get_review_notification_emails(self) -> List[str]:
        emails_str = self.get_updater_setting('review_notification_email', '')
        if not emails_str:
            return []
        return [email.strip() for email in emails_str.split(',')]

    def get_filtered_molecule_ids_for_pairs(self) -> List[int]:
        return self.get_updater_setting('filter_reaction_pairs_for_molecules', [])

    def get_ignored_metabolites_for_pairs(self) -> List[str]:
        return self.get_updater_setting('ignore_metabolites_in_pairs', [])

    def get(self, key: str, default: Any = None) -> Any:
        """
        Generic method to get any configuration value
        """
        return self.config.get(key, default)