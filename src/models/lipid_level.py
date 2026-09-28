# In src/models/lipid_level.py
from enum import Enum
import re


class LipidLevel(Enum):
    UNDEFINED = 0
    CATEGORY = 1
    CLASS = 2
    SPECIES = 3
    MOLECULAR_SPECIES = 4
    SN_POSITION = 5
    STRUCTURE_DEFINED = 6

    @classmethod
    def determine_level(cls, lipid_name: str, components: 'LipidComponents') -> 'LipidLevel':
        """
        Determine the structural level of a lipid based on its name and components.
        """
        if components.is_headgroup_only:
            return cls.CLASS

        # Check for species level notation (sum composition)
        if '(' in lipid_name and ')' in lipid_name:
            inner_content = lipid_name[lipid_name.index('(') + 1:lipid_name.index(')')]

            # Species level check (e.g., PC(36:0))
            if re.match(r'^\d+:\d+(?:;\d+)?$', inner_content):
                return cls.SPECIES

            # Check for molecular species level (chain separation with '_')
            if '_' in inner_content:
                return cls.MOLECULAR_SPECIES

            # Check for sn-position level (chain separation with '/')
            if '/' in inner_content:
                return cls.SN_POSITION

            # Check for structure defined level (contains specific modifications)
            if ';O' in inner_content or ';P' in inner_content or '(OH)' in inner_content:
                return cls.STRUCTURE_DEFINED

        return cls.UNDEFINED