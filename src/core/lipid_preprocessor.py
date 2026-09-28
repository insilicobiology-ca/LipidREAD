from dataclasses import dataclass
from typing import Dict, Set, List, Tuple
import logging
from src.parsing.lipid_parser import LipidParser
from src.models.lipid_components import LipidComponents

logger = logging.getLogger(__name__)

@dataclass
class PreprocessingResult:
    """Contains the results of lipid preprocessing."""
    headgroup_map: Dict[str, Set[str]]  # Maps headgroup to set of lipids
    unparseable_lipids: List[str]
    all_valid_lipids: Set[str]  # Union of all parseable lipids

class LipidPreprocessor:
    """Service class responsible for preprocessing lipids."""

    def __init__(self, lipid_parser: LipidParser):
        self.lipid_parser = lipid_parser

    def preprocess_lipids(self, lipids: List[str]) -> PreprocessingResult:
        """
        Preprocess a list of lipids, categorizing them by headgroup.
        :param lipids (List[str]): List of lipid names to preprocess.
        :return PreprocessingResult: Contains categorized lipids and processing metadata.:
        """

        headgroup_map: Dict[str, Set[str]] = {}
        unparseable_lipids: List[str] = []
        all_valid_lipids: Set[str] = set()

        for lipid in lipids:
            try:
                components = self.lipid_parser.parse_lipid(lipid)
                headgroup = components.headgroup

                if headgroup not in headgroup_map:
                    headgroup_map[headgroup] = set()

                headgroup_map[headgroup].add(lipid)
                all_valid_lipids.add(lipid)

            except ValueError as error:
                logger.warning(f"Could not parse lipid {lipid}: {str(error)}")
                unparseable_lipids.append(lipid)
            except Exception as error:
                logger.error(f"Unexpected error processing lipid {lipid}: {str(error)}")
                unparseable_lipids.append(lipid)

        return PreprocessingResult(
            headgroup_map=headgroup_map,
            unparseable_lipids=unparseable_lipids,
            all_valid_lipids=all_valid_lipids
        )