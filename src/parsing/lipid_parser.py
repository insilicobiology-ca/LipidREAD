from src.data.config_loader import ConfigLoader
from src.parsing.lipid_parser_strategies import *
from src.models.lipid_components import LipidComponents
from typing import Optional, MutableMapping


class LipidType(Enum):
    SPHINGOLIPID = 1
    GLYCEROPHOSPHOLIPID = 2
    LYSOGLYCEROPHOSPHOLIPID = 3
    ACYL_CERAMIDE = 4
    TRIACYLS = 5
    SPHINGOID_BASE = 6
    MONOACYL = 7
    CARDIOLIPIN = 9

class LipidLevel(Enum):
    # Adapted from Goslin project: https://github.com/lifs-tools/pygoslin
    # Undefined / non-inferable lipid level
    UNDEFINED = 0

    # Sphingolipids, Glycerophospholipids, Glycerolipids
    CATEGORY = 1

    # Glyerophospholipids -> Glycerophosphocholine (PC)
    CLASS = 2

    # Phosphatidylcholine (36:0) / PC(36:0)
    SPECIES = 3

    # Phosphatidylinositol (8:0-8:0) or PC(8:0-8:0)
    MOLECULAR_SPECIES = 4

    # Phosphatidylinositol (8:0;O2/8:0) or PI(8:0;O2/8:0)
    SN_POSITION = 5

    # Phosphatidylinositol (8:0;(OH)2/8:0) or PI(8:0;(OH)2/8:0)
    STRUCTURE_DEFINED = 6


class LipidParserFactory:
    def __init__(self, config_loader: ConfigLoader):
        self.lipid_prefixes = config_loader.get_lipid_prefixes()
        self.known_headgroups = set(config_loader.get_known_headgroups())

        self.strategies = {
            LipidType.SPHINGOLIPID: SphingoLipidParserStrategy(self.known_headgroups),
            LipidType.LYSOGLYCEROPHOSPHOLIPID: LysoGlycerophosphoLipidParserStrategy(self.known_headgroups),
            LipidType.GLYCEROPHOSPHOLIPID: GlycerophosphoLipidParserStrategy(self.known_headgroups),
            LipidType.ACYL_CERAMIDE: AcylCeramideParserStrategy(self.known_headgroups),
            LipidType.TRIACYLS: TriAcylsParserStrategy(self.known_headgroups),
            LipidType.SPHINGOID_BASE: SphingoidBaseParserStrategy(self.known_headgroups),
            LipidType.MONOACYL: MonoAcylParserStrategy(self.known_headgroups),
            LipidType.CARDIOLIPIN: CardiolipinParserStrategy(self.known_headgroups)
        }

    def get_parser(self, lipid_name: str) -> LipidParserStrategy:
        for lipid_type, prefixes in self.lipid_prefixes.items():
            if any(lipid_name.startswith(prefix) for prefix in prefixes):
                return self.strategies[LipidType[lipid_type]]
        raise ValueError(f"Unknown lipid type: {lipid_name}")

class LipidComponentCache:
    """
        Centralized cache for parsed lipid components with thread-safe implementation.
    """
    def __init__(self):
        self._cache: MutableMapping[str, LipidComponents] = {}

    def get(self, lipid: str) -> Optional[LipidComponents]:
        """Get cached lipid components if they exist."""
        return self._cache.get(lipid)

    def set(self, lipid: str, components: LipidComponents) -> None:
        """Cache parsed lipid components."""
        self._cache[lipid] = components

    def clear(self) -> None:
        """Clear the cache."""
        self._cache.clear()


class LipidParser:
    def __init__(self, factory: LipidParserFactory, components_cache: LipidComponentCache):
        self.factory = factory
        self.components_cache = components_cache

    def parse_lipid(self, lipid_name: str) -> LipidComponents:
        """
        Parse lipid with caching. Returns cached result if available.
        """
        # Check cache first
        cached_components = self.components_cache.get(lipid_name)
        if cached_components is not None:
            return cached_components

        # Parse if not in cache
        parser = self.factory.get_parser(lipid_name)
        components = parser.parse(lipid_name)

        # Cache the result
        self.components_cache.set(lipid_name, components)
        return components



