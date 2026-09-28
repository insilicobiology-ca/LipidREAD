from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Dict, Callable, Optional
from enum import Enum
from dataclasses import dataclass, field
from src.models.lipid_level import LipidLevel
from src.models.models import ChainPosition
import re


@dataclass
class LipidComponents(ABC):
    headgroup: str
    is_headgroup_only: bool = False
    level: LipidLevel = LipidLevel.UNDEFINED
    species_composition: Optional[str] = None

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def get_chains(self) -> List[str]:
        pass


class SphingoLipidComponents(LipidComponents):
    def __init__(self, headgroup: str, backbone: str, n_acyl_chain: Optional[str], hydroxy: Optional[str] = None,
                 is_headgroup_only: bool = False, is_backbone_only: bool = False):
        super().__init__(headgroup, is_headgroup_only)
        self.backbone = backbone
        self.n_acyl_chain = n_acyl_chain
        self.hydroxy = hydroxy
        self.is_backbone_only = is_backbone_only

    def __str__(self):
        if self.is_headgroup_only:
            return self.headgroup
        hydroxy_str = f"({self.hydroxy})" if self.hydroxy else ""
        return f"{self.headgroup}({self.backbone}/{self.n_acyl_chain}{hydroxy_str})"

    def get_chains(self) -> List[str]:
        return [self.backbone, self.n_acyl_chain]

    def get_backbone_family(self) -> str:
        return f"{self.headgroup}({self.backbone})"


class SphingoidBaseComponents(LipidComponents):
    def __init__(self, headgroup: str, backbone: str, is_headgroup_only: bool = False):
        super().__init__(headgroup, is_headgroup_only)
        self.backbone = backbone

    def __str__(self):
        if self.is_headgroup_only:
            return self.headgroup
        return f"{self.headgroup}({self.backbone})"

    def get_chains(self) -> List[str]:
        return [self.backbone]


class GlycerophosphoLipidComponents(LipidComponents):
    def __init__(self, headgroup: str, sn1_chain: str, sn2_chain: str, hydroxy: Optional[str] = None,
                 is_headgroup_only: bool = False):
        super().__init__(headgroup, is_headgroup_only)
        self.sn1_chain = sn1_chain
        self.sn2_chain = sn2_chain
        self.hydroxy = hydroxy
        # NEW ATTRIBUTE
        self.is_template = 'X:Y' in sn1_chain or 'X:Y' in sn2_chain

    def __str__(self):
        if self.is_headgroup_only:
            return self.headgroup
        hydroxy_str = f"({self.hydroxy})" if self.hydroxy else ""
        return f"{self.headgroup}({self.sn1_chain}/{self.sn2_chain}{hydroxy_str})"

    def get_chains(self) -> List[str]:
        return [self.sn1_chain, self.sn2_chain]


class LysoGlycerophospholipidComponents(LipidComponents):
    def __init__(self, headgroup: str, active_chain: str, active_chain_position: ChainPosition,
                 modification: Optional[str] = None, is_headgroup_only: bool = False):
        super().__init__(headgroup, is_headgroup_only)
        self.active_chain = active_chain
        self.modification = modification
        # For modified lipids (O- or P-), force SN1 position
        self.active_chain_position = ChainPosition.SN1 if modification else active_chain_position

    def __str__(self):
        if self.is_headgroup_only:
            return self.headgroup
        if '_' in self.headgroup:  # if lyso lipid has _sn1 or _sn2 info, get rid of them make str representation
            base_headgroup = self.headgroup.split('_')[0]
        else:
            base_headgroup = self.headgroup.split('(')[0]  # otherwise lyso lipid has (O-) or (P-)

        # For modified lipids, the active chain must be at sn1
        if self.modification:
            return f"{base_headgroup}({self.modification}{self.active_chain[2:]}/0:0)"

        # Place 0:0 and active chain in correct positions based on active_chain_position
        sn1_chain = self.active_chain if self.active_chain_position == ChainPosition.SN1 else '0:0'
        sn2_chain = self.active_chain if self.active_chain_position == ChainPosition.SN2 else '0:0'

        return f"{base_headgroup}({sn1_chain}/{sn2_chain})"

    def get_chains(self) -> List[str]:
        """Return chains in sn1/sn2 order."""
        if self.active_chain_position == ChainPosition.SN1:
            return [self.active_chain, '0:0']
        return ['0:0', self.active_chain]

    @property
    def present_chain(self):
        """Return the position of the active chain."""
        return self.active_chain_position


class CardiolipinComponents(LipidComponents):
    def __init__(self, headgroup: str, sn1_chain: str, sn2_chain: str, sn1_prime_chain: str, sn2_prime_chain: str,
                 hydroxy: Optional[str] = None, is_headgroup_only: bool = False):
        super().__init__(headgroup, is_headgroup_only)
        self.sn1_chain = sn1_chain
        self.sn2_chain = sn2_chain
        self.sn1_prime_chain = sn1_prime_chain
        self.sn2_prime_chain = sn2_prime_chain
        self.hydroxy = hydroxy

    def __str__(self):
        if self.is_headgroup_only:
            return self.headgroup
        hydroxy_str = f"({self.hydroxy})" if self.hydroxy else ""
        return f"{self.headgroup}({self.sn1_chain}/{self.sn2_chain}/{self.sn1_prime_chain}/{self.sn2_prime_chain}{hydroxy_str})"

    def get_chains(self) -> List[str]:
        return [self.sn1_chain, self.sn2_chain, self.sn1_prime_chain, self.sn2_prime_chain]


class AcylCeramideComponents(SphingoLipidComponents):
    def __init__(self, headgroup: str, o_acyl_chain: str, backbone: str, n_acyl_chain: str,
                 hydroxy: Optional[str] = None, is_headgroup_only: bool = False):
        super().__init__(headgroup, backbone, n_acyl_chain, hydroxy, is_headgroup_only)
        self.o_acyl_chain = o_acyl_chain

    def __str__(self):
        if self.is_headgroup_only:
            return self.headgroup
        hydroxy_str = f"({self.hydroxy})" if self.hydroxy else ""
        return f"1-O-{self.o_acyl_chain}-Cer({self.backbone}/{self.n_acyl_chain}{hydroxy_str})"

    def get_chains(self) -> List[str]:
        return [self.o_acyl_chain, self.backbone, self.n_acyl_chain]


class TriAcylComponents(LipidComponents):
    def __init__(self, headgroup: str, sn1_chain: str, sn2_chain: str, sn3_chain: str, hydroxy: str,
                 is_headgroup_only: bool = False):
        super().__init__(headgroup, is_headgroup_only)
        self.sn1_chain = sn1_chain
        self.sn2_chain = sn2_chain
        self.sn3_chain = sn3_chain
        self.hydroxy = hydroxy

    def __str__(self):
        if self.is_headgroup_only:
            return self.headgroup
        return f"{self.headgroup}({self.sn1_chain}/{self.sn2_chain}/{self.sn3_chain})"

    def get_chains(self) -> List[str]:
        return [self.sn1_chain, self.sn2_chain, self.sn3_chain]


class MonoAcylComponents(LipidComponents):
    def __init__(self, headgroup: str, acyl_chain: str,
                 position: Optional[str] = None,  # NEW: e.g., "sn1", "sn2"
                 modification: Optional[str] = None,  # NEW: e.g., "O-"
                 is_headgroup_only: bool = False):
        super().__init__(headgroup, is_headgroup_only)
        self.acyl_chain = acyl_chain
        self.position = position
        self.modification = modification

    def __str__(self):
        """
        Generates a string representation of the mono-acyl lipid.

        MODIFIED: Now explicitly shows "0:0" for empty positions if the
        acyl chain position is known.
        """
        if self.is_headgroup_only:
            return self.headgroup

        # Get the base headgroup, stripping any potential _snX suffix or (O-) modification
        base_headgroup = self.headgroup.split('_')[0].split('(')[0]
        mod_str = self.modification or ""
        full_acyl_chain = f"{mod_str}{self.acyl_chain}"

        # If position is known, construct the full three-part representation
        if self.position == "sn1":
            return f"{base_headgroup}({full_acyl_chain}/0:0/0:0)"
        elif self.position == "sn2":
            return f"{base_headgroup}(0:0/{full_acyl_chain}/0:0)"
        elif self.position == "sn3":
            return f"{base_headgroup}(0:0/0:0/{full_acyl_chain})"
        else:
            # If position is unknown or not applicable, fall back to the simple representation.
            # This is important for cases like "CE(18:1)" where there is only one position.
            return f"{base_headgroup}({full_acyl_chain})"

    def get_chains(self) -> List[str]:
        return [self.acyl_chain]

    def get_chains(self) -> List[str]:
        return [self.acyl_chain]
