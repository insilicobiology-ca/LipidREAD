from src.models.lipid_components import *
from src.models.models import ChainPosition
from typing import Set


class LipidParserStrategy(ABC):
    @abstractmethod
    def parse(self, lipid_name: str) -> LipidComponents:
        pass




class SphingoLipidParserStrategy(LipidParserStrategy):
    def __init__(self, known_headgroups: Set[str]):
        self.known_headgroups = known_headgroups

    def parse(self, lipid_name: str) -> SphingoLipidComponents:
        if lipid_name in self.known_headgroups:
            components = SphingoLipidComponents(lipid_name, "", "", None,
                                          is_headgroup_only=True, is_backbone_only=False)
            components.level = LipidLevel.CLASS
            return components

        try:
            # First try to match species-level notation
            species_match = re.match(r"""
                            ([^(]+)                # Headgroup
                            \(
                            (\d+:\d+(?:;\d+)?)    # Species composition (e.g., 36:1;2)
                            \)
                        """, lipid_name, re.VERBOSE)

            if species_match:
                headgroup, species_comp = species_match.groups()
                components = SphingoLipidComponents(headgroup, "", "", None)
                components.level = LipidLevel.SPECIES
                components.species_composition = species_comp
                return components

            # If not species-level, try regular parsing with updated regex to include deuteration
            match = re.match(r"""
                ([^(]+)                # Headgroup
                \(
                ([\w:]+(?:-D\d{1,2})?) # Backbone with optional deuteration
                (?:
                    /
                    ([\w:]+(?:-D\d{1,2})?) # N-acyl chain with optional deuteration
                )?
                (?:\((OH)\))?          # Optional hydroxylation
                \)
            """, lipid_name, re.VERBOSE)
            if not match:
                raise ValueError(f"Invalid Sphingolipid format: {lipid_name}")
            headgroup, backbone, n_acyl_chain, hydroxy = match.groups()
            if n_acyl_chain is None:
                components = SphingoLipidComponents(headgroup, backbone, n_acyl_chain, hydroxy, is_backbone_only=True)
                components.level = LipidLevel.SPECIES
                components.species_composition = backbone
            else:
                components = SphingoLipidComponents(headgroup, backbone, n_acyl_chain, hydroxy)
                components.level = LipidLevel.SN_POSITION

            return components
        except Exception:
            raise ValueError(f"Invalid Sphingolipid format: {lipid_name}")


class SphingoidBaseParserStrategy(LipidParserStrategy):
    def __init__(self, known_headgroups: Set[str]):
        self.known_headgroups = known_headgroups

    def parse(self, lipid_name: str) -> SphingoidBaseComponents:
        if lipid_name in self.known_headgroups:
            return SphingoidBaseComponents(lipid_name, "", is_headgroup_only=True)
        try:
            # Updated regex to include deuteration
            match = re.match(r"""
                ([^(]+)                # Headgroup
                \(
                ([\w:]+(?:-D\d{1,2})?) # Backbone with optional deuteration
                \)
            """, lipid_name, re.VERBOSE)
            if not match:
                raise ValueError(f"Invalid Sphingoid Base format: {lipid_name}")
            headgroup, backbone = match.groups()
            return SphingoidBaseComponents(headgroup, backbone)
        except Exception:
            raise ValueError(f"Invalid Sphingoid Base format: {lipid_name}")


class GlycerophosphoLipidParserStrategy(LipidParserStrategy):
    LYSO_HEADGROUPS = {
        'PC': 'LPC',
        'PE': 'LPE',
        'PS': 'LPS',
        'PI': 'LPI',
        'PG': 'LPG',
        'PA': 'LPA'
    }

    def __init__(self, known_headgroups: Set[str]):
        self.known_headgroups = known_headgroups

    def parse(self, lipid_name: str) -> GlycerophosphoLipidComponents:
        if lipid_name in self.known_headgroups:
            components = GlycerophosphoLipidComponents(lipid_name, "", "", None,
                                                       is_headgroup_only=True)
            components.level = LipidLevel.CLASS
            return components
        try:
            # First try to match species-level notation
            species_match = re.match(r"""
                            (\w+)                      # Headgroup
                            \(
                            ((?:O-|P-)?               # Optional O- or P- modification
                            \d+:\d+(?:;\d+)?          # Species composition (e.g., 36:1)
                            (?:_\d+:\d+(?:;\d+)?)?)   # Optional second part for molecular species
                            \)
                        """, lipid_name, re.VERBOSE)

            if species_match:
                headgroup, species_comp = species_match.groups()
                # Extract modification if present
                if species_comp.startswith(('O-', 'P-')):
                    modification = species_comp[:2]
                    headgroup = f"{headgroup}({modification})"
                components = GlycerophosphoLipidComponents(headgroup, "", "", None)
                components.species_composition = species_comp
                components.level = LipidLevel.MOLECULAR_SPECIES if '_' in species_comp else LipidLevel.SPECIES
                return components

            # Updated regex to include deuteration
            match = re.match(r"""
                (\w+)                  # Headgroup
                \(
                ((?:O-|P-)?[\w:]+(?:-D\d{1,2})?) # sn1 chain with optional deuteration
                /
                ([\w:]+(?:-D\d{1,2})?) # sn2 chain with optional deuteration
                (?:\((OH)\))?          # Optional hydroxylation
                \)
            """, lipid_name, re.VERBOSE)
            if not match:
                raise ValueError(f"Invalid Glycerophospholipid format: {lipid_name}")

            headgroup, sn1_chain, sn2_chain, hydroxy = match.groups()
            base_headgroup = headgroup
            modification = None

            # Check for O- or P- modification in sn1
            if sn1_chain.startswith(('O-', 'P-')):
                modification = sn1_chain[:2]
                # sn1_chain = sn1_chain[2:]  # Remove modification from chain
                # headgroup = f"{headgroup}({modification})"

            # Check if this is actually a lyso lipid
            if sn1_chain == "0:0" or sn2_chain == "0:0":
                if base_headgroup in self.LYSO_HEADGROUPS:
                    lyso_headgroup = self.LYSO_HEADGROUPS[base_headgroup]
                    if sn1_chain == "0:0":
                        active_chain = sn2_chain
                        active_position = ChainPosition.SN2
                    else:  # sn2_chain == "0:0"
                        active_chain = sn1_chain
                        active_position = ChainPosition.SN1

                    components = LysoGlycerophospholipidComponents(
                        headgroup=lyso_headgroup,
                        active_chain=active_chain,
                        active_chain_position=active_position,
                        modification=modification
                    )
                    components.level = LipidLevel.SN_POSITION
                    return components

            # If not lyso, return regular glycerophospholipid components
            headgroup = f"{base_headgroup}({modification})" if modification else base_headgroup
            components = GlycerophosphoLipidComponents(headgroup, sn1_chain, sn2_chain, hydroxy)
            components.level = LipidLevel.SN_POSITION
            return components
        except Exception:
            raise ValueError(f"Invalid Glycerophospholipid format: {lipid_name}")


class LysoGlycerophosphoLipidParserStrategy(LipidParserStrategy):
    def __init__(self, known_headgroups: Set[str]):
        self.known_headgroups = known_headgroups

    def parse(self, lipid_name: str) -> LysoGlycerophospholipidComponents:
        if lipid_name in self.known_headgroups:
            components = LysoGlycerophospholipidComponents(lipid_name, "", ChainPosition.SN1,
                                                           is_headgroup_only=True)
            components.level = LipidLevel.CLASS
            return components
        try:
            # First try to match species-level notation (no sn position specified)
            # Updated to handle O- and P- modifications in species-level notation
            species_match = re.match(r"""
                           (\w+)                        # Headgroup
                           \(
                           ((?:O-|P-)?                  # Optional O- or P- modification
                           \d+:\d+(?:;\d+)?)            # Chain composition
                           \)
                       """, lipid_name, re.VERBOSE)

            if species_match:
                headgroup, species_comp = species_match.groups()

                # Extract modification if present
                if species_comp.startswith(('O-', 'P-')):
                    modification = species_comp[:2]
                    headgroup = f"{headgroup}({modification})"
                components = LysoGlycerophospholipidComponents(
                    headgroup=headgroup,
                    active_chain=species_comp,
                    active_chain_position=ChainPosition.SN1)
                components.species_composition = species_comp
                components.level = LipidLevel.SPECIES
                return components

            # If not species-level, try regular parsing
            match = re.match(r"""
                (\w+)                  # Headgroup
                \(
                ((?:O-|P-)?[\w:]+(?:-D\d{1,2})?) # sn1 chain with optional deuteration
                /
                ([\w:]+(?:-D\d{1,2})?) # sn2 chain with optional deuteration
                \)
            """, lipid_name, re.VERBOSE)
            if not match:
                raise ValueError(f"Invalid LysoGlycerophospholipid format: {lipid_name}")
            headgroup, sn1_chain, sn2_chain = match.groups()
            modification = None

            # Extract modification if present
            if sn1_chain.startswith(('O-', 'P-')):
                active_chain = sn1_chain
                active_position = ChainPosition.SN1
                modification = sn1_chain[:2]
                lyso_headgroup = f"{headgroup}({sn1_chain[:2]})"

                # sn1_chain = sn1_chain[2:]

            elif sn1_chain == "0:0":
                lyso_headgroup = f"{headgroup}_sn2"
                active_chain = sn2_chain
                active_position = ChainPosition.SN2
            elif sn2_chain == "0:0":
                lyso_headgroup = f"{headgroup}_sn1"
                active_chain = sn1_chain
                active_position = ChainPosition.SN1
            else:
                raise ValueError(f"Invalid LysoGlycerophospholipid: neither chain is 0:0")

            components = LysoGlycerophospholipidComponents(
                headgroup=lyso_headgroup,
                active_chain=active_chain,
                active_chain_position=active_position,
                modification=modification
            )
            components.level = LipidLevel.SN_POSITION
            return components
        except Exception:
            raise ValueError(f"Invalid LysoGlycerophospholipid format: {lipid_name}")


class AcylCeramideParserStrategy(LipidParserStrategy):
    def __init__(self, known_headgroups: Set[str]):
        self.known_headgroups = known_headgroups

    def parse(self, lipid_name: str) -> AcylCeramideComponents:
        if lipid_name in self.known_headgroups:
            return AcylCeramideComponents(lipid_name, "", "", "", None,
                                          is_headgroup_only=True)
        try:
            # Updated regex to include deuteration
            match = re.match(r"""
                1-O-
                ([\w:]+(?:-D\d{1,2})?) # O-acyl chain with optional deuteration
                -Cer
                \(
                ([\w:]+(?:-D\d{1,2})?) # Backbone with optional deuteration
                /
                ([\w:]+(?:-D\d{1,2})?) # N-acyl chain with optional deuteration
                (?:\((OH)\))?          # Optional hydroxylation
                \)
            """, lipid_name, re.VERBOSE)
            if not match:
                raise ValueError(f"Invalid Acyl Ceramide format: {lipid_name}")
            o_acyl_chain, backbone, n_acyl_chain, hydroxy = match.groups()
            return AcylCeramideComponents("1-O-acyl-Cer", o_acyl_chain, backbone, n_acyl_chain, hydroxy)
        except Exception:
            raise ValueError(f"Invalid Acyl Ceramide format: {lipid_name}")


class TriAcylsParserStrategy(LipidParserStrategy):
    def __init__(self, known_headgroups: Set[str]):
        self.known_headgroups = known_headgroups

    def parse(self, lipid_name: str) -> TriAcylComponents:
        if lipid_name in self.known_headgroups:
            components = TriAcylComponents(lipid_name, "", "", "", hydroxy="",
                                           is_headgroup_only=True)
            components.level = LipidLevel.CLASS
            return components
        try:
            # First try to match species-level notation
            species_match = re.match(r"""
                            ([\w-]+)                   # Headgroup
                            \(
                            ((?:O-|P-)?               # Optional O- or P- modification
                            \d+:\d+(?:;\d+)?          # Species composition (e.g., 54:3)
                            (?:_\d+:\d+(?:;\d+)?      # Optional second part for molecular species
                            (?:_\d+:\d+(?:;\d+)?)?)?) # Optional third part for molecular species
                            \)
                        """, lipid_name, re.VERBOSE)

            if species_match:
                headgroup, species_comp = species_match.groups()
                # Extract modification if present
                if species_comp.startswith(('O-', 'P-')):
                    modification = species_comp[:2]
                    headgroup = f"{headgroup}({modification})"
                components = TriAcylComponents(headgroup, "", "", "", "")
                components.species_composition = species_comp
                components.level = LipidLevel.MOLECULAR_SPECIES if '_' in species_comp else LipidLevel.SPECIES
                return components

            # If not species-level, try regular parsing
            match = re.match(r"""
                ([\w-]+)               # Headgroup
                \(
                ((?:O-|P-)?[\w:]+(?:-D\d{1,2})?) # sn1 chain with optional modification and deuteration
                /
                ([\w:]+(?:-D\d{1,2})?) # sn2 chain with optional deuteration
                (?:
                    /
                    ([\w:]+(?:-D\d{1,2})?) # Optional sn3 chain with optional deuteration
                )?
                \)
            """, lipid_name, re.VERBOSE)
            if not match:
                raise ValueError(f"Invalid TriAcyl format: {lipid_name}")
            headgroup, sn1_chain, sn2_chain, sn3_chain = match.groups()
            # Set sn3_chain to "0:0" if it's None
            sn3_chain = sn3_chain or "0:0"
            components = TriAcylComponents(headgroup, sn1_chain, sn2_chain, sn3_chain, "")
            components.level = LipidLevel.SN_POSITION
            return components
        except Exception:
            raise ValueError(f"Invalid TriAcyl format: {lipid_name}")


class MonoAcylParserStrategy(LipidParserStrategy):
    def __init__(self, known_headgroups: Set[str]):
        self.known_headgroups = known_headgroups

    def parse(self, lipid_name: str) -> MonoAcylComponents:
        # Step 1: Handle class-level names like "MG", "MG_sn1", "MG(O-)"
        if lipid_name in self.known_headgroups:
            # The full name (e.g., "MG_sn1") becomes the headgroup, which is correct for matching.
            return MonoAcylComponents(lipid_name, "", is_headgroup_only=True)

        try:
            # Step 2: Parse fully specified names using a robust regex.
            # This regex captures the base headgroup, optional modification, optional position, and full chain part.
            # e.g., "MG(O-16:0/0:0/0:0)" or "MG(18:1)" or "MG_sn2(18:1)"
            match = re.match(r"(\w+)(\(O-\))?(_sn[123])?\((.*)\)", lipid_name)
            if not match:
                raise ValueError(f"Initial format match failed for: {lipid_name}")

            base_headgroup, modification, position, full_chain_part = match.groups()

            # Clean up the position string (e.g., from "_sn1" to "sn1")
            if position:
                position = position.strip('_')

            # --- Logic to find the single active chain ---
            chains = [chain.strip() for chain in full_chain_part.split('/')]
            non_zero_chains = [chain for chain in chains if chain != "0:0"]

            acyl_chain = ""
            if len(non_zero_chains) == 1:
                acyl_chain = non_zero_chains[0]
            elif len(chains) == 1 and chains[0] != "0:0":
                acyl_chain = chains[0]
            else:
                raise ValueError(f"Expected exactly one non-zero acyl chain, found {len(non_zero_chains)}")

            # Reconstruct the full headgroup name for consistency (e.g., "MG(O-)" or "MG_sn2")
            # This is the value that will be used by the matcher.
            final_headgroup = base_headgroup
            if modification:
                final_headgroup += modification
            if position:
                final_headgroup += f"_{position}"

            return MonoAcylComponents(
                headgroup=final_headgroup,
                acyl_chain=acyl_chain,
                position=position,
                modification=modification,
                is_headgroup_only=False
            )

        except Exception as e:
            raise ValueError(f"Invalid MonoAcyl lipid format for '{lipid_name}'") from e

class CardiolipinParserStrategy(LipidParserStrategy):
    def __init__(self, known_headgroups: Set[str]):
        self.known_headgroups = known_headgroups

    def parse(self, lipid_name: str) -> CardiolipinComponents:
        if lipid_name in self.known_headgroups:
            return CardiolipinComponents(lipid_name, "", "", "", "",
                                         is_headgroup_only=True)
        try:
            # Updated regex to include deuteration
            match = re.match(r"""
                            ([\w-]+)               # Headgroup
                            \(
                            ([\w:]+(?:-D\d{1,2})?) # sn1 chain with optional deuteration
                            /
                            ([\w:]+(?:-D\d{1,2})?) # sn2 chain with optional deuteration
                            /
                            ([\w:]+(?:-D\d{1,2})?) # sn1-prime chain with optional deuteration
                            /
                            ([\w:]+(?:-D\d{1,2})?) # sn2-prime chain with optional deuteration
                            \)
                        """, lipid_name, re.VERBOSE)
            if not match:
                raise ValueError(f"Invalid Cardiolipin format: {lipid_name}")
            headgroup, sn1_chain, sn2_chain, sn1_prime_chain, sn2_prime_chain = match.groups()
            return CardiolipinComponents(headgroup, sn1_chain, sn2_chain, sn1_prime_chain, sn2_prime_chain)
        except Exception:
            raise ValueError(f"Invalid Cardiolipin format: {lipid_name}")
