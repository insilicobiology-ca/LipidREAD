from typing import Set, Tuple
import math
import colorsys
from src.models.network_models import NetworkNode, NetworkEdge
from src.parsing.lipid_parser import LipidParser, LipidParserFactory
from src.models.models import Reaction
from src.data.config_loader import ConfigLoader
from src.models.lipid_components import *
import logging
import random

@dataclass
class CircleParameters:
    center_y: float
    center_z: float
    radius: float

class LipidNetworkManager:
    def __init__(self, json_data: Dict, reaction_dict: Dict, super_reaction_dict: Dict, config_loader: ConfigLoader,
                 lipid_parser: LipidParser):
        self.json_data = json_data
        self.reaction_dict = reaction_dict
        self.super_reaction_dict = super_reaction_dict
        self.config_loader = config_loader
        self.lipid_parser = LipidParser(LipidParserFactory(config_loader))
        # self.node_placer = NodePlacer()
        self.color_manager = ColorManager()
        self.lipid_families: Dict[str, Set[str]] = {}
        self.family_circles: Dict[str, CircleParameters] = {}
        self.nodes: Dict[str, NetworkNode] = {}
        self.edges: List[NetworkEdge] = []
        self.enzyme_dict: Dict[Tuple[Tuple[str, ...], str], Set[str]] = {}
        # self.node_placer = NodePlacer()
        self.lipid_components: Dict[str, LipidComponents] = {}
        self.family_circle_params = config_loader.get_family_circle_params()
        self.enzyme_coordinates = config_loader.get_enzyme_coordinates()
        self.node_placer = NodePlacer(self.family_circle_params)
        self.enzyme_position_manager = EnzymePositionManager(self.nodes)
        self.processed_reaction_edges = set()

    def process_network(self):
        self._process_reactions(self.reaction_dict, is_super=False)
        self._process_reactions(self.super_reaction_dict, is_super=True)
        self._place_nodes()
        self._create_edges()
        self._update_json()

    def _process_reactions(self, reaction_dict: Dict, is_super: bool):
        for key, value in reaction_dict.items():
            reaction = value.reaction
            if is_super:
                reactant1, reactant2, product = key[:3]
                reactants = [reactant1, reactant2]
            else:
                reactant, product = key[:2]
                reactants = [reactant]

            for r in reactants:
                self._process_lipid(r)
            self._process_lipid(product)

            self._update_enzyme_dict(reactants, product, reaction)

    def _normalize_family_name(self, family: str) -> str:
        # Remove _sn1, _sn2, (O-), and (P-) suffixes
        normalized = re.sub(r'_(sn1|sn2)|\(O-\)|\(P-\)$', '', family)

        # List of lysoglycerophospholipid prefixes
        lyso_prefixes = ['LPC', 'LPS', 'LPG', 'LPE', 'LPA', 'LPI']

        # If the normalized family starts with any of these prefixes, return just the prefix
        for prefix in lyso_prefixes:
            if normalized.startswith(prefix):
                return prefix

        return normalized

    def _process_lipid(self, lipid: str):
        try:
            components = self.lipid_parser.parse_lipid(lipid)
            self.lipid_components[lipid] = components
            original_family = components.headgroup
            family = self._normalize_family_name(original_family)

            if family not in self.lipid_families:
                self.lipid_families[family] = set()
            self.lipid_families[family].add(lipid)

            layer = self._determine_layer(components)
            chains = components.get_chains()
            for i, chain in enumerate(chains):
                if components.is_headgroup_only:
                    break
                # Construct node_name without the suffix when i == 0
                node_name = f"{lipid}" if i == 0 else f"{lipid}_{i}"
                color = self.color_manager.assign_color(chain)
                self.nodes[node_name] = NetworkNode(node_name, layer, 0, 0, 0, color)  # Placeholder position

            # Add a central node for the headgroup if it doesn't exist
            headgroup_node_name = f"{family}"
            if headgroup_node_name not in self.nodes:
                self.nodes[headgroup_node_name] = NetworkNode(headgroup_node_name, layer, 0, 0, 0, "#FFFFFF")  # Placeholder position and color

            logging.info(f"Processed lipid: {lipid}, Family: {family}")
        except Exception as e:
            logging.error(f"Error processing lipid {lipid}: {str(e)}")

    def _determine_layer(self, components: LipidComponents) -> str:
        if isinstance(components, (SphingoLipidComponents, SphingoidBaseComponents, AcylCeramideComponents)):
            return "Sphingolipids"
        elif isinstance(components,
                        (GlycerophosphoLipidComponents, LysoGlycerophospholipidComponents, TriAcylComponents,
                         MonoAcylComponents)):
            return "Glycerophospholipids"
        else:
            raise ValueError(f"Unknown lipid type: {type(components)}")

    def _update_enzyme_dict(self, reactants: List[str], product: str, reaction: Reaction):
        # Parse and normalize reactant family names
        reactants_family = [
            self._normalize_family_name(self.lipid_parser.parse_lipid(r).headgroup)
            for r in reactants
        ]

        # Parse and normalize product family name
        product_family = self._normalize_family_name(self.lipid_parser.parse_lipid(product).headgroup)

        # Create the key using normalized family names
        key = (tuple(reactants_family), product_family)

        if reactants_family[0] == 'Cer' and product_family == 'Cer':
            key = (tuple(reactants), product)

        if key not in self.enzyme_dict:
            self.enzyme_dict[key] = set()
        self.enzyme_dict[key].add(reaction.enzyme_name)

    def _update_enzyme_dict_regular(self, reactant: str, product: str, reaction: Reaction):
        reactant_components = self.lipid_parser.parse_lipid(reactant)
        product_components = self.lipid_parser.parse_lipid(product)
        key = (reactant_components.headgroup, product_components.headgroup)

        if key not in self.enzyme_dict:
            self.enzyme_dict[key] = set()
        self.enzyme_dict[key].add(f"{reaction.enzyme_name}")

    def calculate_enzyme_position(
            self, reactants: Tuple[str, ...], product: str,
            family_circle_params: Dict[str, Tuple[float, float, float, float]]) \
            -> Tuple[float, float]:
        if len(reactants) > 1:
            return self._calculate_super_reaction_enzyme_position(reactants, product, family_circle_params)
        else:
            return self._calculate_regular_reaction_enzyme_position(reactants[0], product, family_circle_params)

    def _calculate_regular_reaction_enzyme_position(
            self, reactant: str, product: str, family_circle_params: Dict[str, Tuple[float, float, float, float]]) \
            -> Tuple[float, float]:
        reactant_pos_y = self.nodes[reactant].position_y
        reactant_pos_z = self.nodes[reactant].position_z
        product_pos_y = self.nodes[product].position_y
        product_pos_z = self.nodes[product].position_z
        # reactant_pos = family_circle_params.get(reactant, (0, 0, 10, 20))
        # product_pos = family_circle_params.get(product, (0, 0, 10, 20))

        mid_y = (reactant_pos_y + product_pos_y) / 2
        mid_z = (reactant_pos_z + product_pos_z) / 2

        dy = product_pos_y - reactant_pos_y
        dz = product_pos_z - reactant_pos_z

        length = math.sqrt(dy ** 2 + dz ** 2)
        if length == 0:
            return mid_y, mid_z

        dy, dz = dy / length, dz / length

        perp_dy, perp_dz = -dz, dy

        offset = 20  # Adjust this value to change the distance from the midpoint
        side = 1 if reactant < product else -1

        offset_y = side * perp_dy * offset
        offset_z = side * perp_dz * offset

        return mid_y + offset_y, mid_z + offset_z

    def _calculate_super_reaction_enzyme_position(
            self, reactants: Tuple[str, ...], product: str,
            family_circle_params: Dict[str, Tuple[float, float, float, float]]) \
            -> Tuple[float, float]:
        # Get positions from self.nodes instead of family_circle_params
        positions = [
            (self.nodes[lipid].position_y, self.nodes[lipid].position_z)
            for lipid in reactants + (product,)
        ]

        # Calculate the centroid
        centroid_y = sum(pos[0] for pos in positions) / len(positions)
        centroid_z = sum(pos[1] for pos in positions) / len(positions)

        # Calculate the average distance from the centroid to each point
        avg_distance = sum(
            math.sqrt((pos[0] - centroid_y) ** 2 + (pos[1] - centroid_z) ** 2) for pos in positions) / len(positions)

        # Generate a random angle
        angle = random.uniform(0, 2 * math.pi)

        # Calculate the enzyme position at a fixed distance from the centroid
        offset_distance = avg_distance * 0.5  # Adjust this factor to change how far the enzyme is from the centroid
        enzyme_y = centroid_y + offset_distance * math.cos(angle)
        enzyme_z = centroid_z + offset_distance * math.sin(angle)

        return enzyme_y, enzyme_z

    def _place_nodes(self):
        for family, lipids in self.lipid_families.items():
            circle_params = self.node_placer.calculate_circle_parameters(family, len(lipids))
            self.family_circles[family] = circle_params

            family_components = {lipid: self.lipid_components[lipid] for lipid in lipids}
            positions = self.node_placer.place_nodes(list(lipids), family_components, circle_params)

            for lipid, lipid_positions in positions.items():
                components = self.lipid_components[lipid]
                if components.is_headgroup_only:
                    break
                for i, position in enumerate(lipid_positions):
                    # Construct node_name without the suffix when i == 0
                    node_name = f"{lipid}" if i == 0 else f"{lipid}_{i}"

                    (self.nodes[node_name].position_y, self.nodes[node_name].position_y,
                     self.nodes[node_name].position_z) = position

            # Place the headgroup node at the center of the circle
            headgroup_node_name = f"{family}"
            self.nodes[headgroup_node_name].position_x = 0
            self.nodes[headgroup_node_name].position_y = circle_params.center_y
            self.nodes[headgroup_node_name].position_z = circle_params.center_z


    def _create_edges(self):
        processed_enzymes = set()
        for (reactants_family, product_family), enzymes in self.enzyme_dict.items():
            enzyme_node = f"{';'.join(sorted(enzymes))}"

            # Calculate enzyme position if not already processed
            if enzyme_node not in processed_enzymes:
                positions = self.enzyme_position_manager.get_or_calculate_positions(
                    reactants_family, product_family, self.family_circle_params)
                is_forward = self.enzyme_position_manager.is_forward_reaction(reactants_family, product_family)
                enzyme_y, enzyme_z = positions[:2] if is_forward else positions[2:]
                self.nodes[enzyme_node] = NetworkNode(enzyme_node, "Enzymes", 0, enzyme_y, enzyme_z, "#377eb8")
                processed_enzymes.add(enzyme_node)

            # Get channel for this reaction
            edge_channel = self.enzyme_position_manager.get_reaction_channel(reactants_family, product_family, enzyme_node)
            edge_color = '#FF0000' if edge_channel == "Forward" else '#FFFFFF'

            # Create a unique key for this reaction including the enzyme
            reaction_key = (tuple(sorted(reactants_family)), enzyme_node)

            # Create edges from reactants to enzyme if not already processed
            if reaction_key not in self.processed_reaction_edges:
                for reactant in reactants_family:
                    # src_layer = self._determine_layer(self.lipid_parser.parse_lipid(reactant))
                    self.edges.append(NetworkEdge(f"{self.nodes[reactant].name}_{self.nodes[reactant].layer}",
                                                  f"{self.nodes[enzyme_node].name}_{self.nodes[enzyme_node].layer}",
                                                  channel=edge_channel, color=edge_color))
                self.processed_reaction_edges.add(reaction_key)

            # Always add edge from enzyme to product
            # trg_layer = self._determine_layer(self.lipid_parser.parse_lipid(product_family))
            self.edges.append(NetworkEdge(f"{self.nodes[enzyme_node].name}_{self.nodes[enzyme_node].layer}",
                                          f"{self.nodes[product_family].name}_{self.nodes[product_family].layer}",
                                          channel=edge_channel,
                                          color=edge_color))


    def _update_json(self):
        self.json_data['nodes'] = [
            {
                "name": node.name,
                "layer": node.layer,
                "position_x": node.position_x,
                "position_y": node.position_y,
                "position_z": node.position_z,
                "scale": node.scale,
                "color": node.color,
                "url": node.url,
                "descr": node.descr
            } for node in self.nodes.values()
        ]
        self.json_data['edges'] = [edge.__dict__ for edge in self.edges]


class NodePlacer:
    def __init__(self, family_circle_params: Dict[str, Tuple[float, float, float, float]]):
        self.family_circle_params = family_circle_params

    def calculate_circle_parameters(self, family: str, num_lipids: int) -> CircleParameters:
        y, z, min_radius, max_radius = self.family_circle_params.get(family, (0, 0, 10, 20))

        # Calculate radius based on number of lipids
        circumference = num_lipids * 20  # Assuming 20 is the node_spacing
        radius = circumference / (2 * math.pi)

        # Clamp radius between min and max
        radius = max(min_radius, min(max_radius, radius))

        return CircleParameters(center_y=y, center_z=z, radius=radius)

    def place_nodes(self, lipids: List[str], components_dict: Dict[str, LipidComponents], circle_params: CircleParameters) -> Dict[str, List[Tuple[float, float, float]]]:
        positions = {}
        num_lipids = len(lipids)

        for i, lipid in enumerate(lipids):
            angle = 2 * math.pi * i / num_lipids
            components = components_dict[lipid]
            num_chains = len(components.get_chains())

            lipid_positions = self._calculate_lipid_positions(
                circle_params.center_y,
                circle_params.center_z,
                circle_params.radius,
                angle,
                num_chains
            )

            positions[lipid] = lipid_positions

        return positions

    def _calculate_lipid_positions(self, center_y: float, center_z: float, radius: float, angle: float, num_chains: int,
                                   offset: float = 6) -> List[Tuple[float, float, float]]:
        center_y_on_circle = center_y + radius * math.cos(angle)
        center_z_on_circle = center_z + radius * math.sin(angle)

        if num_chains == 1:
            return [(0, center_y_on_circle, center_z_on_circle)]

        positions = []
        if num_chains == 2:
            z1 = center_z_on_circle - offset
            z2 = center_z_on_circle + offset
            positions = [(0, center_y_on_circle, z1), (0, center_y_on_circle, z2)]
        elif num_chains == 3:
            z1 = center_z_on_circle - offset
            z2 = center_z_on_circle
            z3 = center_z_on_circle + offset
            positions = [(0, center_y_on_circle, z1), (0, center_y_on_circle, z2), (0, center_y_on_circle, z3)]

        return positions

class ColorManager:
    def __init__(self):
        self.common_chains = {
            "16:0": "#FF0000",  # Red
            "18:0": "#00FF00",  # Green
            "18:1": "#0000FF",  # Blue
            "20:4": "#FFFF00",  # Yellow
            "0:0": "#808080"    # Grey
        }

    def assign_color(self, chain: str) -> str:
        if chain.startswith(('m', 'd', 't')):  # Removing d from Sphingolipid backbones
            chain = chain[1:]
        if '-D' in chain:
            chain = chain.split('-D')[0]
        if chain.startswith(('P-', 'O-')):
            chain = chain[2:]
        if chain in self.common_chains:
            return self.common_chains[chain]

        A, B = map(int, chain.split(':'))
        hue = (A % 36) / 36.0
        saturation = 0.5 + (B / (2 * A))
        value = 1.0 - (B / (2 * A))

        r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
        return f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"


class EnzymePositionManager:
    FORWARD = "Forward"
    REVERSE = "Reverse"

    def __init__(self, nodes: Dict[str, NetworkNode]):
        self.positions: Dict[Tuple[str, str], Tuple[float, float, float, float]] = {}
        self.processed_edge_colors: Dict[Tuple[Tuple[str, ...], str], str] = {}
        self.nodes: Dict[str, NetworkNode] = nodes

    def get_or_calculate_positions(
            self, reactants: Tuple[str, ...], product: str,
            family_circle_params: Dict[str, Tuple[float, float, float, float]]) -> Tuple[float, float, float, float]:

        key = tuple(sorted(reactants + (product,)))
        if key in self.positions:
            return self.positions[key]

        if len(reactants) > 1:
            positions = self._calculate_super_reaction_enzyme_position(reactants, product, family_circle_params)
        else:
            positions = self._calculate_regular_reaction_enzyme_position(reactants[0], product, family_circle_params)

        self.positions[key] = positions
        return positions

    def _calculate_regular_reaction_enzyme_position(
            self, reactant: str, product: str,
            family_circle_params: Dict[str, Tuple[float, float, float, float]]) -> Tuple[float, float, float, float]:
        reactant_pos_y = self.nodes[reactant].position_y
        reactant_pos_z = self.nodes[reactant].position_z
        product_pos_y = self.nodes[product].position_y
        product_pos_z = self.nodes[product].position_z

        mid_y = (reactant_pos_y + product_pos_y) / 2
        mid_z = (reactant_pos_z + product_pos_z) / 2

        dy = product_pos_y - reactant_pos_y
        dz = product_pos_z - reactant_pos_z

        length = math.sqrt(dy ** 2 + dz ** 2)
        if length == 0:
            return mid_y, mid_z, mid_y, mid_z

        dy, dz = dy / length, dz / length

        perp_dy, perp_dz = -dz, dy

        offset = 20  # Adjust this value to change the distance from the midpoint

        forward_y = mid_y + perp_dy * offset
        forward_z = mid_z + perp_dz * offset
        reverse_y = mid_y - perp_dy * offset
        reverse_z = mid_z - perp_dz * offset

        return forward_y, forward_z, reverse_y, reverse_z

    def _calculate_super_reaction_enzyme_position(
            self, reactants: Tuple[str, ...], product: str,
            family_circle_params: Dict[str, Tuple[float, float, float, float]]) -> Tuple[float, float, float, float]:
        # Get positions from self.nodes instead of family_circle_params
        positions = [
            (self.nodes[lipid].position_y, self.nodes[lipid].position_z)
            for lipid in reactants + (product,)
        ]

        # Calculate the centroid
        centroid_y = sum(pos[0] for pos in positions) / len(positions)
        centroid_z = sum(pos[1] for pos in positions) / len(positions)

        # Calculate the average distance from the centroid to each point
        avg_distance = sum(
            math.sqrt((pos[0] - centroid_y) ** 2 + (pos[1] - centroid_z) ** 2) for pos in positions) / len(positions)

        # Generate two opposite angles
        angle = random.uniform(0, 2 * math.pi)
        opposite_angle = (angle + math.pi) % (2 * math.pi)

        # Calculate the enzyme positions at a fixed distance from the centroid
        offset_distance = avg_distance * 0.5  # Adjust this factor to change how far the enzyme is from the centroid
        forward_y = centroid_y + offset_distance * math.cos(angle)
        forward_z = centroid_z + offset_distance * math.sin(angle)
        reverse_y = centroid_y + offset_distance * math.cos(opposite_angle)
        reverse_z = centroid_z + offset_distance * math.sin(opposite_angle)

        return forward_y, forward_z, reverse_y, reverse_z

    def get_reaction_channel(self, reactants: Tuple[str, ...], product: str, enzyme_node: str) -> str:
        reaction_key = (tuple(sorted(reactants)), enzyme_node)
        if reaction_key not in self.processed_edge_colors:
            # Determine color based on the reactants and first product
            is_forward = min(reactants) < product
            self.processed_edge_colors[reaction_key] = self.FORWARD if is_forward else self.REVERSE
        return self.processed_edge_colors[reaction_key]

    def is_forward_reaction(self, reactants: Tuple[str, ...], product: str) -> bool:
        return min(reactants) < product
