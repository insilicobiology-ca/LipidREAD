from typing import List, Dict, Tuple, Set, Optional
from dataclasses import dataclass


@dataclass
class NetworkNode:
    name: str
    layer: str
    position_x: float
    position_y: float
    position_z: float
    color: str
    scale: float = 1
    url: str = ""
    descr: str = ""

@dataclass
class NetworkEdge:
    def __init__(self, src: str, trg: str, opacity: float = 1, channel: str = "", color: str = "#ffffff"):
        self.src = src
        self.trg = trg
        self.opacity = opacity
        self.channel = channel
        self.color = color
