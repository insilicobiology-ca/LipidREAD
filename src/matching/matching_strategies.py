from typing import Dict, List, Callable
import yaml
from enum import Enum
from abc import ABC, abstractmethod


class MatchingStrategy(Enum):
    COMPLETE = 1
    PC_TO_LPC = 2
    UNMARKED_LPC_TO_PC = 3
    LPC_TO_PC = 4
    N_ACYL_MATCH = 5
    BACKBONE_MATCH = 6
    N_ACYL_MATCH_DESATURASE = 7
    NO_MATCH = 8
    SUPER_MATCH = 9
    CER_PC_TO_ACER_LPC = 10
    LYSOGPL_GPL_TO_GPL = 11
    NAPE_PA_MATCH = 12
    NAPE_GPNAE_MATCH = 13
    NAPE_NAE_MATCH = 14
    NAPE_NALPE_MATCH = 15
    PC_PE_TO_NAPE = 16
    PC_PE_TO_PC = 17
    PC_PE_TO_PE = 18
    PCO_SPB_TO_CER = 19
    DG_TO_TG = 20
    PCO_SPH_TO_LPCO_CER = 23








