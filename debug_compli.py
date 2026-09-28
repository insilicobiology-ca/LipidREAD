from src.core.lipid_analysis_workflow import LipidAnalysisWorkflow
import logging
import os
from logging.handlers import RotatingFileHandler
from pathlib import Path


def setup_logging(log_file_path, console_level=logging.INFO, file_level=logging.DEBUG):
    # Create a logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)  # Capture all levels

    # Create formatters
    detailed_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_formatter = logging.Formatter('%(levelname)s - %(message)s')

    # Console Handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(console_level)
    console_handler.setFormatter(console_formatter)

    # File Handler
    file_handler = RotatingFileHandler(log_file_path, maxBytes=300*1024*1024, backupCount=5)
    file_handler.setLevel(file_level)
    file_handler.setFormatter(detailed_formatter)

    # Add handlers to logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

if __name__ == "__main__":
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    # log_file_path = os.path.join(log_dir, "lexicon_analysis_Aug1_run3.log")

    log_file_path = os.path.join(log_dir, "translator_debug.log")
    setup_logging(log_file_path)
    config_file_path = "config/lipid_config.yaml"
    workflow = LipidAnalysisWorkflow(config_file_path, 'debug_files', debug_mode=True)
    # files = ['output_tests/quickTest.xlsx']
    # files = ['debug_files/super.xlsx']
    # files = ['debug_files/circleTest2.xlsx']
    # files = ['debug_files/wrong_name_test.xlsx']
    # files = ['debug_files/list_of_lipids_Jorg_for_compli.xlsx']
    # files = ['debug_files/network_error_test.xlsx']
    # files = ['debug_files/Lipid CRED_SPH Species List.xlsx']
    # files = ['C:\MastersVault\complitest\lexicon\Aug2024CompliList.xlsx']
    # files = ['output_tests/lexicon_list.xlsx']
    file_folder = Path("../www/")
    file = file_folder / "exampleinput.xlsx"
    files = [file]
    # files = ['C:\MastersVault\complitest/bench/C:\MastersVault\complitest\benchmarking\benchmark\linex_comple.xlsx']
    # files = ['C:\MastersVault\complitest/debug_files/complitest_mira2.xlsx']
    organism_id = 9606
    translate = True
    network = False
    lipids = ["Cer(d14:1/16:0)"]
    enzyme = "ASAH1"
    workflow.analyze_single_lipid(lipids, organism_id)
    # workflow.search_reactions_by_enzyme(enzyme, organism_id)
    # workflow.process_files(files, organism_id, translate, network)


# TODO:
#   - DONE dd C1P to molecules, and to matching strategies in lipid_config
#   - DONE Add parsing strategy for SphingolipidBackbone e.g. Sph(d18:1), SM(d18:1), Cer(d18:1)
#   - DONE Clean up lipids (especialy after being translated) from double bond locations
#   - DONE Check LipidMatcher when it returns SUPER_MATCH or NO_MATCH, needs to be a function I believe
#   - DONE Add NAPE parsing strategy
#   - EFFICIENCY: seeing logs sometimes the same reaction tries to be added, maybe can keep cache of them, to skip



# TODO - LONG TERM
#   - Add ability to parse lipids based on both old and newer nomenclature
