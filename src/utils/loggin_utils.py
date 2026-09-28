# Function to set up logging (to be used in the main script)
import logging
import sys # Needed for StreamHandler to output to console
from pathlib import Path
from datetime import datetime

def setup_logging(log_dir: str, log_level=logging.INFO):
    """Configures logging to file and console."""
    log_dir_path = Path(log_dir)
    log_dir_path.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir_path / f"db_update_{timestamp}.log"

    # Configure root logger
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'), # Log to file
            logging.StreamHandler(sys.stdout) # Log to console
        ]
    )
    logging.getLogger('psycopg2').setLevel(logging.WARNING) # Quiet down noisy libraries if needed
    logging.info(f"Logging initialized. Log file: {log_file}")