import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
import signal


# --- Path Setup ---
def setup_project_path():
    """Setup the project path to enable imports from src/"""
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir  # Assuming script is in project root
    src_path = project_root / "src"

    if not src_path.exists():
        raise FileNotFoundError(f"Source directory not found: {src_path}")

    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))

    return project_root


# Setup path before imports
try:
    PROJECT_ROOT = setup_project_path()
except FileNotFoundError as e:
    print(f"Error: {e}", file=sys.stderr)
    sys.exit(1)

from src.data.config_loader import ConfigLoader
from src.utils.db_updater import LipidDatabaseUpdater, DryRunRollbackException

# Global logger for the CLI script
logger = logging.getLogger(__name__)


def setup_cli_logging(log_level: str, project_root: Path):
    numeric_log_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_log_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    log_file_dir = project_root / "logs"
    log_file_dir.mkdir(parents=True, exist_ok=True)
    log_file_path = log_file_dir / f"updater_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    # Ensure handlers use UTF-8
    handlers = []

    # Console Handler (StreamHandler)
    # Try to make stdout UTF-8 compatible if it isn't
    # This is a common issue on Windows.
    try:
        # If stdout is a tty and encoding is not UTF-8, try to reconfigure
        if sys.stdout.isatty() and getattr(sys.stdout, 'encoding', '').lower() != 'utf-8':
            # This can be tricky and platform-dependent.
            # PyCharm usually handles this well. If running from cmd.exe on Windows,
            # you might need `chcp 65001` before running the script.
            # Forcing sys.stdout.reconfigure might work in some Python versions/OS combos.
            # A simpler approach for logging is to specify encoding on the handler.
            pass  # Let PyCharm handle its console, FileHandler is more critical for encoding
    except Exception:
        pass  # Best effort

    console_handler = logging.StreamHandler(sys.stdout)
    # For StreamHandler, Python tries to use sys.stdout.encoding.
    # If that's problematic (like cp1252), the error occurs.
    # We can't directly set encoding on StreamHandler easily if the underlying stream doesn't support it.
    # The FileHandler is where we have more control.

    file_handler = logging.FileHandler(log_file_path, encoding='utf-8')  # Explicitly use UTF-8 for file

    handlers.append(console_handler)
    handlers.append(file_handler)

    logging.basicConfig(
        level=numeric_log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=handlers,
    )
    # For good measure, also ensure the root logger's handlers are set if basicConfig was called before
    # logging.getLogger().handlers = handlers
    # This line above might be too aggressive, basicConfig should suffice if called early.

    logger.info(f"Logging to console and to file: {log_file_path} (File encoding: UTF-8)")


def validate_file_exists(file_path: str, description: str) -> Path:
    """Validate that a file exists and return Path object"""
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {path}")
    if not path.is_file():
        raise ValueError(f"{description} is not a file: {path}")
    return path


def validate_directory_exists(dir_path: str, description: str) -> Path:
    """Validate that a directory exists and return Path object"""
    path = Path(dir_path)
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {path}")
    if not path.is_dir():
        raise ValueError(f"{description} is not a directory: {path}")
    return path


def handle_generate_command(args, updater: LipidDatabaseUpdater) -> None:
    """Handle the generate command"""
    # Validate input files
    lipids_file = validate_file_exists(args.lipids_file, "Lipids file")
    enzymes_file = validate_file_exists(args.enzymes_file, "Enzymes file")

    lipids2uniprot_file = None
    if args.lipids2uniprot_file:
        lipids2uniprot_file = validate_file_exists(
            args.lipids2uniprot_file, "Lipids2UniProt file"
        )

    # Determine run timestamp
    effective_run_timestamp = (
        args.run_timestamp
        if args.run_timestamp
        else datetime.now().strftime("%Y%m%d_%H%M%S_auto")
    )

    logger.info(
        f"Executing 'generate' command. Run timestamp: {effective_run_timestamp}"
    )
    logger.info(f"Input files - Lipids: {lipids_file}, Enzymes: {enzymes_file}")
    if lipids2uniprot_file:
        logger.info(f"Lipids2UniProt: {lipids2uniprot_file}")

    summary = updater.process_files_for_review(
        lipids_file_path=str(lipids_file),
        enzymes_file_path=str(enzymes_file),
        lipids2uniprot_file_path=(
            str(lipids2uniprot_file) if lipids2uniprot_file else None
        ),
        current_data_run_timestamp=effective_run_timestamp,
    )

    output_dir = updater.output_dir / effective_run_timestamp
    logger.info(f"'generate' process completed. Summary: {summary}")
    logger.info(f"Review files generated in: {output_dir}")
    logger.info("Next steps:")
    logger.info(f"1. Review the delta files in: {output_dir}")
    logger.info("2. Edit/approve the files as needed")
    logger.info("3. Run the 'apply' command with the approved files")


def handle_apply_command(args, updater: LipidDatabaseUpdater) -> None:
    """Handle the apply command"""
    approved_dir = validate_directory_exists(args.approved_dir, "Approved directory")

    logger.info(f"Executing 'apply' command")
    logger.info(f"Approved files from: {approved_dir}")
    logger.info(f"Original run timestamp: {args.run_timestamp}")
    logger.info(f"Commit mode: {'COMMIT' if args.commit else 'DRY RUN'}")

    if not args.commit:
        logger.warning("Running in DRY RUN mode - no changes will be committed")

    updater.apply_all_approved_updates(
        approved_delta_dir_path=approved_dir,
        run_timestamp_of_deltas=args.run_timestamp,
        commit_changes=args.commit,
    )

    if args.commit:
        logger.info("'apply' process completed and changes COMMITTED to database.")
    else:
        logger.info(
            "'apply' process completed in DRY RUN mode. No changes were committed."
        )

def handle_report_command(args, updater: LipidDatabaseUpdater):
    """Handle the report command"""
    logger.info("Executing 'report' command...")
    if args.report_type == 'missing-abbr':
        updater.report_missing_abbreviations(args.output_file)
    else:
        logger.error(f"Unknown report type: {args.report_type}")


def create_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser"""
    parser = argparse.ArgumentParser(
        description="LipidCRED Database Updater CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Examples:
          # Generate delta files for review
          python update_lipidcred_db.py generate --lipids-file data/lipids.tsv --enzymes-file data/enzymes.tsv

          # Apply approved changes (dry run)
          python update_lipidcred_db.py apply --approved-dir output/20250523_143022_auto --run-timestamp 20250523_143022_auto

          # Apply approved changes (commit to database)
          python update_lipidcred_db.py apply --approved-dir output/20250523_143022_auto --run-timestamp 20250523_143022_auto --commit
                """,
    )

    parser.add_argument(
        "--config-file",
        type=str,
        default="config/lipid_config.yaml",
        help="Path to the configuration file (default: config/lipid_config.yaml)",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO)",
    )

    subparsers = parser.add_subparsers(
        dest="command", required=True, help="Available commands", metavar="COMMAND"
    )

    # Generate command
    parser_generate = subparsers.add_parser(
        "generate",
        help="Ingest SwissLipids data and generate delta files for review",
        description="Process SwissLipids TSV files and generate delta files for manual review",
    )
    parser_generate.add_argument(
        "--lipids-file",
        type=str,
        required=True,
        help="Path to the SwissLipids lipids.tsv file",
    )
    parser_generate.add_argument(
        "--enzymes-file",
        type=str,
        required=True,
        help="Path to the SwissLipids enzymes.tsv file",
    )
    parser_generate.add_argument(
        "--lipids2uniprot-file",
        type=str,
        help="Optional: Path to the SwissLipids lipids2uniprot.tsv file",
    )
    parser_generate.add_argument(
        "--run-timestamp",
        type=str,
        help="Optional: Custom timestamp for output directory (e.g., YYYYMMDD_HHMMSS_custom). If not provided, current datetime is used",
    )

    # Apply command
    parser_apply = subparsers.add_parser(
        "apply",
        help="Apply approved delta files to the database",
        description="Apply manually reviewed and approved delta files to the database",
    )
    parser_apply.add_argument(
        "--approved-dir",
        type=str,
        required=True,
        help="Path to directory containing approved delta files (approved_delta_*.csv/tsv)",
    )
    parser_apply.add_argument(
        "--run-timestamp",
        type=str,
        required=True,
        help="Timestamp of the original 'generate' run (used for snapshot continuity)",
    )
    parser_apply.add_argument(
        "--commit",
        action="store_true",
        help="Commit changes to database. Without this flag, performs a dry run (default: dry run)",
    )

    # Report command
    parser_report = subparsers.add_parser(
        "report",
        help="Generate data quality reports from the production database",
    )
    parser_report.add_argument(
        "report_type",
        choices=['missing-abbr'],
        help="The type of report to generate.",
    )
    parser_report.add_argument(
        "--output-file",
        type=str,
        default="reports/missing_abbreviations.tsv",
        help="Path to save the report file.",
    )

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    # Setup logging
    try:
        setup_cli_logging(args.log_level, PROJECT_ROOT)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    logger.info(f"LipidCRED Database Updater CLI started - Command: {args.command}")
    logger.debug(f"Arguments: {args}")

    updater = None

    def cleanup_handler(sig, frame):
        signal_name = "SIGINT" if sig == signal.SIGINT else "SIGTERM"
        logger.warning(f"Received {signal_name} signal. Cleaning up database connections...")

        if updater and hasattr(updater, 'conn') and updater.conn and not updater.conn.closed:
            try:
                # Try to restore synchronous_commit if it was changed
                try:
                    cursor = updater.conn.cursor()
                    cursor.execute("SET synchronous_commit = ON;")
                except:
                    pass

                updater.conn.close()
                logger.info("Database connection closed gracefully.")
            except Exception as e:
                logger.warning(f"Error closing database connection: {e}")

        logger.info("Cleanup completed. Exiting.")
        sys.exit(130)  # Standard exit code for SIGINT

    # Install signal handlers
    signal.signal(signal.SIGINT, cleanup_handler)  # Ctrl+C
    signal.signal(signal.SIGTERM, cleanup_handler)  # Kill command

    try:
        # Load configuration
        config_path = validate_file_exists(args.config_file, "Configuration file")
        config = ConfigLoader(str(config_path))
        logger.info(f"Configuration loaded from: {config_path}")

        # Execute command with database updater
        with LipidDatabaseUpdater(config) as updater:
            logger.info("Database connection established")

            if args.command == "generate":
                handle_generate_command(args, updater)
            elif args.command == "apply":
                handle_apply_command(args, updater)
            elif args.command == "report":
                handle_report_command(args, updater)

        logger.info("Database connection closed successfully")

    except ConnectionError as e:
        logger.critical(f"Database connection error: {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        logger.critical(f"File not found: {e}")
        sys.exit(1)
    except ValueError as e:
        logger.critical(f"Configuration error: {e}")
        sys.exit(1)
    except DryRunRollbackException as e:
        logger.info(f"Dry run completed successfully: {e}")
    except KeyboardInterrupt:
        logger.warning("Operation interrupted by user")
        sys.exit(130)  # Standard exit code for SIGINT
    except Exception as e:
        logger.critical(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)
    finally:
        logger.info("CLI execution finished")


if __name__ == "__main__":
    main()
