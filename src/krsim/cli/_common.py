"""Helpers shared by the CLI subcommands."""

import os
import sys
import logging


def format_elapsed(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}m"
    else:
        return f"{seconds / 3600:.2f}h"


def require_lineage_files(*paths):
    """Exit with a clear message if --save_lineage outputs are missing."""
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        sys.exit(
            "Missing pedigree file(s): " + ", ".join(missing) + "\n"
            "Re-run `krsim simulate` with --save_lineage to create them."
        )


def setup_run_logging(output_dir, random_seed):
    """
    Full detail (including DEBUG) to KR_sim_<seed>.log in `output_dir`, only
    INFO-and-above to the terminal, so the terminal stays a clean list of
    "what got created" while the file keeps the full record.
    """
    log_path = os.path.join(output_dir, f"KR_sim_{random_seed}.log")

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter("%(asctime)s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(logging.Formatter("%(message)s"))

    logging.basicConfig(level=logging.DEBUG, handlers=[file_handler, stream_handler])
