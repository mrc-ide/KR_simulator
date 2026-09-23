"""Helpers shared by the CLI subcommands: the run log, and user-facing errors."""

import os
import sys
import time
import shlex
import logging
from datetime import datetime

# library modules log through `logging.getLogger(__name__)`, i.e. under "krsim"
log = logging.getLogger("krsim")


class KrsimError(Exception):
    """An error with a message meant for the user (no traceback on screen)."""


def format_elapsed(seconds):
    if seconds < 60:
        return f"{seconds:.1f} s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f} min"
    else:
        return f"{seconds / 3600:.2f} h"


def require_lineage_files(*paths):
    """Stop with a clear message if --save_lineage outputs are missing."""
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise KrsimError(
            "missing pedigree file(s): " + ", ".join(missing)
            + " (re-run `krsim simulate` with --save_lineage)"
        )


class _Formatter(logging.Formatter):
    """Appends a record's `extra` text: extra detail in the log file, the log path on screen."""

    def __init__(self, fmt, datefmt=None, extra_attr=None):
        super().__init__(fmt, datefmt)
        self.extra_attr = extra_attr

    def format(self, record):
        text = super().format(record)
        extra = getattr(record, self.extra_attr, None) if self.extra_attr else None
        return text + extra if extra else text


class Run:
    """
    One krsim command's run. Short messages go to the screen (stderr); the
    full record goes to KR_sim_<seed>.log in the output directory, one block
    per command, starting with a header that holds the exact command line.
    """

    current = None  # the active run, so main() can log errors into it

    def __init__(self, command, output_dir, random_seed):
        from importlib.metadata import version

        os.makedirs(output_dir, exist_ok=True)
        self.log_path = os.path.join(output_dir, f"KR_sim_{random_seed}.log")
        self.start = time.time()

        with open(self.log_path, "a") as f:
            f.write(f"\n===== {datetime.now():%Y-%m-%d %H:%M:%S}  krsim {version('krsim')}  {command} =====\n"
                    f"command: krsim {shlex.join(sys.argv[1:])}\n")

        file_handler = logging.FileHandler(self.log_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(_Formatter("%(asctime)s  %(message)s", "%H:%M:%S", extra_attr="detail"))

        screen_handler = logging.StreamHandler(sys.stderr)
        screen_handler.setLevel(logging.INFO)
        screen_handler.setFormatter(_Formatter("  %(message)s", extra_attr="screen"))

        log.handlers = [file_handler, screen_handler]
        log.setLevel(logging.DEBUG)
        log.propagate = False

        print(file=sys.stderr)  # blank line under the shell prompt
        Run.current = self

    def detail(self, message):
        """Log file only."""
        log.debug(message)

    def wrote(self, path, rows=None):
        """Screen and log file; the row count goes to the log file only."""
        log.info(f"wrote  {path}", extra={"detail": f" ({rows:,} rows)" if rows is not None else None})

    def done(self):
        log.info(f"done in {format_elapsed(time.time() - self.start)}",
                 extra={"screen": f"  (log: {self.log_path})"})
        self._close()

    def fail(self, message):
        """Error: one line on screen, the full traceback in the log file."""
        log.error(f"error: {message}", extra={"screen": f"\n  (details in {self.log_path})"})
        log.debug("traceback:", exc_info=True)
        self._close()

    def _close(self):
        print(file=sys.stderr)
        for handler in log.handlers:
            handler.close()
        log.handlers = []
        Run.current = None
