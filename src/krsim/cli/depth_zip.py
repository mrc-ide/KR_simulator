"""krsim depth-zip — helper: mosdepth per-base files -> one zip of depth profiles for `vcf-host-depth`."""

import os
import sys
import time
import logging
import tempfile

from krsim.cli._common import KrsimError, format_elapsed

HELP = "mosdepth per-base files -> zip for vcf-host-depth"

EPILOG = """\
Turns each sample's mosdepth per-base output (<sample>.per-base.bed.gz) into its
depth profile, the mean and SD of read depth in every 500 bp window
(<sample>.window_stats.tsv), and puts them all in one zip. That zip is the input
for `krsim vcf-host-depth --depth_zip`.

Run mosdepth against the 3D7 reference used by MalariaGEN Pf8 (PlasmoDB-54), the
genome krsim simulates; files mapped to another reference are rejected. Each
sample takes a few seconds and up to about 1 GB of memory.

Outputs:
  always  the zip given by --output

Example:
  mosdepth --fasta PlasmoDB-54-Pfalciparum3D7-Genome.fasta sample1 sample1.cram
  krsim depth-zip --bed sample1.per-base.bed.gz sample2.per-base.bed.gz \\
      --output my_depth.zip
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--bed", type=str, nargs="+", required=True, metavar="FILE",
                          help="mosdepth per-base output, one per sample (<sample>.per-base.bed.gz)")
    required.add_argument("--output", type=str, required=True, metavar="FILE",
                          help="Zip to write (e.g. my_depth.zip)")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    from krsim.output.pf8_profiles import make_depth_zip

    if not args.output.endswith(".zip"):
        raise KrsimError(f"--output must end in .zip, got {args.output}")
    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)

    # progress on screen; this helper has no seed, so no KR_sim_<seed>.log
    log = logging.getLogger("krsim")
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("  %(message)s"))
    log.handlers = [handler]
    log.setLevel(logging.INFO)
    log.propagate = False

    start = time.time()
    print(file=sys.stderr)
    try:
        with tempfile.TemporaryDirectory(dir=out_dir) as work_dir:
            samples = make_depth_zip(args.bed, args.output, work_dir)
    except ValueError as exc:
        raise KrsimError(str(exc)) from exc
    finally:
        if os.path.exists(args.output + ".part"):
            os.remove(args.output + ".part")
    log.info(f"wrote  {args.output} ({len(samples)} sample{'s' if len(samples) != 1 else ''})")
    log.info(f"done in {format_elapsed(time.time() - start)}")
    print(file=sys.stderr)
    log.handlers = []
