"""krsim relate-strain — true pairwise IBD fraction between strains."""

import os
import time
import logging

from krsim.cli._common import format_elapsed

HELP = "True pairwise IBD fraction between strains"

EPILOG = """\
Outputs:
  always  relate_strain_<seed>.csv

Example:
  krsim relate-strain --random_seed 42 --time_boundary 2 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Seed of the run (finds and names its files)")
    required.add_argument("--time_boundary", type=float, required=True, metavar="T",
                          help="Segments that coalesce within T generations count as IBD")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Where to read and write files")

    parser.add_argument("--coal_csv", type=str, default=None, metavar="FILE",
                        help="Coalescence table (default: <output_dir>/coal_T_<seed>.csv)")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import pandas as pd
    from krsim.stats.relate_strain import strain_ibd_fractions

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    coal_csv = args.coal_csv
    if coal_csv is None:
        coal_csv = os.path.join(args.output_dir, f"coal_T_{args.random_seed}.csv")

    pairs = strain_ibd_fractions(pd.read_csv(coal_csv), args.time_boundary)

    out_path = os.path.join(args.output_dir, f"relate_strain_{args.random_seed}.csv")
    pairs.to_csv(out_path, index=False)
    logging.info(f"### IBD fractions - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"krsim relate-strain COMPLETED in {format_elapsed(elapsed)}")
