"""krsim relate-strain — true pairwise IBD fraction between strains."""

import os
import time
import logging

from krsim.cli._common import format_elapsed

HELP = "True pairwise strain IBD fraction from the coalescence time table"


def add_arguments(parser):
    # --- required parameters ---
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility (also used as sim_id)")
    parser.add_argument("--time_boundary", type=float, required=True,
                        help="Coalescence time threshold; segments with coal_time <= this are counted as IBD")

    # --- optional parameters ---
    parser.add_argument("--coal_csv", type=str, default=None,
                        help="Path to coal_T CSV (default: <output_dir>/coal_T_<random_seed>.csv)")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Directory to write output files to (default: current directory).")


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
