"""krsim relate-host — per-host average of the true strain IBD fraction."""

import os
import time
import logging

from krsim.cli._common import format_elapsed

HELP = "True average IBD between the strains in each current host"

EPILOG = """\
Outputs:
  always  relate_host_<seed>.csv

Example:
  krsim relate-host --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Seed of the run (finds and names its files)")

    parser.add_argument("--relate_csv", type=str, default=None, metavar="FILE",
                        help="Strain IBD table (default: <output_dir>/relate_strain_<seed>.csv)")
    parser.add_argument("--ind_node_map_csv", type=str, default=None, metavar="FILE",
                        help="Host-strain map (default: <output_dir>/sim_ind_node_map_<seed>.csv)")
    parser.add_argument("--output_dir", type=str, default=".", metavar="DIR",
                        help="Where to read and write files (default: current directory)")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import pandas as pd
    from krsim.stats.relate_host import host_average_ibd

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    relate_csv = args.relate_csv
    if relate_csv is None:
        relate_csv = os.path.join(args.output_dir, f"relate_strain_{args.random_seed}.csv")

    ind_node_map_csv = args.ind_node_map_csv
    if ind_node_map_csv is None:
        ind_node_map_csv = os.path.join(args.output_dir, f"sim_ind_node_map_{args.random_seed}.csv")

    host_relatedness = host_average_ibd(
        pd.read_csv(relate_csv), pd.read_csv(ind_node_map_csv), sim_id=args.random_seed
    )

    out_path = os.path.join(args.output_dir, f"relate_host_{args.random_seed}.csv")
    host_relatedness.to_csv(out_path, index=False)
    logging.info(f"### host relatedness - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"krsim relate-host COMPLETED in {format_elapsed(elapsed)}")
