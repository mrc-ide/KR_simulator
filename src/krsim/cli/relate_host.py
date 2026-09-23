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
  krsim relate-host --relate_csv results/relate_strain_42.csv \\
      --ind_node_map_csv results/sim_ind_node_map_42.csv \\
      --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--relate_csv", type=str, required=True, metavar="FILE",
                          help="Strain IBD table from `krsim relate-strain` (relate_strain_<seed>.csv)")
    required.add_argument("--ind_node_map_csv", type=str, required=True, metavar="FILE",
                          help="Host-strain map from `krsim simulate` (sim_ind_node_map_<seed>.csv)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import pandas as pd
    from krsim.stats.relate_host import host_average_ibd

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    host_relatedness = host_average_ibd(
        pd.read_csv(args.relate_csv), pd.read_csv(args.ind_node_map_csv), sim_id=args.random_seed
    )

    out_path = os.path.join(args.output_dir, f"relate_host_{args.random_seed}.csv")
    host_relatedness.to_csv(out_path, index=False)
    logging.info(f"### host relatedness - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"krsim relate-host COMPLETED in {format_elapsed(elapsed)}")
