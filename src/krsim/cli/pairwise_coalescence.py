"""krsim pairwise-coalescence — interactive HTML pedigree coalescence trace."""

import os
import time
import logging

from krsim.cli._common import format_elapsed, require_lineage_files

HELP = "Interactive HTML trace of two strains back to shared pedigree ancestors"

EPILOG = """\
Outputs:
  always  pairwise_coalescence_<seed>.html

Example:
  krsim pairwise-coalescence --lineage_csv results/lineage_42.csv \\
      --node_id_map_csv results/node_id_map_42.csv \\
      --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--lineage_csv", type=str, required=True, metavar="FILE",
                          help="Pedigree table from `krsim simulate --save_lineage` (lineage_<seed>.csv)")
    required.add_argument("--node_id_map_csv", type=str, required=True, metavar="FILE",
                          help="Node ID map from `krsim simulate --save_lineage` (node_id_map_<seed>.csv)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Seed of the run (names the output file)")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Where to write the output file")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import pandas as pd
    from krsim.visualization.pairwise_coalescence import save_pairwise_coalescence

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    require_lineage_files(args.lineage_csv, args.node_id_map_csv)

    lineage = pd.read_csv(args.lineage_csv)
    node_map = pd.read_csv(args.node_id_map_csv)

    out_path = save_pairwise_coalescence(lineage, node_map, sim_id=args.random_seed, output_dir=args.output_dir)
    logging.info(f"### pairwise coalescence - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"krsim pairwise-coalescence COMPLETED in {format_elapsed(elapsed)}")
