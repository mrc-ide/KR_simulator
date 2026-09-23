"""krsim pairwise-coalescence — interactive HTML pedigree coalescence trace."""

import os
import time
import logging

from krsim.cli._common import format_elapsed, require_lineage_files

HELP = "Interactive HTML trace of two strains back to shared pedigree ancestors"

EPILOG = """\
Needs `krsim simulate --save_lineage`.

Outputs:
  always  pairwise_coalescence_<seed>.html

Example:
  krsim pairwise-coalescence --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Seed of the run (finds and names its files)")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Where to read and write files")

    parser.add_argument("--lineage_csv", type=str, default=None, metavar="FILE",
                        help="Pedigree table (default: <output_dir>/lineage_<seed>.csv)")
    parser.add_argument("--node_id_map_csv", type=str, default=None, metavar="FILE",
                        help="Node ID map (default: <output_dir>/node_id_map_<seed>.csv)")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import pandas as pd
    from krsim.visualization.pairwise_coalescence import save_pairwise_coalescence

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    lineage_csv = args.lineage_csv or os.path.join(args.output_dir, f"lineage_{args.random_seed}.csv")
    node_id_map_csv = args.node_id_map_csv or os.path.join(args.output_dir, f"node_id_map_{args.random_seed}.csv")

    require_lineage_files(lineage_csv, node_id_map_csv)

    lineage = pd.read_csv(lineage_csv)
    node_map = pd.read_csv(node_id_map_csv)

    out_path = save_pairwise_coalescence(lineage, node_map, sim_id=args.random_seed, output_dir=args.output_dir)
    logging.info(f"### pairwise coalescence - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"krsim pairwise-coalescence COMPLETED in {format_elapsed(elapsed)}")
