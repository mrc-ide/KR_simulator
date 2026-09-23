"""krsim pairwise-strain-pedigree-map — interactive HTML pedigree coalescence trace."""

import os
import time
import logging

from krsim.cli._common import format_elapsed, require_lineage_files

HELP = "Interactive HTML trace of two strains' shared pedigree ancestors"

EPILOG = """\
Outputs:
  always  pairwise_strain_pedigree_map_<seed>.html

Example:
  krsim pairwise-strain-pedigree-map --lineage_csv results/lineage_42.csv \\
      --node_id_map_csv results/node_id_map_42.csv \\
      --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--lineage_csv", type=str, required=True, metavar="FILE",
                          help="Pedigree table (lineage_<seed>.csv, from simulate --save_lineage)")
    required.add_argument("--node_id_map_csv", type=str, required=True, metavar="FILE",
                          help="Node ID map (node_id_map_<seed>.csv, from simulate --save_lineage)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Where to write the output file")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import pandas as pd
    from krsim.visualization.pairwise_strain_pedigree_map import save_pairwise_strain_pedigree_map

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    require_lineage_files(args.lineage_csv, args.node_id_map_csv)

    lineage = pd.read_csv(args.lineage_csv)
    node_map = pd.read_csv(args.node_id_map_csv)

    out_path = save_pairwise_strain_pedigree_map(lineage, node_map, sim_id=args.random_seed, output_dir=args.output_dir)
    logging.info(f"### pairwise strain pedigree map - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"krsim pairwise-strain-pedigree-map COMPLETED in {format_elapsed(elapsed)}")
