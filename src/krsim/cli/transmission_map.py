"""krsim transmission-map — interactive HTML transmission map."""

import os
import time
import logging

from krsim.cli._common import format_elapsed, require_lineage_files

HELP = "Interactive HTML map of who infected whom"

EPILOG = """\
Outputs:
  always  transmission_map_<seed>.html

Example:
  krsim transmission-map --lineage_csv results/lineage_42.csv \\
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
    from krsim.visualization.transmission_map import save_transmission_map

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    require_lineage_files(args.lineage_csv, args.node_id_map_csv)

    lineage = pd.read_csv(args.lineage_csv)
    node_map = pd.read_csv(args.node_id_map_csv)

    out_path = save_transmission_map(lineage, node_map, sim_id=args.random_seed, output_dir=args.output_dir)
    logging.info(f"### transmission map - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"krsim transmission-map COMPLETED in {format_elapsed(elapsed)}")
