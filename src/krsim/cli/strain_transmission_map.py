"""krsim strain-transmission-map — interactive HTML transmission map."""

import os

from krsim.cli._common import Run, require_lineage_files

HELP = "Interactive HTML map of who infected whom"

EPILOG = """\
Outputs:
  always  strain_transmission_map_<seed>.html

Example:
  krsim strain-transmission-map --lineage_csv results/lineage_42.csv \\
      --node_id_map_csv results/node_id_map_42.csv \\
      --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--lineage_csv", type=str, required=True, metavar="FILE",
                          help="Pedigree from `krsim simulate --save_lineage` (lineage_<seed>.csv)")
    required.add_argument("--node_id_map_csv", type=str, required=True, metavar="FILE",
                          help="Node IDs from `krsim simulate --save_lineage` (node_id_map_<seed>.csv)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import pandas as pd
    from krsim.visualization.strain_transmission_map import save_strain_transmission_map

    run = Run("strain-transmission-map", args.output_dir, args.random_seed)

    require_lineage_files(args.lineage_csv, args.node_id_map_csv)
    lineage = pd.read_csv(args.lineage_csv)
    node_map = pd.read_csv(args.node_id_map_csv)

    out_path = save_strain_transmission_map(lineage, node_map, sim_id=args.random_seed, output_dir=args.output_dir)
    run.wrote(out_path)

    run.done()
