"""krsim transmission-map — interactive HTML transmission map."""

import os
import time
import logging

from krsim.cli._common import format_elapsed, require_lineage_files

HELP = "Render an interactive HTML transmission map for one simulation run"


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    # --- required parameters ---
    required.add_argument("--random_seed", type=int, required=True,
                          help="Random seed for reproducibility (also used as sim_id)")

    # --- optional parameters ---
    parser.add_argument("--lineage_csv", type=str, default=None,
                        help="Path to lineage CSV (default: <output_dir>/lineage_<random_seed>.csv)")
    parser.add_argument("--node_id_map_csv", type=str, default=None,
                        help="Path to node ID map CSV (default: <output_dir>/node_id_map_<random_seed>.csv)")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Directory to write output files to (default: current directory).")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import pandas as pd
    from krsim.visualization.transmission_map import save_transmission_map

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    lineage_csv = args.lineage_csv or os.path.join(args.output_dir, f"lineage_{args.random_seed}.csv")
    node_id_map_csv = args.node_id_map_csv or os.path.join(args.output_dir, f"node_id_map_{args.random_seed}.csv")

    require_lineage_files(lineage_csv, node_id_map_csv)

    lineage = pd.read_csv(lineage_csv)
    node_map = pd.read_csv(node_id_map_csv)

    out_path = save_transmission_map(lineage, node_map, sim_id=args.random_seed, output_dir=args.output_dir)
    logging.info(f"### transmission map - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"krsim transmission-map COMPLETED in {format_elapsed(elapsed)}")
