"""krsim relate-host — per-host average of the true strain IBD fraction."""

import os
import time
import logging

from krsim.cli._common import format_elapsed

HELP = ("Per-host average pairwise IBD fraction from strain-level relatedness, "
        "restricted to the most recent (currently sampled) human hosts")


def add_arguments(parser):
    # --- required parameters ---
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility (also used as sim_id)")

    # --- optional parameters ---
    parser.add_argument("--relate_csv", type=str, default=None,
                        help="Path to relate_strain CSV (default: <output_dir>/relate_strain_<random_seed>.csv)")
    parser.add_argument("--ind_node_map_csv", type=str, default=None,
                        help="Path to sim_ind_node_map CSV (default: <output_dir>/sim_ind_node_map_<random_seed>.csv)")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Directory to write output files to (default: current directory).")


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
