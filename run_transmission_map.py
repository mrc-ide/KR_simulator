import os
import time
import logging
import argparse
import pandas as pd

from output.transmission_map import save_transmission_map


def parse_args():
    parser = argparse.ArgumentParser(
        description="Render an interactive HTML transmission map for one simulation run"
    )

    # --- required parameters ---
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility (also used as sim_id)")

    # --- optional parameters ---
    parser.add_argument("--lineage_csv", type=str, default=None,
                        help="Path to lineage CSV (default: <output_dir>/lineage_<random_seed>.csv)")
    parser.add_argument("--node_id_map_csv", type=str, default=None,
                        help="Path to node ID map CSV (default: <output_dir>/node_id_map_<random_seed>.csv)")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Directory to write output files to (default: current directory).")

    return parser.parse_args()


def _format_elapsed(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}m"
    else:
        return f"{seconds / 3600:.2f}h"


def main():
    args = parse_args()
    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    lineage_csv = args.lineage_csv or os.path.join(args.output_dir, f"lineage_{args.random_seed}.csv")
    node_id_map_csv = args.node_id_map_csv or os.path.join(args.output_dir, f"node_id_map_{args.random_seed}.csv")

    lineage = pd.read_csv(lineage_csv)
    node_map = pd.read_csv(node_id_map_csv)

    out_path = save_transmission_map(lineage, node_map, sim_id=args.random_seed, output_dir=args.output_dir)
    logging.info(f"### transmission map - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"run_transmission_map.py COMPLETED in {_format_elapsed(elapsed)}")


if __name__ == "__main__":
    main()
