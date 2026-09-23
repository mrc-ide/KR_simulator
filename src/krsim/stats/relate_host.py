import os
import time
import logging
import argparse
import itertools
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute per-host average pairwise IBD fraction from strain-level relatedness, "
                     "restricted to the most recent (currently sampled) human hosts"
    )

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

    relate_csv = args.relate_csv
    if relate_csv is None:
        relate_csv = os.path.join(args.output_dir, f"relate_strain_{args.random_seed}.csv")

    ind_node_map_csv = args.ind_node_map_csv
    if ind_node_map_csv is None:
        ind_node_map_csv = os.path.join(args.output_dir, f"sim_ind_node_map_{args.random_seed}.csv")

    # --- node_id -> frac_IBD lookup, keyed by sorted (sample_1, sample_2) ---
    relate = pd.read_csv(relate_csv)
    pair_to_frac = {
        (int(row.sample_1), int(row.sample_2)): row.frac_IBD
        for row in relate.itertuples()
    }

    # --- ind_id -> node_ids for this sim_id, restricted to the most recent generation (time=0) ---
    ind_node_map = pd.read_csv(ind_node_map_csv)
    ind_node_map = ind_node_map[
        (ind_node_map["sim_id"] == args.random_seed) & (ind_node_map["time"] == 0)
    ]
    ind_to_nodes = ind_node_map.groupby("ind_id")["node_id"].apply(list).to_dict()

    # --- average pairwise IBD across each host's own strain nodes ---
    # pairs with no recorded coalescence anywhere in the ARG are assumed to have frac_IBD=0
    rows = []
    for ind_id, node_ids in sorted(ind_to_nodes.items()):
        fracs = [
            pair_to_frac.get((node_a, node_b), 0.0)
            for node_a, node_b in itertools.combinations(sorted(node_ids), 2)
        ]

        avg_frac_IBD = sum(fracs) / len(fracs) if fracs else float("nan")
        rows.append({"host": ind_id, "avg_frac_IBD": avg_frac_IBD})

    host_relatedness = pd.DataFrame(rows)

    out_path = os.path.join(args.output_dir, f"relate_host_{args.random_seed}.csv")
    host_relatedness.to_csv(out_path, index=False)
    logging.info(f"### host relatedness - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"relate_t_host.py COMPLETED in {_format_elapsed(elapsed)}")


if __name__ == "__main__":
    main()
