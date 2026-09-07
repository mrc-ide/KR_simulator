import os
import sys
import time
import logging
import argparse
import pandas as pd

# repo root, so `sim` resolves regardless of where this script is invoked from
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sim.params import GENOME_LENGTH


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute pairwise IBD fraction from a coalescence time table"
    )

    # --- required parameters ---
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility (also used as sim_id)")
    parser.add_argument("--time_boundary", type=float, required=True,
                        help="Coalescence time threshold; segments with coal_time <= this are counted as IBD")

    # --- optional parameters ---
    parser.add_argument("--coal_csv", type=str, default=None,
                        help="Path to coal_T CSV (default: <output_dir>/coal_T_<random_seed>.csv)")
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

    coal_csv = args.coal_csv
    if coal_csv is None:
        coal_csv = os.path.join(args.output_dir, f"coal_T_{args.random_seed}.csv")

    coal = pd.read_csv(coal_csv)

    # segments whose coalescence time falls at or below the boundary count as IBD
    ibd = coal[coal["coal_time"] <= args.time_boundary].copy()
    ibd["segment_length"] = ibd["right"] - ibd["left"]

    ibd_length = ibd.groupby(["sample_1", "sample_2"])["segment_length"].sum()

    pairs = coal[["sample_1", "sample_2"]].drop_duplicates().reset_index(drop=True)
    pairs["frac_IBD"] = pairs.apply(
        lambda row: ibd_length.get((row["sample_1"], row["sample_2"]), 0.0) / GENOME_LENGTH,
        axis=1,
    )

    out_path = os.path.join(args.output_dir, f"relate_strain_{args.random_seed}.csv")
    pairs.to_csv(out_path, index=False)
    logging.info(f"### IBD fractions - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"relate_t_strain.py COMPLETED in {_format_elapsed(elapsed)}")


if __name__ == "__main__":
    main()
