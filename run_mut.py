import os
import time
import logging
import argparse
import tskit
import msprime


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add mutations to a tree sequence"
    )

    # --- required parameters ---
    parser.add_argument("--ts", type=str, required=True,
                        help="Path to tree sequence (.trees) file")
    parser.add_argument("--mutation_rate", type=float, required=True,
                        help="Per bp mutation rate (e.g. 1e-6)")
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility")

    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to write output files.")

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

    log_path = os.path.join(args.output_dir, f"KR_sim_{args.random_seed}.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(),
        ]
    )

    logging.info("")
    logging.info("========================================")
    logging.info("========== RUNNING run_mut.py ==========")
    logging.info("========================================")

    # load tree sequence
    ts = tskit.load(args.ts)
    # logging.info(f"loaded tree sequence from {args.ts}")

    # add mutations
    ts = msprime.sim_mutations(ts, rate=args.mutation_rate, random_seed=args.random_seed)
    logging.info(f"mutations added at rate {args.mutation_rate}")

    # save mutated tree sequence
    out_path = os.path.join(args.output_dir, f"mutated_{args.random_seed}.trees")
    ts.dump(out_path)
    logging.info(f"### mutated tree sequence - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"run_mut.py COMPLETED in {_format_elapsed(elapsed)}")


if __name__ == "__main__":
    main()
