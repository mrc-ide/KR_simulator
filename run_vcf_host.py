import os
import time
import logging
import argparse

from output.vcf_host_single import process_host_vcf


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate a per-strain VCF to per-host (individual) level"
    )

    # --- required parameters ---
    parser.add_argument("--strain_vcf", type=str, required=True,
                        help="Path to the per-strain VCF file (output of run_vcf_strain.py)")
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility (also used as sim_id)")
    parser.add_argument("--depth", type=int, required=True,
                        help="Total reads drawn per host per site (Multinomial depth)")

    # --- optional parameters ---
    parser.add_argument("--ind_node_map_csv", type=str, default=None,
                        help="Path to sim_ind_node_map CSV (default: <output_dir>/sim_ind_node_map_<random_seed>.csv)")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Directory to write output files to (default: current directory).")
    parser.add_argument("--save_proportions", action="store_true",
                        help="If set, saves a CSV of per-host Dirichlet strain proportions and read counts.")

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

    ind_node_map_csv = args.ind_node_map_csv
    if ind_node_map_csv is None:
        ind_node_map_csv = os.path.join(args.output_dir, f"sim_ind_node_map_{args.random_seed}.csv")

    host_vcf = os.path.join(args.output_dir, f"output_{args.random_seed}_host.vcf")
    proportions_out = None
    if args.save_proportions:
        proportions_out = os.path.join(args.output_dir, f"host_proportions_{args.random_seed}.csv")

    process_host_vcf(
        sim_id=args.random_seed,
        input_vcf=args.strain_vcf,
        ind_node_map_csv=ind_node_map_csv,
        output_vcf=host_vcf,
        depth=args.depth,
        random_seed=args.random_seed,
        proportions_out=proportions_out,
    )
    logging.info(f"### host VCF - {host_vcf}")
    if proportions_out:
        logging.info(f"### host proportions - {proportions_out}")

    elapsed = time.time() - start_time
    logging.info(f"run_vcf_host.py COMPLETED in {_format_elapsed(elapsed)}")


if __name__ == "__main__":
    main()
