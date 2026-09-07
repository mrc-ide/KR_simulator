import os
import time
import logging
import argparse
import tskit

from output.vcf_strain_single import process_vcf


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write per-strain VCF from a mutated tree sequence"
    )

    # --- required parameters ---
    parser.add_argument("--ts", type=str, required=True,
                        help="Path to mutated tree sequence (.trees) file")
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility")
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

    # load mutated tree sequence
    ts = tskit.load(args.ts)

    raw_vcf = os.path.join(args.output_dir, f"output_temp_{args.random_seed}.vcf")
    final_vcf = os.path.join(args.output_dir, f"output_{args.random_seed}_strain.vcf")

    with open(raw_vcf, "w") as f:
        ts.write_vcf(f)

    process_vcf(ts, input_vcf=raw_vcf, output_vcf=final_vcf)
    logging.info(f"### strain VCF - {final_vcf}")

    elapsed = time.time() - start_time
    logging.info(f"run_vcf_strain.py COMPLETED in {_format_elapsed(elapsed)}")


if __name__ == "__main__":
    main()
