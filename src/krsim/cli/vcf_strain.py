"""krsim vcf-strain — Stage 4: per-strain VCF."""

import os
import time
import logging

from krsim.cli._common import format_elapsed

HELP = "Stage 4: write a per-strain VCF from a mutated tree sequence"


def add_arguments(parser):
    # --- required parameters ---
    parser.add_argument("--ts", type=str, required=True,
                        help="Path to mutated tree sequence (.trees) file")
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Directory to write output files to (default: current directory).")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import tskit
    from krsim.output.vcf_strain_single import process_vcf

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
    logging.info(f"krsim vcf-strain COMPLETED in {format_elapsed(elapsed)}")
