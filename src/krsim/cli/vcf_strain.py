"""krsim vcf-strain — Stage 4: per-strain VCF."""

import os
import time
import logging

from krsim.cli._common import format_elapsed

HELP = "Stage 4: per-strain VCF"

EPILOG = """\
Outputs:
  always  output_<seed>_strain.vcf

Example:
  krsim vcf-strain --ts results/mutated_42.trees --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--ts", type=str, required=True, metavar="FILE",
                          help="Mutated tree sequence from `krsim mutate` (mutated_<seed>.trees)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Seed of the run (names the output files)")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Where to write the output files")


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
