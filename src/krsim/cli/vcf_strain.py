"""krsim vcf-strain — Stage 3: per-strain VCF."""

import os

from krsim.cli._common import Run

HELP = "Stage 3: per-strain VCF"

EPILOG = """\
Outputs:
  always  output_<seed>_strain.vcf

Example:
  krsim vcf-strain --ts results/ts_mut_42.trees --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--ts", type=str, required=True, metavar="FILE",
                          help="Mutated tree sequence from `krsim mutate` (ts_mut_<seed>.trees)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import tskit
    from krsim.output.vcf_strain_single import process_vcf

    run = Run("vcf-strain", args.output_dir, args.random_seed)

    ts = tskit.load(args.ts)
    run.detail(f"read {args.ts} ({ts.num_samples} samples, {ts.num_sites} variant sites)")

    raw_vcf = os.path.join(args.output_dir, f"output_temp_{args.random_seed}.vcf")
    final_vcf = os.path.join(args.output_dir, f"output_{args.random_seed}_strain.vcf")

    with open(raw_vcf, "w") as f:
        ts.write_vcf(f)

    process_vcf(ts, input_vcf=raw_vcf, output_vcf=final_vcf)
    run.wrote(final_vcf)

    run.done()
