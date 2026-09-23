"""krsim vcf-host — Stage 5: per-host VCF."""

import os
import time
import logging

from krsim.cli._common import format_elapsed

HELP = "Stage 5: per-host VCF (mixes each host's strains into reads)"

EPILOG = """\
Outputs:
  always              output_<seed>_host.vcf
  --save_proportions  host_proportions_<seed>.csv

Example:
  krsim vcf-host --strain_vcf results/output_42_strain.vcf \\
      --ind_node_map_csv results/sim_ind_node_map_42.csv \\
      --random_seed 42 --depth 10 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--strain_vcf", type=str, required=True, metavar="FILE",
                          help="Strain VCF from `krsim vcf-strain`")
    required.add_argument("--ind_node_map_csv", type=str, required=True, metavar="FILE",
                          help="Host-strain map from `krsim simulate` (sim_ind_node_map_<seed>.csv)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--depth", type=int, required=True, metavar="N",
                          help="Reads drawn per host per site")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")

    parser.add_argument("--save_proportions", action="store_true",
                        help="Also save per-host strain proportions and read counts")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    from krsim.output.vcf_host_single import process_host_vcf

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    host_vcf = os.path.join(args.output_dir, f"output_{args.random_seed}_host.vcf")
    proportions_out = None
    if args.save_proportions:
        proportions_out = os.path.join(args.output_dir, f"host_proportions_{args.random_seed}.csv")

    process_host_vcf(
        sim_id=args.random_seed,
        input_vcf=args.strain_vcf,
        ind_node_map_csv=args.ind_node_map_csv,
        output_vcf=host_vcf,
        depth=args.depth,
        random_seed=args.random_seed,
        proportions_out=proportions_out,
    )
    logging.info(f"### host VCF - {host_vcf}")
    if proportions_out:
        logging.info(f"### host proportions - {proportions_out}")

    elapsed = time.time() - start_time
    logging.info(f"krsim vcf-host COMPLETED in {format_elapsed(elapsed)}")
