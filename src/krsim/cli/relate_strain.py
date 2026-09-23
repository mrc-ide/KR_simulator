"""krsim relate-strain — true pairwise IBD fraction between strains."""

import os

from krsim.cli._common import Run

HELP = "True pairwise IBD fraction between strains"

EPILOG = """\
Outputs:
  always  relate_strain_<seed>.csv

Example:
  krsim relate-strain --coal_csv results/coal_T_42.csv --time_boundary 2 \\
      --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--coal_csv", type=str, required=True, metavar="FILE",
                          help="Coalescence table from `krsim coal` (coal_T_<seed>.csv)")
    required.add_argument("--time_boundary", type=float, required=True, metavar="T",
                          help="Segments that coalesce within T generations count as IBD")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import pandas as pd
    from krsim.stats.relate_strain import strain_ibd_fractions

    run = Run("relate-strain", args.output_dir, args.random_seed)

    coal = pd.read_csv(args.coal_csv)
    run.detail(f"read {args.coal_csv} ({len(coal):,} rows); IBD if coal_time <= {args.time_boundary}")

    pairs = strain_ibd_fractions(coal, args.time_boundary)
    out_path = os.path.join(args.output_dir, f"relate_strain_{args.random_seed}.csv")
    pairs.to_csv(out_path, index=False)
    run.wrote(out_path, rows=len(pairs))

    run.done()
