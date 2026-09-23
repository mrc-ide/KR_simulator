"""krsim coal — Stage 3: coalescence and mutation time tables."""

import os

HELP = "Stage 3: coalescence-time and mutation-time tables"

EPILOG = """\
Outputs:
  always  coal_T_<seed>.csv, mut_T_<seed>.csv

Example:
  krsim coal --ts results/mutated_42.trees --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--ts", type=str, required=True, metavar="FILE",
                          help="Mutated tree sequence from `krsim mutate` (mutated_<seed>.trees)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import tskit
    from krsim.output.coal_mut_times import mutation_table, coalescence_table

    os.makedirs(args.output_dir, exist_ok=True)

    # load tree sequence
    ts = tskit.load(args.ts)

    mut_csv_path = os.path.join(args.output_dir, f"mut_T_{args.random_seed}.csv")
    mutation_table(ts).to_csv(mut_csv_path, index=False)

    coal_csv_path = os.path.join(args.output_dir, f"coal_T_{args.random_seed}.csv")
    coalescence_table(ts).to_csv(coal_csv_path, index=False)
