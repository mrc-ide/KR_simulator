"""krsim coal — Stage 3: coalescence and mutation time tables."""

import os

from krsim.cli._common import Run

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

    run = Run("coal", args.output_dir, args.random_seed)

    ts = tskit.load(args.ts)
    run.detail(f"read {args.ts} ({ts.num_samples} samples, {ts.num_trees} trees, {ts.num_mutations} mutations)")

    coal = coalescence_table(ts)
    coal_csv_path = os.path.join(args.output_dir, f"coal_T_{args.random_seed}.csv")
    coal.to_csv(coal_csv_path, index=False)
    run.wrote(coal_csv_path, rows=len(coal))
    n_pairs = len(coal[["sample_1", "sample_2"]].drop_duplicates())
    run.detail(f"  {len(coal):,} coalescence events (one per strain pair per tree); "
               f"{n_pairs} of {ts.num_samples * (ts.num_samples - 1) // 2} strain pairs coalesce")

    mut = mutation_table(ts)
    mut_csv_path = os.path.join(args.output_dir, f"mut_T_{args.random_seed}.csv")
    mut.to_csv(mut_csv_path, index=False)
    run.wrote(mut_csv_path, rows=len(mut))

    run.done()
