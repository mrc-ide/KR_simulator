"""krsim coal — Stage 3: coalescence and mutation time tables."""

import os

HELP = "Stage 3: write the coalescence-time (coal_T) and mutation-time (mut_T) tables"


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    # --- required parameters ---
    required.add_argument("--ts", type=str, required=True,
                          help="Path to tree sequence (.trees) file")
    required.add_argument("--random_seed", type=int, required=True,
                          help="Random seed for reproducibility")
    required.add_argument("--output_dir", type=str, required=True,
                          help="Directory to write output files.")


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
