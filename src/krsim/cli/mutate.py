"""krsim mutate — Stage 2: add mutations to a tree sequence."""

import os
import time
import logging

from krsim.cli._common import format_elapsed, setup_run_logging

HELP = "Stage 2: add mutations to a tree sequence"


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    # --- required parameters ---
    required.add_argument("--ts", type=str, required=True,
                          help="Path to tree sequence (.trees) file")
    required.add_argument("--mutation_rate", type=float, required=True,
                          help="Per bp mutation rate (e.g. 1e-6)")
    required.add_argument("--random_seed", type=int, required=True,
                          help="Random seed for reproducibility")

    required.add_argument("--output_dir", type=str, required=True,
                          help="Directory to write output files.")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    import tskit
    import msprime

    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    setup_run_logging(args.output_dir, args.random_seed)

    # load tree sequence
    ts = tskit.load(args.ts)

    # add mutations
    ts = msprime.sim_mutations(ts, rate=args.mutation_rate, random_seed=args.random_seed)
    logging.debug(f"mutations added at rate {args.mutation_rate}")

    # save mutated tree sequence
    out_path = os.path.join(args.output_dir, f"mutated_{args.random_seed}.trees")
    ts.dump(out_path)
    logging.info(f"### mutated tree sequence - {os.path.basename(out_path)}")

    elapsed = time.time() - start_time
    logging.debug(f"krsim mutate COMPLETED in {format_elapsed(elapsed)}")
