"""krsim mutate — Stage 2: add mutations to a tree sequence."""

import os
import time
import logging

from krsim.cli._common import format_elapsed, setup_run_logging
from krsim.defaults import MUTATION_RATE

HELP = "Stage 2: add mutations to a tree sequence"

EPILOG = """\
Outputs:
  always  mutated_<seed>.trees

Example:
  krsim mutate --ts results/ts_42.trees --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--ts", type=str, required=True, metavar="FILE",
                          help="Tree sequence from `krsim simulate` (ts_<seed>.trees)")
    # model parameter, but falls back to the P. falciparum value if not given
    required.add_argument("--mutation_rate", type=float, default=MUTATION_RATE, metavar="R",
                          help=f"Per bp mutation rate (default: {MUTATION_RATE})")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")


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
