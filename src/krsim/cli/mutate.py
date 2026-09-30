"""krsim mutate — Stage 2: add mutations to a tree sequence."""

import os

from krsim.cli._common import Run
from krsim.defaults import MUTATION_RATE

HELP = "Stage 2: add mutations to a tree sequence"

EPILOG = """\
Outputs:
  always  ts_mut_<seed>.trees

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

    run = Run("mutate", args.output_dir, args.random_seed)

    ts = tskit.load(args.ts)
    run.detail(f"read {args.ts} ({ts.num_samples} samples, {ts.num_trees} trees)")

    ts = msprime.sim_mutations(ts, rate=args.mutation_rate, random_seed=args.random_seed)

    out_path = os.path.join(args.output_dir, f"ts_mut_{args.random_seed}.trees")
    ts.dump(out_path)
    run.wrote(out_path)
    run.detail(f"  {ts.num_mutations} mutations at {ts.num_sites} sites across {ts.num_trees} trees "
               f"(rate {args.mutation_rate})")

    run.done()
