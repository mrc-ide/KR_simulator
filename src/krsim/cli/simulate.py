"""krsim simulate — Stage 1: simulate the genealogy."""

import os

from krsim.cli._common import Run
from krsim.defaults import RECOMBINATION_RATE

HELP = "Stage 1: simulate the transmission genealogy"

EPILOG = """\
Outputs:
  always          ts_<seed>.trees, sim_ind_node_map_<seed>.csv, KR_sim_<seed>.log
  --save_lineage  lineage_<seed>.csv, node_id_map_<seed>.csv
  --save_pairs    node_pairs_<seed>.csv

Example:
  krsim simulate --Nh 10 --T 50 --K 2 --chi 0.18 --C 1 \\
      --random_seed 42 --output_dir results/ --save_lineage
"""



def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--Nh", type=int, required=True, metavar="N",
                          help="Number of human hosts")
    required.add_argument("--T", type=int, required=True, metavar="N",
                          help="Number of generations")
    required.add_argument("--K", type=int, required=True, metavar="N",
                          help="Parasites transmitted per mosquito bite")
    required.add_argument("--chi", type=float, required=True, metavar="P",
                          help="Superinfection probability (0-1)")
    required.add_argument("--C", type=int, required=True, metavar="N",
                          help="Oocysts per mosquito")
    # model parameter, but falls back to the P. falciparum value if not given
    required.add_argument("--recombination_rate", type=float, default=RECOMBINATION_RATE, metavar="R",
                          help=f"Per bp recombination rate (default: {RECOMBINATION_RATE})")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")

    parser.add_argument("--save_lineage", action="store_true",
                        help="Also save the true pedigree (for the visualizations)")
    parser.add_argument("--save_pairs", action="store_true",
                        help="Also save same-meiosis strain pairs")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    from krsim.sim.params import SimParams
    from krsim.sim.simulator import KRSimulator
    from krsim.output.sim_ind_node_id import append_ind_node_map
    from krsim.output.node_pairs import save_node_pairs
    from krsim.output.lineage import save_lineage
    from krsim.output.node_id_map import save_node_id_map

    run = Run("simulate", args.output_dir, args.random_seed)

    params = SimParams(
        Nh=args.Nh,
        T=args.T,
        K=args.K,
        chi=args.chi,
        C=args.C,
        recombination_rate=args.recombination_rate,
        random_seed=args.random_seed
    )
    run.detail(f"parameters: Nh={params.Nh} T={params.T} K={params.K} chi={params.chi} C={params.C} "
               f"recombination_rate={params.recombination_rate} random_seed={params.random_seed}")

    ts, raw_tables, node_map = KRSimulator(params).run(capture_lineage=args.save_lineage)
    run.detail(f"simulated {params.T} generations -> {ts.num_nodes} nodes after simplify")

    # always save the tree sequence: every downstream stage starts from it
    trees_path = os.path.join(args.output_dir, f"ts_{args.random_seed}.trees")
    ts.dump(trees_path)
    run.wrote(trees_path)

    # always save the individual-node mapping: vcf-host and relate-host need it
    df = append_ind_node_map(ts, sim_id=args.random_seed, output_dir=args.output_dir)
    run.wrote(os.path.join(args.output_dir, f"sim_ind_node_map_{args.random_seed}.csv"), rows=len(df))

    if args.save_lineage:
        df = save_lineage(raw_tables, sim_id=args.random_seed, output_dir=args.output_dir)
        run.wrote(os.path.join(args.output_dir, f"lineage_{args.random_seed}.csv"), rows=len(df))
        df = save_node_id_map(node_map, sim_id=args.random_seed, output_dir=args.output_dir)
        run.wrote(os.path.join(args.output_dir, f"node_id_map_{args.random_seed}.csv"), rows=len(df))

    if args.save_pairs:
        df = save_node_pairs(ts, sim_id=args.random_seed, output_dir=args.output_dir)
        run.wrote(os.path.join(args.output_dir, f"node_pairs_{args.random_seed}.csv"), rows=len(df))

    run.done()
