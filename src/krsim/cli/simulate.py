"""krsim simulate — Stage 1: simulate the genealogy."""

import os
import time
import logging

from krsim.cli._common import format_elapsed, setup_run_logging
from krsim.defaults import RECOMBINATION_RATE

HELP = "Stage 1: simulate the transmission genealogy (tree sequence)"


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    # --- required parameters ---
    required.add_argument("--Nh", type=int, required=True,
                          help="Number of human hosts")
    required.add_argument("--T", type=int, required=True,
                          help="Number of generations")
    required.add_argument("--K", type=int, required=True,
                          help="Number of parasite strains per transmission")
    required.add_argument("--chi", type=float, required=True,
                          help="Superinfection probability (0.0 - 1.0)")
    required.add_argument("--C", type=int, required=True,
                          help="Number of oocysts per mosquito")
    required.add_argument("--random_seed", type=int, required=True,
                          help="Random seed for reproducibility")

    required.add_argument("--out_dir", type=str, required=True,
                          help="Directory to write output files to.")

    # --- optional parameters ---
    parser.add_argument("--recombination_rate", type=float, default=RECOMBINATION_RATE,
                        help=f"Per bp recombination rate (default: {RECOMBINATION_RATE})")
    parser.add_argument("--save_ts", action="store_true",
                        help="If set, saves the tree sequence object to a .trees file.")
    parser.add_argument("--save_ID", action="store_true",
                        help="If set, saves a CSV of simulation, individual and node mapping across all generations.")
    parser.add_argument("--save_pairs", action="store_true",
                        help="If set, saves a CSV of recombinant node pairs across all generations.")
    parser.add_argument("--save_lineage", action="store_true",
                        help="If set, saves a CSV of per-generation child-parent strain lineage.")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/tskit/msprime)
    from krsim.sim.params import SimParams
    from krsim.sim.simulator import KRSimulator
    from krsim.output.sim_ind_node_id import append_ind_node_map
    from krsim.output.node_pairs import save_node_pairs
    from krsim.output.lineage import save_lineage
    from krsim.output.node_id_map import save_node_id_map

    start_time = time.time()

    # create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    setup_run_logging(args.out_dir, args.random_seed)
    # drop the entries argparse adds for subcommand dispatch
    sim_args = {k: v for k, v in vars(args).items() if k not in ("command", "func")}
    logging.debug(f"Simulation parameters: {sim_args}")

    # run simulation
    params = SimParams(
        Nh=args.Nh,
        T=args.T,
        K=args.K,
        chi=args.chi,
        C=args.C,
        recombination_rate=args.recombination_rate,
        random_seed=args.random_seed
    )

    sim = KRSimulator(params)
    ts, raw_tables, node_map = sim.run(capture_lineage=args.save_lineage)

    elapsed = time.time() - start_time
    logging.debug(f"Simulation COMPLETED in {format_elapsed(elapsed)}")

    ### OUPUT FILES ###

    # save tree sequence
    if args.save_ts:
        trees_path = os.path.join(args.out_dir, f"ts_{args.random_seed}.trees")
        ts.dump(trees_path)
        logging.info(f"### tree sequence - {os.path.basename(trees_path)}")

    # save individual-node mapping
    if args.save_ID:
        append_ind_node_map(ts, sim_id=args.random_seed, output_dir=args.out_dir)
        logging.info(f"### individual-node mapping - sim_ind_node_map_{args.random_seed}.csv")

    # save node pairs
    if args.save_pairs:
        save_node_pairs(ts, sim_id=args.random_seed, output_dir=args.out_dir)
        logging.info(f"### recombinant node pairs - node_pairs_{args.random_seed}.csv")

    # save per-generation child-parent lineage
    if args.save_lineage:
        save_lineage(raw_tables, sim_id=args.random_seed, output_dir=args.out_dir)
        logging.info(f"### strain lineage - lineage_{args.random_seed}.csv")
        save_node_id_map(node_map, sim_id=args.random_seed, output_dir=args.out_dir)
        logging.info(f"### node ID map (pre/post simplify) - node_id_map_{args.random_seed}.csv")
