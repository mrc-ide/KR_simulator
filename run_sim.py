import os
import time
import logging
import argparse
from sim.params import SimParams
from sim.simulator import KRSimulator
from output.sim_ind_node_id import append_ind_node_map
from output.node_pairs import save_node_pairs
from output.lineage import save_lineage
from output.node_id_map import save_node_id_map


def parse_args():
    parser = argparse.ArgumentParser(
        description="KR Malaria Transmission Simulator"
    )

    # --- required parameters ---
    parser.add_argument("--Nh", type=int, required=True,
                        help="Number of human hosts")
    parser.add_argument("--T", type=int, required=True,
                        help="Number of generations")
    parser.add_argument("--K", type=int, required=True,
                        help="Number of parasite strains per transmission")
    parser.add_argument("--chi", type=float, required=True,
                        help="Superinfection probability (0.0 - 1.0)")
    parser.add_argument("--C", type=int, required=True,
                        help="Number of oocysts per mosquito")
    parser.add_argument("--recombination_rate", type=float, required=True,
                        help="Per bp recombination rate")
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility")

    parser.add_argument("--out_dir", type=str, required=True,
                        help="Directory to write output files to.")
    

    # --- optional parameters ---
    parser.add_argument("--save_ts", action="store_true",
                        help="If set, saves the tree sequence object to a .trees file.")
    parser.add_argument("--save_ID", action="store_true",
                        help="If set, saves a CSV of simulation, individual and node mapping across all generations.")
    parser.add_argument("--save_pairs", action="store_true",
                        help="If set, saves a CSV of recombinant node pairs across all generations.")
    parser.add_argument("--save_lineage", action="store_true",
                        help="If set, saves a CSV of per-generation child-parent strain lineage.")


    return parser.parse_args()


def _format_elapsed(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}m"
    else:
        return f"{seconds / 3600:.2f}h"


def main():
    args = parse_args()
    start_time = time.time()

    # create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    # set up logging: full detail (including DEBUG) to the log file, only
    # INFO-and-above to the terminal, so the terminal stays a clean list of
    # "what got created" while the file keeps the full record
    log_path = os.path.join(args.out_dir, f"KR_sim_{args.random_seed}.log")

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter("%(asctime)s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(logging.Formatter("%(message)s"))

    logging.basicConfig(level=logging.DEBUG, handlers=[file_handler, stream_handler])

    logging.debug(f"Simulation parameters: {vars(args)}")

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
    logging.debug(f"Simulation COMPLETED in {_format_elapsed(elapsed)}")

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


if __name__ == "__main__":
    main()