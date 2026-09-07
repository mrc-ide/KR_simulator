import os
import time
import logging
import argparse
import itertools
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute pairwise strain relatedness by tracing the pedigree back to a common ancestor"
    )

    # --- required parameters ---
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility (also used as sim_id)")

    # --- optional parameters ---
    parser.add_argument("--lineage_csv", type=str, default=None,
                        help="Path to lineage CSV (default: <output_dir>/lineage_<random_seed>.csv)")
    parser.add_argument("--node_id_map_csv", type=str, default=None,
                        help="Path to node ID map CSV (default: <output_dir>/node_id_map_<random_seed>.csv)")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Directory to write output files to (default: current directory).")

    return parser.parse_args()


def _format_elapsed(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}m"
    else:
        return f"{seconds / 3600:.2f}h"


def _is_mirror_pair(recomb_a, recomb_b):
    """
    Two recomb IDs are reciprocal products of the same crossover event if
    they are adjacent with the lower one even (e.g. 0 & 1, but not 1 & 2).
    Reciprocal products share 0% of the genome by descent through that
    event, even though they came from the same two immediate parents.
    """
    if pd.isna(recomb_a) or pd.isna(recomb_b):
        return False
    lower, upper = sorted((int(recomb_a), int(recomb_b)))
    return lower % 2 == 0 and upper == lower + 1


def compute_relatedness(pre_a, pre_b, parents_of, recomb_of, max_steps):
    """
    Trace two strains' pedigrees back generation by generation until they
    share a common ancestor node, or run out of ancestry.

    Parameters
    ----------
    pre_a, pre_b : int
        Pre-simplify node IDs of the two strains being compared.
    parents_of : dict[int, list[int]]
        Pre-simplify node ID -> list of its parent node ID(s).
    recomb_of : dict[int, int]
        Pre-simplify node ID -> the recomb_ID it was produced by.
    max_steps : int
        Upper bound on how many generations to climb (safety cap).

    Returns
    -------
    relatedness : float
        0.5 ** n, where n is the number of generations back to the first
        shared ancestor; 0.0 if the pair never coalesces, or if they are
        reciprocal recombination products of the same meiosis.
    n_gens : int or None
        Number of generations back to coalescence, 0 for the reciprocal
        special case, or None if the pair never coalesces.
    is_mirror : bool
        True if pre_a and pre_b are reciprocal recombination products of
        the same meiosis (relatedness is forced to 0.0 without climbing).
    """
    if _is_mirror_pair(recomb_of.get(pre_a), recomb_of.get(pre_b)):
        return 0.0, 0, True

    frontier_a = {pre_a}
    frontier_b = {pre_b}

    for n in range(1, max_steps + 1):
        frontier_a = {p for node in frontier_a for p in parents_of.get(node, [])}
        frontier_b = {p for node in frontier_b for p in parents_of.get(node, [])}

        if not frontier_a or not frontier_b:
            return 0.0, None, False

        if frontier_a & frontier_b:
            return 0.5 ** n, n, False

    return 0.0, None, False


def main():
    args = parse_args()
    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    lineage_csv = args.lineage_csv or os.path.join(args.output_dir, f"lineage_{args.random_seed}.csv")
    node_id_map_csv = args.node_id_map_csv or os.path.join(args.output_dir, f"node_id_map_{args.random_seed}.csv")

    lineage = pd.read_csv(lineage_csv)
    node_map = pd.read_csv(node_id_map_csv).dropna(subset=["post_simplify_node_ID"])
    node_map["post_simplify_node_ID"] = node_map["post_simplify_node_ID"].astype(int)

    # pedigree lookups, in pre-simplify node ID space
    parents_of = lineage.groupby("child_node_ID")["parent_node_ID"].unique().apply(list).to_dict()
    recomb_of = lineage.drop_duplicates("child_node_ID").set_index("child_node_ID")["recomb_ID"].to_dict()

    # restrict to the final generation (gen == 0); node_id_map also retains
    # ancestral root nodes (kept via keep_input_roots=True) which aren't
    # actual strains to compare
    pre_to_post = dict(zip(node_map["pre_simplify_node_ID"], node_map["post_simplify_node_ID"]))
    final_gen_pre = lineage.loc[lineage["gen"] == 0, "child_node_ID"].unique()

    post_to_pre = {pre_to_post[pre_id]: pre_id for pre_id in final_gen_pre}
    samples = sorted(post_to_pre)

    max_steps = int(lineage["gen"].max()) + 1

    rows = []
    for sample_1, sample_2 in itertools.combinations(samples, 2):
        relatedness, n_gens, is_mirror = compute_relatedness(
            post_to_pre[sample_1], post_to_pre[sample_2], parents_of, recomb_of, max_steps
        )
        rows.append({
            "sim_id": args.random_seed,
            "sample_1": sample_1,
            "sample_2": sample_2,
            "relatedness": relatedness,
            "n_gens_to_coalescence": n_gens,
            "mirror_pair": "Y" if is_mirror else "NA",
        })

    result = pd.DataFrame(rows)

    out_path = os.path.join(args.output_dir, f"pedigree_relate_strain_{args.random_seed}.csv")
    result.to_csv(out_path, index=False, na_rep="NA")
    logging.info(f"### pedigree relatedness - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"run_relate_strain.py COMPLETED in {_format_elapsed(elapsed)}")


if __name__ == "__main__":
    main()
