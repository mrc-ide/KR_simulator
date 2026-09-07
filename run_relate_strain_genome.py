import os
import time
import logging
import argparse
import pandas as pd

from run_relate_strain import _is_mirror_pair


def parse_args():
    parser = argparse.ArgumentParser(
        description="Genome-wide pairwise strain relatedness: per-interval generation of "
                     "coalescence (from coal_T), annotated with same-meiosis mirror pairs "
                     "(from the pedigree lineage)."
    )

    # --- required parameters ---
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility (also used as sim_id)")

    # --- optional parameters ---
    parser.add_argument("--coal_csv", type=str, default=None,
                        help="Path to coalescence-time CSV (default: <output_dir>/coal_T_<random_seed>.csv)")
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


def main():
    args = parse_args()
    start_time = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    coal_csv = args.coal_csv or os.path.join(args.output_dir, f"coal_T_{args.random_seed}.csv")
    lineage_csv = args.lineage_csv or os.path.join(args.output_dir, f"lineage_{args.random_seed}.csv")
    node_id_map_csv = args.node_id_map_csv or os.path.join(args.output_dir, f"node_id_map_{args.random_seed}.csv")

    coal = pd.read_csv(coal_csv)
    lineage = pd.read_csv(lineage_csv)
    node_map = pd.read_csv(node_id_map_csv).dropna(subset=["post_simplify_node_ID"])
    node_map["post_simplify_node_ID"] = node_map["post_simplify_node_ID"].astype(int)

    # recomb_ID / parent lookups per strain, in pre-simplify node ID space
    # (same lookups as run_relate_strain.py)
    recomb_of = lineage.drop_duplicates("child_node_ID").set_index("child_node_ID")["recomb_ID"].to_dict()
    parents_of = lineage.groupby("child_node_ID")["parent_node_ID"].unique().apply(list).to_dict()
    post_to_pre = dict(zip(node_map["post_simplify_node_ID"], node_map["pre_simplify_node_ID"]))

    # one mirror-pair label per (sample_1, sample_2) pair -- a static fact about
    # which meiosis produced them, not something that varies by genomic interval.
    # "Y"   = reciprocal recombination products of a genuine two-parent meiosis:
    #         they share 0% of the genome through that event, so at every
    #         position a per-locus trace (via coal_T's tree.mrca/tree.time)
    #         looks straight past it to whatever real ancestry the two distinct
    #         immediate parents may or may not share.
    # "Y_C" = reciprocal products of a clonal self-cross (both drawn from the
    #         same single parent strain): adjacent recomb_IDs, but there's only
    #         one immediate parent, so they still coalesce normally -- coal_T's
    #         real MRCA times already reflect this correctly, this label is
    #         purely informational.
    # "NA"  = not a mirror pair at all.
    pairs = coal[["sample_1", "sample_2"]].drop_duplicates()
    mirror_lookup = {}
    for sample_1, sample_2 in pairs.itertuples(index=False):
        pre_1 = post_to_pre.get(sample_1)
        pre_2 = post_to_pre.get(sample_2)
        if not _is_mirror_pair(recomb_of.get(pre_1), recomb_of.get(pre_2)):
            mirror_lookup[(sample_1, sample_2)] = "NA"
        else:
            distinct_parents = set(parents_of.get(pre_1, [])) | set(parents_of.get(pre_2, []))
            mirror_lookup[(sample_1, sample_2)] = "Y" if len(distinct_parents) > 1 else "Y_C"

    result = coal.copy()
    result.insert(0, "sim_id", args.random_seed)
    result = result.rename(columns={"coal_time": "n_gens_to_coalescence"})
    result["relatedness"] = 0.5 ** result["n_gens_to_coalescence"]
    result["mirror_pair"] = [
        mirror_lookup[(s1, s2)] for s1, s2 in zip(result["sample_1"], result["sample_2"])
    ]

    out_path = os.path.join(args.output_dir, f"pedigree_relate_strain_genome_{args.random_seed}.csv")
    result.to_csv(out_path, index=False)
    logging.info(f"### genome-wide pedigree relatedness - {out_path}")

    elapsed = time.time() - start_time
    logging.info(f"run_relate_strain_genome.py COMPLETED in {_format_elapsed(elapsed)}")


if __name__ == "__main__":
    main()
