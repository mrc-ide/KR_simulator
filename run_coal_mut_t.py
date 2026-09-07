import os
import time
import logging
import argparse
import itertools
import numpy as np
import pandas as pd
import tskit
import msprime

from sim.params import COMP_BREAKPOINTS, GENOME_LENGTH

# start/end (bp) of each chromosome in the concatenated genome coordinate system
CHR_EDGES = np.concatenate(([0], COMP_BREAKPOINTS, [GENOME_LENGTH]))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Output information on coalescence and mutation time"
    )

    # --- required parameters ---
    parser.add_argument("--ts", type=str, required=True,
                        help="Path to tree sequence (.trees) file")
    # parser.add_argument("--pos", type=float, required=True,
    #                     help="Information across a specific position in the genome otherwise defaults to whole genome")
    # parser.add_argument("--mut_info", type=bool, action="store_true",
    #                     help="If set, outputs information on mutation times and positions")
    
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed for reproducibility")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to write output files.")

    return parser.parse_args()


def _format_elapsed(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}m"
    else:
        return f"{seconds / 3600:.2f}h"



args = parse_args()
os.makedirs(args.output_dir, exist_ok=True)


############ MUTATIONS INFORMATION TABLE ############

# load tree sequence
ts = tskit.load(args.ts)

sites = ts.tables.sites
mutations = ts.tables.mutations

genome_position = sites[mutations.site].position
# chromosome number (1-indexed) each position falls in
chromosome = np.searchsorted(CHR_EDGES, genome_position, side="right")
chr_position = genome_position - CHR_EDGES[chromosome - 1]

mutation_summary = pd.DataFrame({
    "node_id": mutations.node,
    "chromosome": chromosome,
    "position": chr_position,
    "time": mutations.time,
})

mut_csv_path = os.path.join(args.output_dir, f"mut_T_{args.random_seed}.csv")
mutation_summary.to_csv(mut_csv_path, index=False)


############ COALESCENCE INFORMATION TABLE ############

samples = ts.samples()
sample_pairs = list(itertools.combinations(samples, 2))
rows = []

for tree in ts.trees():
    interval_left, interval_right = tree.interval
    # chromosome number (1-indexed) the tree interval falls in
    chromosome = np.searchsorted(CHR_EDGES, interval_left, side="right")
    chr_start = CHR_EDGES[chromosome - 1]

    for sample_a, sample_b in sample_pairs:
        mrca = tree.mrca(sample_a, sample_b)

        # tskit.NULL means the pair has no common ancestor in this tree; skip it
        if mrca == tskit.NULL:
            continue

        rows.append({
            "sample_1": sample_a,
            "sample_2": sample_b,
            "ancestor": mrca,
            "chromosome": chromosome,
            "left": interval_left - chr_start,
            "right": interval_right - chr_start,
            "coal_time": ts.node(mrca).time,
        })

pair_coalescence_times = pd.DataFrame(rows)

coal_csv_path = os.path.join(args.output_dir, f"coal_T_{args.random_seed}.csv")
pair_coalescence_times.to_csv(coal_csv_path, index=False)


