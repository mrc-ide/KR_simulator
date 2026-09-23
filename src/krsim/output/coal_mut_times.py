"""Mutation-time and pairwise coalescence-time tables from a tree sequence."""

import numpy as np
import pandas as pd

from krsim.sim.params import COMP_BREAKPOINTS, GENOME_LENGTH

# start/end (bp) of each chromosome in the concatenated genome coordinate system
CHR_EDGES = np.concatenate(([0], COMP_BREAKPOINTS, [GENOME_LENGTH]))


def mutation_table(ts):
    """
    One row per mutation: the node it sits on, its chromosome, its
    chromosome-local position and its time.
    """
    sites = ts.tables.sites
    mutations = ts.tables.mutations

    genome_position = sites[mutations.site].position
    # chromosome number (1-indexed) each position falls in
    chromosome = np.searchsorted(CHR_EDGES, genome_position, side="right")
    chr_position = genome_position - CHR_EDGES[chromosome - 1]

    return pd.DataFrame({
        "node_id": mutations.node,
        "chromosome": chromosome,
        "position": chr_position,
        "time": mutations.time,
    })


def coalescence_table(ts):
    """
    One row per (sample pair, tree interval) where the pair has a common
    ancestor: the MRCA node, its time, and the interval's chromosome-local
    coordinates.
    """
    # For every tree, find each pair's MRCA and its time. A naive implementation
    # calls tree.mrca(a, b) for every one of the O(n^2) sample pairs in every
    # tree; at a few hundred samples that's `100k`s of Python-level calls per
    # tree, and collecting results as a list of per-row dicts (rather than
    # columnar arrays) multiplies memory several-fold on top of that — together
    # they OOM well before Nh reaches the low hundreds.
    #
    # Instead, for each tree we do a single postorder pass (children before
    # parents, handles multiple roots/incomplete coalescence and polytomies)
    # tracking which samples fall under each node so far. Two samples' MRCA is
    # exactly the node where their subtrees first merge, so the new pairs
    # coalescing at a node are the cross product of its not-yet-merged children's
    # sample sets — every pair is generated exactly once, in vectorized numpy
    # array ops rather than one Python call each.
    samples = ts.samples()
    node_times = ts.tables.nodes.time
    sample_flags = np.zeros(ts.num_nodes, dtype=bool)
    sample_flags[samples] = True

    s1_chunks, s2_chunks, anc_chunks = [], [], []
    chr_chunks, left_chunks, right_chunks, time_chunks = [], [], [], []

    for tree in ts.trees():
        interval_left, interval_right = tree.interval
        # chromosome number (1-indexed) the tree interval falls in
        chromosome = np.searchsorted(CHR_EDGES, interval_left, side="right")
        chr_start = CHR_EDGES[chromosome - 1]

        under = {}  # node -> array of sample ids under it, not yet consumed by a parent
        tree_s1, tree_s2, tree_anc = [], [], []
        for u in tree.postorder():
            if sample_flags[u]:
                under[u] = np.array([u], dtype=np.int64)
                continue
            child_arrays = [under.pop(c) for c in tree.children(u) if c in under]
            # every pair of samples drawn from two *different* children of u
            # coalesces here, for the first time, at u
            for i in range(len(child_arrays)):
                for j in range(i + 1, len(child_arrays)):
                    a, b = child_arrays[i], child_arrays[j]
                    aa = np.repeat(a, b.size)
                    bb = np.tile(b, a.size)
                    tree_s1.append(np.minimum(aa, bb))
                    tree_s2.append(np.maximum(aa, bb))
                    tree_anc.append(np.full(aa.size, u, dtype=np.int64))
            if child_arrays:
                under[u] = np.concatenate(child_arrays)

        if not tree_s1:
            continue

        s1 = np.concatenate(tree_s1)
        s2 = np.concatenate(tree_s2)
        anc = np.concatenate(tree_anc)
        n = s1.size

        s1_chunks.append(s1)
        s2_chunks.append(s2)
        anc_chunks.append(anc)
        chr_chunks.append(np.full(n, chromosome, dtype=np.int64))
        left_chunks.append(np.full(n, interval_left - chr_start))
        right_chunks.append(np.full(n, interval_right - chr_start))
        time_chunks.append(node_times[anc])

    pair_coalescence_times = pd.DataFrame({
        "sample_1": np.concatenate(s1_chunks),
        "sample_2": np.concatenate(s2_chunks),
        "ancestor": np.concatenate(anc_chunks),
        "chromosome": np.concatenate(chr_chunks),
        "left": np.concatenate(left_chunks),
        "right": np.concatenate(right_chunks),
        "coal_time": np.concatenate(time_chunks),
    })
    # match the original per-tree row order (sample_1 asc, then sample_2 asc)
    return pair_coalescence_times.sort_values(
        ["chromosome", "left", "sample_1", "sample_2"], kind="stable"
    ).reset_index(drop=True)
