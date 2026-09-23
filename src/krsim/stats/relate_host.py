"""Per-host average of the true pairwise strain IBD fraction."""

import itertools

import pandas as pd


def host_average_ibd(relate, ind_node_map, sim_id):
    """
    For each host in the most recent generation (time == 0) of simulation
    `sim_id`, the average frac_IBD over every pair of that host's strains.
    `relate` is the relate_strain table, `ind_node_map` the sim_ind_node_map
    table. Returns a DataFrame with columns host, avg_frac_IBD.
    """
    # --- node_id -> frac_IBD lookup, keyed by sorted (sample_1, sample_2) ---
    pair_to_frac = {
        (int(row.sample_1), int(row.sample_2)): row.frac_IBD
        for row in relate.itertuples()
    }

    # --- ind_id -> node_ids for this sim_id, restricted to the most recent generation (time=0) ---
    ind_node_map = ind_node_map[
        (ind_node_map["sim_id"] == sim_id) & (ind_node_map["time"] == 0)
    ]
    ind_to_nodes = ind_node_map.groupby("ind_id")["node_id"].apply(list).to_dict()

    # --- average pairwise IBD across each host's own strain nodes ---
    # pairs with no recorded coalescence anywhere in the ARG are assumed to have frac_IBD=0
    rows = []
    for ind_id, node_ids in sorted(ind_to_nodes.items()):
        fracs = [
            pair_to_frac.get((node_a, node_b), 0.0)
            for node_a, node_b in itertools.combinations(sorted(node_ids), 2)
        ]

        avg_frac_IBD = sum(fracs) / len(fracs) if fracs else float("nan")
        rows.append({"host": ind_id, "avg_frac_IBD": avg_frac_IBD})

    return pd.DataFrame(rows)
