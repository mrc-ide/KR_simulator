import os
import csv
import pandas as pd
from collections import defaultdict


def _build_recomb_to_nodes(ts):
    recomb_to_nodes = defaultdict(list)
    for node_id, node in enumerate(ts.tables.nodes):
        recomb_id = node.metadata.get("recomb_ID")
        if recomb_id is not None:
            recomb_to_nodes[recomb_id].append(node_id)
    return recomb_to_nodes


def get_node_pairs(ts):
    recomb_to_nodes = _build_recomb_to_nodes(ts)
    node_pairs = {}
    for node_id, node in enumerate(ts.tables.nodes):
        recomb_id = node.metadata.get("recomb_ID")
        if recomb_id is None:
            continue
        pair_recomb_id = recomb_id + 1 if recomb_id % 2 == 0 else recomb_id - 1
        node_pairs[node_id] = recomb_to_nodes.get(pair_recomb_id, [])
    return node_pairs


def save_node_pairs(ts, sim_id, output_dir="."):
    node_pairs = get_node_pairs(ts)

    rows = [
        {"sim_id": sim_id, "node_id": node_id, "paired_node_id": paired_node_id}
        for node_id, paired_nodes in sorted(node_pairs.items())
        for paired_node_id in paired_nodes
        if node_id < paired_node_id
    ]
    df = pd.DataFrame(rows, columns=["sim_id", "node_id", "paired_node_id"])

    output_path = os.path.join(output_dir, f"node_pairs_{sim_id}.csv")
    if os.path.exists(output_path):
        existing = pd.read_csv(output_path)
        existing = existing[existing["sim_id"] != sim_id]
        df = pd.concat([existing, df], ignore_index=True)

    df.to_csv(output_path, index=False)
    return df
