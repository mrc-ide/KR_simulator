import os
import pandas as pd
from collections import defaultdict

def append_ind_node_map(ts, sim_id, output_dir="."):

    # build {ind_id: [(node_id, time), ...]} from nodes table
    ind_node_map = defaultdict(list)
    for node_id, node in enumerate(ts.tables.nodes):
        if node.individual != -1:
            ind_node_map[node.individual].append((node_id, node.time))

    # flatten into rows
    rows = []
    for ind_id, node_id_times in sorted(ind_node_map.items()):
        for node_id, time in node_id_times:
            rows.append({
                "sim_id":  sim_id,
                "ind_id":  ind_id,
                "node_id": node_id,
                "time":    time,
            })

    df = pd.DataFrame(rows)

    output_path = os.path.join(output_dir, f"sim_ind_node_map_{sim_id}.csv")
    if os.path.exists(output_path):
        existing = pd.read_csv(output_path)
        existing = existing[existing["sim_id"] != sim_id]
        df = pd.concat([existing, df], ignore_index=True)

    df.to_csv(output_path, index=False)
    return df