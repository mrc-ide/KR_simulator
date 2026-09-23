import os
import numpy as np
import pandas as pd
import tskit


def save_node_id_map(node_map, sim_id, output_dir="."):
    """
    Save the mapping from pre-simplify node ID to post-simplify node ID.

    Parameters
    ----------
    node_map : numpy.ndarray
        Array returned by TableCollection.simplify(), where node_map[i] is
        the ID of raw node i in the final tree sequence, or tskit.NULL (-1)
        if that node did not survive simplification.
    sim_id : int
        Simulation / random seed identifier.
    output_dir : str
        Directory to write the CSV to.
    """
    df = pd.DataFrame({
        "sim_id": sim_id,
        "pre_simplify_node_ID": np.arange(len(node_map)),
        "post_simplify_node_ID": node_map,
    })
    df.loc[df["post_simplify_node_ID"] == tskit.NULL, "post_simplify_node_ID"] = pd.NA

    output_path = os.path.join(output_dir, f"node_id_map_{sim_id}.csv")
    if os.path.exists(output_path):
        existing = pd.read_csv(output_path)
        existing = existing[existing["sim_id"] != sim_id]
        df = pd.concat([existing, df], ignore_index=True)

    df.to_csv(output_path, index=False)
    return df
