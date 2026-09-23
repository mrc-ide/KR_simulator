import os
import logging
import pandas as pd


def build_lineage(tables):
    """
    Build a table of child strain, the recombination event it was produced
    by, and the parent strain(s) it inherited genome segments from.

    A child node can have more than one parent node if its recombinant
    product carries segments from both parents of that transmission event,
    so each unique (child, parent) pair is written as its own row.

    Parameters
    ----------
    tables : tskit.TableCollection
        Tables to extract lineage information from. Must be the raw
        (pre-simplify) tables to get one row per real transmission event;
        a simplified tree sequence's tables can skip generations.

    Returns
    -------
    pd.DataFrame
        Columns: gen, child_node_ID, child_ind_ID, recomb_ID, parent_node_ID, parent_ind_ID.
    """
    nodes = tables.nodes
    edges = tables.edges

    lineage = pd.DataFrame({
        "child_node_ID": edges.child,
        "parent_node_ID": edges.parent,
    }).drop_duplicates().reset_index(drop=True)

    lineage["gen"] = nodes.time[lineage["child_node_ID"]]
    lineage["child_ind_ID"] = nodes.individual[lineage["child_node_ID"]]
    lineage["parent_ind_ID"] = nodes.individual[lineage["parent_node_ID"]]
    lineage["recomb_ID"] = [
        nodes[child_id].metadata.get("recomb_ID")
        for child_id in lineage["child_node_ID"]
    ]

    lineage = lineage[[
        "gen", "child_node_ID", "child_ind_ID", "recomb_ID", "parent_node_ID", "parent_ind_ID",
    ]]
    lineage = lineage.sort_values(["gen", "child_node_ID", "parent_node_ID"]).reset_index(drop=True)

    return lineage


def save_lineage(tables, sim_id, output_dir="."):
    lineage = build_lineage(tables)
    lineage.insert(0, "sim_id", sim_id)

    output_path = os.path.join(output_dir, f"lineage_{sim_id}.csv")
    if os.path.exists(output_path):
        existing = pd.read_csv(output_path)
        existing = existing[existing["sim_id"] != sim_id]
        lineage = pd.concat([existing, lineage], ignore_index=True)

    lineage.to_csv(output_path, index=False)
    logging.debug(f"sim_id={sim_id} | wrote {len(lineage)} lineage rows")
    return lineage
