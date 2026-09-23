"""Human-mosquito-human transmission for the KR malaria model."""

import numpy as np

from sim.recombination import recomb_poisson

__all__ = ["human_mosq_human"]


def human_mosq_human(tables, parasites_in_human, parent, rng, recomb_event_ID, t,
                     population_id, parasites_in_human_next, child, child_ind_id, params):
    """
    Simulate a single human-mosquito-human transmission event.

    For each oocyst, two parent parasites are sampled and recombined to
    produce 4 recombinant products. K products are then sampled and
    passed to the child host.

    Parameters
    ----------
    tables : tskit.TableCollection
        Tree sequence tables to record nodes and edges into.
    parasites_in_human : dict
        Mapping of host index to list of parasite node IDs.
    parent : int
        Index of the parent host.
    rng : np.random.Generator
        Random number generator.
    recomb_event_ID : int
        Current recombination event counter.
    t : int
        Current generation (time).
    population_id : int
        tskit population ID.
    parasites_in_human_next : dict
        Mapping of child host index to list of parasite node IDs for next generation.
    child : int
        Index of the child host.
    child_ind_id : int
        tskit individual ID for the child.
    params : SimParams
        Simulation parameters.

    Returns
    -------
    parasites_in_human_next : dict
        Updated mapping with child parasite node IDs appended.
    recomb_event_ID : int
        Updated recombination event counter.
    """

    # store all recombinant products and IDs across all oocysts
    recomb_prod_ALL = []
    recomb_ID_ls_ALL = []

    for oocyst in range(params.C):

        # sample 2 parasites for this oocyst (with replacement)
        parent_parasites = rng.choice(parasites_in_human[parent], 2, replace=True)

        # produce 4 recombinant products from the 2 parent parasites
        recomb_prod, recomb_ID_ls = recomb_poisson(
            parent_parasites=parent_parasites,
            params=params,
            recomb_event_ID=recomb_event_ID,
            rng=rng
        )

        recomb_prod_ALL.extend(recomb_prod)

        # offset recombination IDs by oocyst index
        recomb_ID_ls = [x + (oocyst * 4) for x in recomb_ID_ls]
        recomb_ID_ls_ALL.extend(recomb_ID_ls)

    # sample K products from all oocyst recombinant products
    chosen_products = rng.choice(params.C * 4, size=params.K, replace=True)

    for i, prod_idx in enumerate(chosen_products):

        recomb_ID = recomb_ID_ls_ALL[prod_idx]

        # record child parasite node in the tree sequence
        child_node_id = tables.nodes.add_row(
            time=t,
            individual=child_ind_id,
            population=population_id,
            metadata={"recomb_ID": recomb_ID}
        )

        parasites_in_human_next[child].append(child_node_id)

        # record edges from parent segments to child node
        for left, right, parent_node in recomb_prod_ALL[prod_idx]:
            tables.edges.add_row(left, right, parent=parent_node, child=child_node_id)

    # increment recombination event counter for next oocyst group
    recomb_event_ID += (params.C * 4)

    return parasites_in_human_next, recomb_event_ID