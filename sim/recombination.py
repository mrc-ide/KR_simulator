"""Recombination functions for the KR malaria transmission model."""

import numpy as np

from sim.params import SimParams, GENOME_LENGTH

__all__ = ["recomb_event", "recomb_poisson"]


def recomb_event(params, parent_parasites, rng):
    """
    Simulate a single recombination event between two parent parasites.

    Combines compulsory chromosome boundary breaks (Bernoulli 0.5) with
    Poisson-distributed random breakpoints across the genome.

    Parameters
    ----------
    params : SimParams
        Simulation parameters.
    parent_parasites : tuple
        Two parent parasite node IDs (p0, p1).
    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    tuple
        Two recombinant products (product1, product2), each a list of
        (start, end, parent) segments.
    """

    # Bernoulli(0.5) for each chromosome boundary
    keep = rng.binomial(1, 0.5, size=len(params.comp_breakpoints)).astype(bool)
    compulsory_breaks = params.comp_breakpoints[keep]

    # parents
    p0, p1 = parent_parasites

    # Poisson distributed random breakpoints across the genome
    num_recombinations = rng.poisson(params.recombination_rate * GENOME_LENGTH)
    breakpoints = rng.uniform(0, GENOME_LENGTH, size=num_recombinations)

    # combine compulsory + random breakpoints + start and end then sort
    breakpoints = np.sort(np.concatenate((breakpoints, compulsory_breaks, [0], [GENOME_LENGTH])))

    # build segments from consecutive breakpoint pairs
    segments = list(zip(breakpoints[:-1], breakpoints[1:]))

    # alternate segments between parents to build two recombinant products
    product1 = []
    product2 = []

    for i, (start, end) in enumerate(segments):
        if i % 2 == 0:
            product1.append((start, end, p0))
            product2.append((start, end, p1))
        else:
            product1.append((start, end, p1))
            product2.append((start, end, p0))

    return product1, product2


def recomb_poisson(parent_parasites, params, recomb_event_ID, rng):
    """
    Simulate two independent recombination events producing 4 recombinant products.

    Parameters
    ----------
    parent_parasites : tuple
        Two parent parasite node IDs.
    params : SimParams
        Simulation parameters.
    recomb_event_ID : int
        Current recombination event counter.
    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    tuple
        products : list of 4 recombinant products.
        ids : list of 4 recombination IDs.
    """

    # assign 4 sequential recombination IDs for this event group
    ids = list(range(recomb_event_ID, recomb_event_ID + 4))

    products = []

    # two independent shuffling events each producing 2 recombinant products
    for _ in range(2):
        prod1, prod2 = recomb_event(params, parent_parasites, rng)
        products.extend([prod1, prod2])

    return products, ids