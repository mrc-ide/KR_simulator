"""Core simulation engine for the KR malaria transmission model."""

import numpy as np
import tskit
import collections

from krsim.sim.params import SimParams, GENOME_LENGTH
from krsim.sim.transmission import human_mosq_human

__all__ = ["KRSimulator"]


class KRSimulator:
    """
    Forward simulator for P. falciparum transmission and recombination.

    Parameters
    ----------
    params : SimParams
        Simulation parameters.
    """

    def __init__(self, params):
        """Initialise the simulator with the given parameters."""
        self.params = params
        self.rng = params.make_rng()

        # initialise tskit table collection
        self.tables = tskit.TableCollection(GENOME_LENGTH)
        self.tables.time_units = "generations"

        # set metadata schemas
        schema = tskit.MetadataSchema.permissive_json()
        self.tables.populations.metadata_schema = schema
        self.tables.individuals.metadata_schema = schema
        self.tables.nodes.metadata_schema = schema

        # add single population
        self.population_id = self.tables.populations.add_row(metadata={"Region": "London"})

        # recombination event counter
        self.recomb_event_ID = 0

    def _initialise_population(self):
        """
        Initialise the founder population at time T.

        Returns
        -------
        parasites_in_human : collections.defaultdict
            Mapping of host index to list of parasite node IDs.
        """
        params = self.params

        # determine number of parasites per host based on superinfection probability
        parasites_per_human = np.where(
            self.rng.random(params.Nh) > params.chi,
            params.K,
            2 * params.K
        )

        parasites_in_human = collections.defaultdict(list)

        for i, num_parasites in enumerate(parasites_per_human):
            human_id = self.tables.individuals.add_row(metadata={"type": "founder"})
            for _ in range(num_parasites):
                parasites_in_human[i].append(
                    self.tables.nodes.add_row(
                        time=params.T,
                        individual=human_id,
                        population=self.population_id
                    )
                )

        return parasites_in_human

    def _run_generation(self, t, parasites_in_human):
        """
        Simulate a single generation of transmission.

        Parameters
        ----------
        t : int
            Current generation (time).
        parasites_in_human : collections.defaultdict
            Parasite node IDs for the current generation.

        Returns
        -------
        parasites_in_human_next : collections.defaultdict
            Parasite node IDs for the next generation.
        """
        params = self.params

        parasites_in_human_next = collections.defaultdict(list)

        # determine superinfection status for each child host
        parasites_per_human = np.where(
            self.rng.random(params.Nh) > params.chi,
            params.K,
            2 * params.K
        )

        for child in range(params.Nh):

            child_ind_id = self.tables.individuals.add_row(metadata={"type": "child"})

            if parasites_per_human[child] == params.K:
                # monoinfected - single parent, single transmission
                parent = self.rng.choice(params.Nh)
                parasites_in_human_next, self.recomb_event_ID = human_mosq_human(
                    tables=self.tables,
                    parasites_in_human=parasites_in_human,
                    parent=parent,
                    rng=self.rng,
                    recomb_event_ID=self.recomb_event_ID,
                    t=t,
                    population_id=self.population_id,
                    parasites_in_human_next=parasites_in_human_next,
                    child=child,
                    child_ind_id=child_ind_id,
                    params=params
                )

            else:
                # superinfected - two parents, two transmissions
                parent1, parent2 = self.rng.choice(params.Nh, size=2, replace=False)

                for parent in [parent1, parent2]:
                    parasites_in_human_next, self.recomb_event_ID = human_mosq_human(
                        tables=self.tables,
                        parasites_in_human=parasites_in_human,
                        parent=parent,
                        rng=self.rng,
                        recomb_event_ID=self.recomb_event_ID,
                        t=t,
                        population_id=self.population_id,
                        parasites_in_human_next=parasites_in_human_next,
                        child=child,
                        child_ind_id=child_ind_id,
                        params=params
                    )

        return parasites_in_human_next

    def run(self, capture_lineage=False):
        """
        Run the full forward simulation.

        Parameters
        ----------
        capture_lineage : bool
            If True, also return a snapshot of the tables taken before
            simplify() runs. simplify() drops any node that isn't a branch
            point relevant to the final generation, which can make an edge
            in the final tree sequence skip over several generations of
            real transmission events. The raw snapshot keeps every edge
            from every generation, at the cost of its own node numbering
            not matching the final (simplified) tree sequence.

        Returns
        -------
        ts : tskit.TreeSequence
            The simulated tree sequence.
        raw_tables : tskit.TableCollection or None
            Snapshot of the tables before simplify(), or None if
            capture_lineage is False.
        node_map : numpy.ndarray or None
            Mapping from raw_tables node ID to its ID in the final,
            simplified tree sequence (tskit.NULL for nodes that did not
            survive simplification), or None if capture_lineage is False.
        """
        params = self.params

        # initialise founder population at time T
        parasites_in_human = self._initialise_population()

        # forward simulation from T to 0
        t = params.T
        while t > 0:
            t -= 1
            parasites_in_human = self._run_generation(t, parasites_in_human)

        raw_tables = self.tables.copy() if capture_lineage else None

        # sort and simplify with respect to current generation
        self.tables.sort()
        node_map = self.tables.simplify(
            np.where(self.tables.nodes.time == 0)[0],
            keep_input_roots=True
        )

        return self.tables.tree_sequence(), raw_tables, (node_map if capture_lineage else None)