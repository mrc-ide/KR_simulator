"""Simulation parameter definitions for the KR malaria model."""

from dataclasses import dataclass, field
import numpy as np

__all__ = ["COMP_BREAKPOINTS", "GENOME_LENGTH", "SimParams"]

# P. falciparum chromosome boundary positions in base pairs: the cumulative ends
# of chromosomes 1-13 in the 3D7 reference used by MalariaGEN Pf8
# (PlasmoDB-54, Pf3D7_01_v3 ... Pf3D7_14_v3; https://zenodo.org/records/18681980)
COMP_BREAKPOINTS = np.array([
    640_851, 1_587_953, 2_655_924, 3_856_414, 5_199_971, 6_618_213, 8_063_420,
    9_536_225, 11_077_960, 12_765_616, 14_803_956, 17_075_450, 20_000_686,
], dtype=np.int64)

# P. falciparum 3D7 nuclear genome length in base pairs (14 chromosomes, PlasmoDB-54)
GENOME_LENGTH = 23_292_622

@dataclass
class SimParams:
    """
    Container for stochastic simulation parameters.

    Attributes
    ----------
    Nh : int
        Number of human hosts in the population (must be >= 1).
    T : int
        Number of generations to simulate (must be >= 1).
    K : int
        Number of parasites transmitted from mosquito to human host per bite (must be >= 1).
    chi : float
        Probability of superinfection (must be 0.0 <= chi <= 1.0).
    C : int
        Number of oocysts per mosquito (must be >= 1).
    recombination_rate : float
        Per base pair recombination rate (must be >= 0).
    random_seed : int
        Seed for reproducibility (must be >= 0).
    """

    Nh: int
    T: int
    K: int
    chi: float
    C: int
    recombination_rate: float
    random_seed: int
    genome_length: int = GENOME_LENGTH
    comp_breakpoints: np.ndarray = field(default_factory=lambda: COMP_BREAKPOINTS.copy())

    def make_rng(self):
        """Create a reproducible random number generator for the simulation."""
        return np.random.default_rng(self.random_seed)

    def __post_init__(self):
        """Validate all parameters after initialisation."""
        if not isinstance(self.Nh, int) or self.Nh <= 0:
            raise ValueError(f"Nh (number of human hosts) must be an integer >= 1, got {self.Nh}")

        if not isinstance(self.T, int) or self.T <= 0:
            raise ValueError(f"T (number of generations) must be an integer >= 1, got {self.T}")

        if not isinstance(self.K, int) or self.K <= 0:
            raise ValueError(f"K (number of parasites transmitted per mosquito bite) must be an integer >= 1, got {self.K}")

        if not isinstance(self.chi, float) or not 0.0 <= self.chi <= 1.0:
            raise ValueError(f"chi (superinfection probability) must be a float between 0.0 and 1.0, got {self.chi}")

        if not isinstance(self.C, int) or self.C <= 0:
            raise ValueError(f"C (number of oocysts per mosquito) must be an integer >= 1, got {self.C}")

        if not isinstance(self.recombination_rate, float) or self.recombination_rate < 0:
            raise ValueError(f"recombination_rate (per bp recombination rate) must be a float >= 0.0, got {self.recombination_rate}")

        if not isinstance(self.random_seed, int) or self.random_seed < 0:
            raise ValueError(f"random_seed must be an integer >= 0, got {self.random_seed}")
        
        if not isinstance(self.genome_length, int) or self.genome_length <= 0:
            raise ValueError(f"genome_length must be an integer >= 1, got {self.genome_length}")

        if not isinstance(self.comp_breakpoints, np.ndarray) or len(self.comp_breakpoints) == 0:
            raise ValueError(f"comp_breakpoints must be a non-empty numpy array")

        if not np.all(self.comp_breakpoints[:-1] < self.comp_breakpoints[1:]):
            raise ValueError("comp_breakpoints must be strictly increasing")

        if self.comp_breakpoints[-1] > self.genome_length:
            raise ValueError(f"comp_breakpoints must not exceed genome_length={self.genome_length}")