"""
Default P. falciparum rates used when the user doesn't pass their own.

Kept free of numpy/tskit imports so the CLI can read them without slowing
down `krsim --help`.
"""

# per bp per generation recombination rate
RECOMBINATION_RATE = 7.4e-7

# per bp per generation mutation rate
MUTATION_RATE = 1e-6
