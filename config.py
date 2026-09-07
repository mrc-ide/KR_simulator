"""Load simulation parameters from a YAML config file."""

import yaml
import numpy as np
from sim.params import SimParams, GENOME_LENGTH, COMP_BREAKPOINTS

def load_config(path):
    """
    Load simulation parameters from a YAML config file.
    
    Parameters
    ----------
    path : str
        Path to the YAML config file.
    
    Returns
    -------
    SimParams
        Simulation parameters.
    """
    with open(path) as f:
        cfg = yaml.safe_load(f)

    return SimParams(
        Nh=cfg["Nh"],
        T=cfg["T"],
        K=cfg["K"],
        chi=cfg["chi"],
        C=cfg["C"],
        recombination_rate=cfg["recombination_rate"],
        random_seed=cfg["random_seed"],
        genome_length=cfg.get("genome_length", GENOME_LENGTH),
        comp_breakpoints=np.array(cfg["comp_breakpoints"], dtype=np.int64)
                        if "comp_breakpoints" in cfg
                        else COMP_BREAKPOINTS.copy()
    )