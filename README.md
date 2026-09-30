# krsim — Malaria Transmission Simulator

A forward simulator that brings together malaria modelling and population genetics to model parasite transmission. Whole-genome information is tracked through every recombination event as an ancestral recombination graph, stored efficiently with [tskit](https://tskit.dev). From the recorded genealogy (a tree sequence), krsim overlays mutations and outputs coalescence times across the genome and VCF files at the strain level and at the (observed) host level.

**Documentation: <https://mrc-ide.github.io/KR_simulator>**

## Installation

You need [conda](https://conda-forge.org/download/) (Miniforge, Miniconda or Anaconda). Download the environment file and create the environment:

```bash
curl -O https://raw.githubusercontent.com/mrc-ide/KR_simulator/main/environment-user.yml
conda env create -f environment-user.yml
conda activate krsim
```

Check it worked with `krsim --help`, which lists every command.

## Getting started

- [Installation](https://mrc-ide.github.io/KR_simulator/installation.html): install krsim.
- [Tutorial](https://mrc-ide.github.io/KR_simulator/tutorial.html): a small example with a Jupyter notebook to play with, and the command-line scripts to run the pipeline.
- [Command line](https://mrc-ide.github.io/KR_simulator/commands.html): every command and option.

## Development

To work on the code itself, see [DEVELOPMENT.md](DEVELOPMENT.md).
