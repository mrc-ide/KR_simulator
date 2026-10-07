<h1 align="center">krsim — Malaria Transmission Simulator</h1>

<p align="center">
  <a href="https://mrc-ide.github.io/KR_simulator">
    <img alt="Go to the krsim documentation" width="720" src=".github/readme/documentation.svg">
  </a>
</p>

<h3 align="center">👉 Installation, tutorials and the full reference all live on the <a href="https://mrc-ide.github.io/KR_simulator">krsim website</a>.</h3>

<p align="center">
  <a href="https://mrc-ide.github.io/KR_simulator/background.html"><img alt="Background" src=".github/readme/background.svg"></a>
  <a href="https://mrc-ide.github.io/KR_simulator/installation.html"><img alt="Installation" src=".github/readme/installation.svg"></a>
  <a href="https://mrc-ide.github.io/KR_simulator/tutorial-small.html"><img alt="Tutorial" src=".github/readme/tutorial.svg"></a>
  <a href="https://mrc-ide.github.io/KR_simulator/commands.html"><img alt="Command line" src=".github/readme/command-line.svg"></a>
  <a href="https://mrc-ide.github.io/KR_simulator/scaling.html"><img alt="Running at scale" src=".github/readme/scaling.svg"></a>
</p>

---

A forward simulator that brings together malaria modelling and population genetics to model parasite transmission. Whole-genome information is tracked through every recombination event as an ancestral recombination graph, stored efficiently with [tskit](https://tskit.dev). From the recorded genealogy (a tree sequence), krsim overlays mutations and outputs coalescence times across the genome and VCF files at the strain level and at the (observed) host level.

**To install and use krsim, follow the [installation guide on the website](https://mrc-ide.github.io/KR_simulator/installation.html).**

## Development

To work on the code itself, see [DEVELOPMENT.md](DEVELOPMENT.md).
