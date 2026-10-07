<h1 align="center">krsim — Malaria Transmission Simulator</h1>

<p align="center">
  <a href="https://mrc-ide.github.io/KR_simulator">
    <img alt="Go to the krsim documentation" width="720" src="https://img.shields.io/badge/%F0%9F%93%96%20Go%20to%20the%20documentation%20%E2%86%92-93358F?style=for-the-badge">
  </a>
</p>

<h3 align="center">👉 Installation, tutorials and the full reference all live on the <a href="https://mrc-ide.github.io/KR_simulator">krsim website</a>.</h3>

<p align="center">
  <a href="https://mrc-ide.github.io/KR_simulator/background.html"><img alt="Background" height="34" src="https://img.shields.io/badge/Background-93358F?style=for-the-badge"></a>
  <a href="https://mrc-ide.github.io/KR_simulator/installation.html"><img alt="Installation" height="34" src="https://img.shields.io/badge/Installation-93358F?style=for-the-badge"></a>
  <a href="https://mrc-ide.github.io/KR_simulator/tutorial-small.html"><img alt="Tutorial" height="34" src="https://img.shields.io/badge/Tutorial-93358F?style=for-the-badge"></a>
  <a href="https://mrc-ide.github.io/KR_simulator/commands.html"><img alt="Command line" height="34" src="https://img.shields.io/badge/Command%20line-93358F?style=for-the-badge"></a>
  <a href="https://mrc-ide.github.io/KR_simulator/scaling.html"><img alt="Running at scale" height="34" src="https://img.shields.io/badge/Running%20at%20scale-93358F?style=for-the-badge"></a>
</p>

---

A forward simulator that brings together malaria modelling and population genetics to model parasite transmission. Whole-genome information is tracked through every recombination event as an ancestral recombination graph, stored efficiently with [tskit](https://tskit.dev). From the recorded genealogy (a tree sequence), krsim overlays mutations and outputs coalescence times across the genome and VCF files at the strain level and at the (observed) host level.

**To install and use krsim, follow the [installation guide on the website](https://mrc-ide.github.io/KR_simulator/installation.html).**

## Development

To work on the code itself, see [DEVELOPMENT.md](DEVELOPMENT.md).
