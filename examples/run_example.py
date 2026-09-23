"""
Run the full KR Malaria Transmission Simulator pipeline for the worked example.

Five stages, each producing files consumed by the next:
  1. krsim simulate   — simulate genealogy → tree sequence + ID map + node pairs + lineage
  2. krsim mutate     — add mutations to tree sequence → mutated tree sequence
  3. krsim coal       — extract coalescence & mutation times → coal_T, mut_T CSVs
  4. krsim vcf-strain — write per-strain VCF from mutated tree sequence
  5. krsim vcf-host   — aggregate to host level with Dirichlet-Multinomial sampling

Plus two pedigree-based visualization tools, both run right after simulation
since they only need Stage 1's --save_lineage output:
  - krsim transmission-map — renders an interactive HTML pedigree diagram.
  - krsim pairwise-coalescence — renders an interactive HTML page for
    tracing any two strains back to their shared pedigree ancestor(s).

Requires krsim to be installed (`pip install -e .` from the repo root).

Usage (from repo root):
    python examples/run_example.py
"""

import os
import sys
import subprocess

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")
SEED       = 123

os.makedirs(OUTPUT_DIR, exist_ok=True)

# `python -m krsim` is the same as the `krsim` command, but guarantees the
# krsim installed in *this* Python environment is the one that runs
KRSIM = [sys.executable, "-m", "krsim"]

def run(label, cmd):
    print(f"\n=== {label} ===")
    print("  krsim " + " ".join(cmd[3:]))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)


TS_FILE      = os.path.join(OUTPUT_DIR, f"ts_{SEED}.trees")
MUT_TS_FILE  = os.path.join(OUTPUT_DIR, f"mutated_{SEED}.trees")
STRAIN_VCF   = os.path.join(OUTPUT_DIR, f"output_{SEED}_strain.vcf")

run("Stage 1: simulate", [
    *KRSIM, "simulate",
    "--Nh",               "3",
    "--T",                "10",
    "--K",                "2",
    "--chi",              "0.18",
    "--C",                "1",
    "--recombination_rate", "7.4e-7",
    "--random_seed",      str(SEED),
    "--output_dir",       OUTPUT_DIR,
    "--save_pairs",
    "--save_lineage",
])

run("Transmission map: visualize simulated pedigree", [
    *KRSIM, "transmission-map",
    "--random_seed", str(SEED),
    "--output_dir",  OUTPUT_DIR,
])

run("Pairwise coalescence: trace two strains to a shared ancestor", [
    *KRSIM, "pairwise-coalescence",
    "--random_seed", str(SEED),
    "--output_dir",  OUTPUT_DIR,
])

run("Stage 2: add mutations", [
    *KRSIM, "mutate",
    "--ts",           TS_FILE,
    "--mutation_rate", "1e-6",
    "--random_seed",  str(SEED),
    "--output_dir",   OUTPUT_DIR,
])

run("Stage 3: coalescence & mutation times", [
    *KRSIM, "coal",
    "--ts",          MUT_TS_FILE,
    "--random_seed", str(SEED),
    "--output_dir",  OUTPUT_DIR,
])

run("Stage 4: strain VCF", [
    *KRSIM, "vcf-strain",
    "--ts",          MUT_TS_FILE,
    "--random_seed", str(SEED),
    "--output_dir",  OUTPUT_DIR,
])

run("Stage 5: host VCF + strain proportions", [
    *KRSIM, "vcf-host",
    "--strain_vcf",  STRAIN_VCF,
    "--random_seed", str(SEED),
    "--depth",       "10",
    "--output_dir",  OUTPUT_DIR,
    "--save_proportions",
])

print(f"\nDone. Outputs written to {OUTPUT_DIR}/")
