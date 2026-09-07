"""
Run the full KR Malaria Transmission Simulator pipeline for the worked example.

Five stages, each producing files consumed by the next:
  1. run_sim.py       — simulate genealogy → tree sequence + ID map + node pairs + lineage
  2. run_mut.py       — add mutations to tree sequence → mutated tree sequence
  3. run_coal_mut_t.py — extract coalescence & mutation times → coal_T, mut_T CSVs
  4. run_vcf_strain.py — write per-strain VCF from mutated tree sequence
  5. run_vcf_host.py  — aggregate to host level with Dirichlet-Multinomial sampling

Plus two pedigree-based extras:
  - run_transmission_map.py — renders an interactive HTML pedigree diagram
    from Stage 1's --save_lineage output alone; runs immediately after
    simulation, before mutations/VCFs, so you can look at the raw
    transmission history right away.
  - run_relate_strain_genome.py — genome-wide pairwise relatedness, one row
    per (sample pair, tree interval), read straight off Stage 3's coal_T
    (annotated with mirror-pair info from Stage 1's lineage); run last,
    since it needs coal_T to exist.

Usage (from repo root):
    python examples/run_example.py
"""

import os
import sys
import subprocess

REPO_ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")
SEED       = 123

os.makedirs(OUTPUT_DIR, exist_ok=True)

def run(label, cmd):
    print(f"\n=== {label} ===")
    print("  " + " ".join(cmd[2:]))
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode != 0:
        sys.exit(result.returncode)


TS_FILE      = os.path.join(OUTPUT_DIR, f"ts_{SEED}.trees")
MUT_TS_FILE  = os.path.join(OUTPUT_DIR, f"mutated_{SEED}.trees")
STRAIN_VCF   = os.path.join(OUTPUT_DIR, f"output_{SEED}_strain.vcf")

run("Stage 1: simulate", [
    sys.executable, "run_sim.py",
    "--Nh",               "3",
    "--T",                "10",
    "--K",                "2",
    "--chi",              "0.18",
    "--C",                "1",
    "--recombination_rate", "7.4e-7",
    "--random_seed",      str(SEED),
    "--out_dir",          OUTPUT_DIR,
    "--save_ts",
    "--save_ID",
    "--save_pairs",
    "--save_lineage",
])

run("Transmission map: visualize simulated pedigree", [
    sys.executable, "run_transmission_map.py",
    "--random_seed", str(SEED),
    "--output_dir",  OUTPUT_DIR,
])

run("Stage 2: add mutations", [
    sys.executable, "run_mut.py",
    "--ts",           TS_FILE,
    "--mutation_rate", "1e-6",
    "--random_seed",  str(SEED),
    "--output_dir",   OUTPUT_DIR,
])

run("Stage 3: coalescence & mutation times", [
    sys.executable, "run_coal_mut_t.py",
    "--ts",          MUT_TS_FILE,
    "--random_seed", str(SEED),
    "--output_dir",  OUTPUT_DIR,
])

run("Stage 4: strain VCF", [
    sys.executable, "run_vcf_strain.py",
    "--ts",          MUT_TS_FILE,
    "--random_seed", str(SEED),
    "--output_dir",  OUTPUT_DIR,
])

run("Stage 5: host VCF + strain proportions", [
    sys.executable, "run_vcf_host.py",
    "--strain_vcf",  STRAIN_VCF,
    "--random_seed", str(SEED),
    "--depth",       "10",
    "--output_dir",  OUTPUT_DIR,
    "--save_proportions",
])

run("Pedigree relatedness (genome-wide)", [
    sys.executable, "run_relate_strain_genome.py",
    "--random_seed", str(SEED),
    "--output_dir",  OUTPUT_DIR,
])

print(f"\nDone. Outputs written to {OUTPUT_DIR}/")
