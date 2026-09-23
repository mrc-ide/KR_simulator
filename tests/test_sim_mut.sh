#!/usr/bin/env bash
# Smoke test for the full simulation pipeline.
#
# Runs `krsim simulate` (with --save_lineage) to produce a tree sequence plus
# the pedigree-tracing outputs, then the strain transmission map + pairwise
# strain pedigree map (which only need Stage 1), then mutations,
# coalescence/mutation time tables, strain/host VCFs, and the coalescence-time
# relatedness extras -- printing the resulting files so you can eyeball that
# each stage worked.
#
# Requires krsim to be installed (`pip install -e .` from the repo root).
#
# Usage: bash tests/test_sim_mut.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# run from the repo root so paths print short (tests/test_output/...)
cd "$REPO_ROOT"
OUT_DIR="tests/test_output"
SEED=123

# simulation parameters for this smoke test (the only place these live --
# nothing else in this script duplicates these values in a comment)
NH=3
T=10
K=2
CHI=0.18
C=1
RECOMB_RATE=7.4e-7

# colors auto-disable when not a real terminal (piped to a log file, CI, etc.)
if [ -t 1 ]; then
    BOLD=$'\033[1m'; GREEN=$'\033[32m'; CYAN=$'\033[36m'; YELLOW=$'\033[33m'; RESET=$'\033[0m'
else
    BOLD=""; GREEN=""; CYAN=""; YELLOW=""; RESET=""
fi

# one-line content check printed under a stage, e.g. for the binary .trees
# files the regression check can't compare
check() {
    python -c "$1" | sed "s/^/  ${GREEN}✓${RESET} /"
    echo
}

# runs a labeled stage under a colored banner (krsim itself reports what it
# wrote and how long it took)
run_stage() {
    local label="$1"; shift
    echo "${BOLD}${CYAN}=== ${label} ===${RESET}"
    "$@"
}

# delete and recreate the output directory
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

echo
echo "Outputs: $OUT_DIR"
echo "Log:     $OUT_DIR/KR_sim_${SEED}.log  (every command adds its details here)"
echo

run_stage "Stage 1: Simulate genealogy" \
    krsim simulate \
    --Nh "$NH" --T "$T" --K "$K" --chi "$CHI" --C "$C" \
    --recombination_rate "$RECOMB_RATE" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR" \
    --save_pairs --save_lineage

TS_FILE="$OUT_DIR/ts_${SEED}.trees"
check "
import tskit
ts = tskit.load('$TS_FILE')
print(f'ts_${SEED}.trees: {ts.num_trees} trees, {ts.num_samples} samples, {ts.num_nodes} nodes')"

run_stage "Strain transmission map (only needs Stage 1's --save_lineage)" \
    krsim strain-transmission-map \
    --lineage_csv "$OUT_DIR/lineage_${SEED}.csv" \
    --node_id_map_csv "$OUT_DIR/node_id_map_${SEED}.csv" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

run_stage "Pairwise strain pedigree map (only needs Stage 1's --save_lineage)" \
    krsim pairwise-strain-pedigree-map \
    --lineage_csv "$OUT_DIR/lineage_${SEED}.csv" \
    --node_id_map_csv "$OUT_DIR/node_id_map_${SEED}.csv" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

run_stage "Stage 2: Add mutations" \
    krsim mutate \
    --ts "$TS_FILE" \
    --mutation_rate 1e-6 \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

MUT_TS_FILE="$OUT_DIR/mutated_${SEED}.trees"
check "
import tskit
ts = tskit.load('$MUT_TS_FILE')
print(f'mutated_${SEED}.trees: {ts.num_mutations} mutations at {ts.num_sites} sites across {ts.num_trees} trees')"

run_stage "Stage 3: Coalescence & mutation time tables" \
    krsim coal \
    --ts "$MUT_TS_FILE" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"
check "
import pandas as pd
coal = pd.read_csv('$OUT_DIR/coal_T_${SEED}.csv')
mut = pd.read_csv('$OUT_DIR/mut_T_${SEED}.csv')
pairs = coal[['sample_1', 'sample_2']].drop_duplicates()
print(f'coal_T_${SEED}.csv: {len(coal):,} coalescence events (one per strain pair per tree; {len(pairs)} strain pairs coalesce)')
print(f'mut_T_${SEED}.csv: {len(mut):,} mutations')"

run_stage "Stage 4: Strain VCF" \
    krsim vcf-strain \
    --ts "$MUT_TS_FILE" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

STRAIN_VCF="$OUT_DIR/output_${SEED}_strain.vcf"

run_stage "Stage 5: Host VCF" \
    krsim vcf-host \
    --strain_vcf "$STRAIN_VCF" \
    --ind_node_map_csv "$OUT_DIR/sim_ind_node_map_${SEED}.csv" \
    --random_seed "$SEED" \
    --depth 10 \
    --output_dir "$OUT_DIR" \
    --save_proportions

run_stage "Extras: pairwise IBD fractions (krsim relate-strain)" \
    krsim relate-strain \
    --coal_csv "$OUT_DIR/coal_T_${SEED}.csv" \
    --random_seed "$SEED" \
    --time_boundary 2 \
    --output_dir "$OUT_DIR"

run_stage "Extras: host relatedness (krsim relate-host)" \
    krsim relate-host \
    --relate_csv "$OUT_DIR/relate_strain_${SEED}.csv" \
    --ind_node_map_csv "$OUT_DIR/sim_ind_node_map_${SEED}.csv" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

echo
echo "${GREEN}${BOLD}Pipeline completed successfully.${RESET} Output in: $OUT_DIR"

echo
echo "${BOLD}=== Regression check ===${RESET}"
compare_exit=0
python tests/compare_outputs.py \
    --actual  "$OUT_DIR" \
    --expected "tests/expected" || compare_exit=$?

if [ -t 1 ]; then
    echo
    read -r -p "${BOLD}Update expected results with these outputs? [y/N]${RESET} " reply
    if [[ "$reply" =~ ^[Yy]$ ]]; then
        python tests/compare_outputs.py \
            --actual   "$OUT_DIR" \
            --expected "tests/expected" \
            --bless
    fi
fi

exit $compare_exit
