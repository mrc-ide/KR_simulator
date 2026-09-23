#!/usr/bin/env bash
# Smoke test for the full simulation pipeline.
#
# Runs `krsim simulate` (with --save_lineage) to produce a tree sequence plus
# the pedigree-tracing outputs, then the transmission map + pairwise
# coalescence tools (which only need Stage 1), then mutations,
# coalescence/mutation time tables, strain/host VCFs, and the coalescence-time
# relatedness extras -- printing the resulting files so you can eyeball that
# each stage worked.
#
# Requires krsim to be installed (`pip install -e .` from the repo root).
#
# Usage: bash tests/test_sim_mut.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$REPO_ROOT/tests/test_output"
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

# runs a labeled stage, timing it and printing a colored banner + elapsed time
run_stage() {
    local label="$1"; shift
    echo
    echo "${BOLD}${CYAN}=== ${label} ===${RESET}"
    local start end elapsed
    start=$(date +%s.%N)
    "$@"
    end=$(date +%s.%N)
    elapsed=$(awk -v s="$start" -v e="$end" 'BEGIN{printf "%.1f", e-s}')
    echo "${GREEN}✓ done in ${elapsed}s${RESET}"
}

# delete and recreate the output directory
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

run_stage "Stage 1: Simulate genealogy" \
    krsim simulate \
    --Nh "$NH" --T "$T" --K "$K" --chi "$CHI" --C "$C" \
    --recombination_rate "$RECOMB_RATE" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR" \
    --save_ts --save_ID --save_pairs --save_lineage

TS_FILE="$OUT_DIR/ts_${SEED}.trees"

run_stage "Transmission map (only needs Stage 1's --save_lineage)" \
    krsim transmission-map \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

run_stage "Pairwise coalescence (only needs Stage 1's --save_lineage)" \
    krsim pairwise-coalescence \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

run_stage "Stage 2: Add mutations" \
    krsim mutate \
    --ts "$TS_FILE" \
    --mutation_rate 1e-6 \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

MUT_TS_FILE="$OUT_DIR/mutated_${SEED}.trees"

run_stage "Stage 3: Coalescence & mutation time tables" \
    krsim coal \
    --ts "$MUT_TS_FILE" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

run_stage "Stage 4: Strain VCF" \
    krsim vcf-strain \
    --ts "$MUT_TS_FILE" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

STRAIN_VCF="$OUT_DIR/output_${SEED}_strain.vcf"

run_stage "Stage 5: Host VCF" \
    krsim vcf-host \
    --strain_vcf "$STRAIN_VCF" \
    --random_seed "$SEED" \
    --depth 10 \
    --output_dir "$OUT_DIR" \
    --save_proportions

run_stage "Extras: pairwise IBD fractions (krsim relate-strain)" \
    krsim relate-strain \
    --random_seed "$SEED" \
    --time_boundary 2 \
    --output_dir "$OUT_DIR"

run_stage "Extras: host relatedness (krsim relate-host)" \
    krsim relate-host \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

# prints a bounded preview of a CSV: row count + first N data rows
preview_csv() {
    local label="$1" path="$2" n="${3:-5}"
    echo
    if [ ! -f "$path" ]; then
        echo "${BOLD}=== ${label} (missing) ===${RESET}"
        return
    fi
    echo "${BOLD}=== ${label} ($(($(wc -l < "$path") - 1)) rows) ===${RESET}"
    head -n $((n + 1)) "$path"
}

# prints a bounded preview of a VCF, skipping the noisy ## meta-header lines
# (writes to a temp file rather than piping into head, since `grep | head`
# makes grep exit via SIGPIPE once head has enough lines, which aborts the
# whole script under `set -o pipefail`)
preview_vcf() {
    local label="$1" path="$2" n="${3:-5}"
    echo
    echo "${BOLD}=== ${label} ===${RESET}"
    local tmp
    tmp=$(mktemp)
    grep -v "^##" "$path" > "$tmp"
    head -n $((n + 1)) "$tmp"
    rm -f "$tmp"
}

echo
echo "${BOLD}=== Output files ===${RESET}"
ls -lh "$OUT_DIR"

preview_csv "sim_ind_node_map_${SEED}.csv" "$OUT_DIR/sim_ind_node_map_${SEED}.csv"
preview_csv "node_pairs_${SEED}.csv"       "$OUT_DIR/node_pairs_${SEED}.csv"
preview_csv "lineage_${SEED}.csv"          "$OUT_DIR/lineage_${SEED}.csv"
preview_csv "node_id_map_${SEED}.csv"      "$OUT_DIR/node_id_map_${SEED}.csv"

echo
echo "${BOLD}=== transmission_map_${SEED}.html ===${RESET}"
ls -lh "$OUT_DIR/transmission_map_${SEED}.html"

echo
echo "${BOLD}=== pairwise_coalescence_${SEED}.html ===${RESET}"
ls -lh "$OUT_DIR/pairwise_coalescence_${SEED}.html"

echo
echo "${BOLD}=== tree sequence summary (pre-mutation) ===${RESET}"
python -c "
import tskit
ts = tskit.load('$TS_FILE')
print(f'samples={ts.num_samples}, trees={ts.num_trees}, sequence_length={ts.sequence_length:.0f}')
"

echo
echo "${BOLD}=== mutated tree sequence summary ===${RESET}"
python -c "
import tskit
ts = tskit.load('$MUT_TS_FILE')
print(f'samples={ts.num_samples}, sites={ts.num_sites}, mutations={ts.num_mutations}')
"

preview_csv "mut_T_${SEED}.csv"            "$OUT_DIR/mut_T_${SEED}.csv"
preview_csv "coal_T_${SEED}.csv"           "$OUT_DIR/coal_T_${SEED}.csv"
preview_vcf "output_${SEED}_strain.vcf"    "$STRAIN_VCF"
preview_vcf "output_${SEED}_host.vcf"      "$OUT_DIR/output_${SEED}_host.vcf"
preview_csv "host_proportions_${SEED}.csv" "$OUT_DIR/host_proportions_${SEED}.csv"
preview_csv "relate_strain_${SEED}.csv"    "$OUT_DIR/relate_strain_${SEED}.csv"
preview_csv "relate_host_${SEED}.csv"      "$OUT_DIR/relate_host_${SEED}.csv"

echo
echo "${GREEN}${BOLD}Pipeline completed successfully.${RESET} Output in: $OUT_DIR"

echo
echo "${BOLD}=== Regression check ===${RESET}"
compare_exit=0
python "$REPO_ROOT/tests/compare_outputs.py" \
    --actual  "$OUT_DIR" \
    --expected "$REPO_ROOT/tests/expected" || compare_exit=$?

if [ -t 1 ]; then
    echo
    read -r -p "${BOLD}Update expected results with these outputs? [y/N]${RESET} " reply
    if [[ "$reply" =~ ^[Yy]$ ]]; then
        python "$REPO_ROOT/tests/compare_outputs.py" \
            --actual   "$OUT_DIR" \
            --expected "$REPO_ROOT/tests/expected" \
            --bless
    fi
fi

exit $compare_exit
