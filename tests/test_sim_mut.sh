#!/usr/bin/env bash
# Smoke test for the full simulation pipeline.
#
# Runs `krsim simulate` (with --save_lineage) to produce a tree sequence plus
# the pedigree-tracing outputs, then the strain transmission map + pairwise
# strain pedigree map (which only need Stage 1), then mutations,
# coalescence/mutation time tables, strain/host VCFs (including the host VCF with
# real read depth, from the three small depth profiles in tests/data/test_depth.zip, and
# with modelled sWGA and gDNA read depth),
# the depth-zip helper (on a tiny made-up mosdepth file), and the coalescence-time
# relatedness extras -- printing the resulting files so you can eyeball that
# each stage worked.
#
# Requires krsim to be installed (`pip install -e .` from the repo root).
#
# Usage: bash tests/test_sim_mut.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# run from the repo root so paths print short (tests/output/...)
cd "$REPO_ROOT"
OUT_DIR="tests/output"
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

MUT_TS_FILE="$OUT_DIR/ts_mut_${SEED}.trees"
check "
import tskit
ts = tskit.load('$MUT_TS_FILE')
print(f'ts_mut_${SEED}.trees: {ts.num_mutations} mutations at {ts.num_sites} sites across {ts.num_trees} trees')"

run_stage "Strain transmission map (needs Stage 1's --save_lineage and Stage 2's mutations)" \
    krsim strain-transmission-map \
    --lineage_csv "$OUT_DIR/lineage_${SEED}.csv" \
    --node_id_map_csv "$OUT_DIR/node_id_map_${SEED}.csv" \
    --ts "$MUT_TS_FILE" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

run_stage "Stage 3: Coalescence & mutation time tables" \
    krsim coal-mut \
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

# depth profiles of three real Pf8 sWGA samples, kept in tests/data/test_depth.zip so nothing is downloaded
run_stage "Stage 5: Host VCF, read depth from a zip of depth files" \
    krsim vcf-host-depth \
    --strain_vcf "$STRAIN_VCF" \
    --ind_node_map_csv "$OUT_DIR/sim_ind_node_map_${SEED}.csv" \
    --depth_zip "tests/data/test_depth.zip" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

# modelled depth, one run per sample type: both write output_<seed>_host_depth.vcf etc., so each runs in its
# own folder and its files are moved into $OUT_DIR renamed *_model_<type>.* for the regression check
for TYPE in sWGA gDNA; do
    MODEL_TMP="$OUT_DIR/model_$TYPE"
    run_stage "Stage 5: Host VCF, modelled $TYPE read depth" \
        krsim vcf-host-depth \
        --strain_vcf "$STRAIN_VCF" \
        --ind_node_map_csv "$OUT_DIR/sim_ind_node_map_${SEED}.csv" \
        --depth_model "$TYPE" \
        --random_seed "$SEED" \
        --output_dir "$MODEL_TMP"
    mv "$MODEL_TMP/output_${SEED}_host_depth.vcf" "$OUT_DIR/output_${SEED}_host_depth_model_${TYPE}.vcf"
    mv "$MODEL_TMP/host_depth_profiles_${SEED}.csv" "$OUT_DIR/host_depth_profiles_${SEED}_model_${TYPE}.csv"
    mv "$MODEL_TMP/host_depth_proportions_${SEED}.csv" "$OUT_DIR/host_depth_proportions_${SEED}_model_${TYPE}.csv"
    rm -rf "$MODEL_TMP"
done
check "
import numpy as np
from krsim.output.depth_model import model_settings, relative_track, host_depths
from krsim.output.pf8_profiles import CHROMS
for t in ['sWGA', 'gDNA']:
    s = model_settings(t)
    tr = relative_track(np.random.default_rng(1), s)
    allw = np.concatenate([tr[c] for c in CHROMS])
    assert allw.max() <= s['cap'] + 1e-9, 'cap not applied'
    k = int(s['end_width'] // 500)
    if k:
        ends = np.concatenate([np.r_[tr[c][:k], tr[c][-k:]] for c in CHROMS])
        assert np.median(ends) < 0.5 * np.median(allw), 'end zone not lowered'
    assert np.all(host_depths(np.random.default_rng(1), model_settings(t, depth_spread=0.0), 5) == s['depth_mean'])
    flat = model_settings(t, depth_spread=0.0)
    flat.update(scale=0.0, end_width=0.0)
    ft = np.concatenate([relative_track(np.random.default_rng(1), flat)[c] for c in CHROMS])
    assert np.allclose(ft, 1.0), 'scale 0 should give a flat track'
    print(f'{t} model: median window {np.median(allw):.2f}x (expect about 1), highest {allw.max():.0f}x '
          f'(cap {s[\"cap\"]:.0f}x), end zone {s[\"end_width\"]:,.0f} bp; --depth_spread 0 gives one typical depth; scale 0 gives flat depth')"

# a tiny made-up mosdepth per-base file: depth 30 along every 3D7 chromosome
# (written outside $OUT_DIR, so the regression check doesn't compare it)
DEPTH_TMP="$(mktemp -d)"
python -c "
import gzip
from krsim.output.pf8_profiles import CHROM_LENGTHS
with gzip.open('$DEPTH_TMP/toy.per-base.bed.gz', 'wt') as f:
    for chrom, length in CHROM_LENGTHS.items():
        f.write(f'{chrom}\t0\t{length}\t30\n')"

run_stage "Helper: zip of depth files from mosdepth output (krsim depth-zip)" \
    krsim depth-zip \
    --bed "$DEPTH_TMP/toy.per-base.bed.gz" \
    --output "$DEPTH_TMP/toy_depth.zip"
check "
import zipfile
names = zipfile.ZipFile('$DEPTH_TMP/toy_depth.zip').namelist()
assert names == ['toy.window_stats.tsv'], names
print(f'toy_depth.zip: {len(names)} depth file ({names[0]})')"
rm -rf "$DEPTH_TMP"

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
    --output   "$OUT_DIR" \
    --baseline "tests/baseline" || compare_exit=$?

if [ -t 1 ]; then
    echo
    read -r -p "${BOLD}Update the baseline with these outputs? [y/N]${RESET} " reply
    if [[ "$reply" =~ ^[Yy]$ ]]; then
        python tests/compare_outputs.py \
            --output   "$OUT_DIR" \
            --baseline "tests/baseline" \
            --bless
    fi
fi

exit $compare_exit
