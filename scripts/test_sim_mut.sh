#!/usr/bin/env bash
# Smoke test for the full simulation pipeline.
#
# Runs run_sim.py (with --save_lineage) to produce a tree sequence plus the
# pedigree-tracing outputs, then the transmission map + pedigree relatedness
# tools (which only need Stage 1), then mutations, coalescence/mutation time
# tables, strain/host VCFs, and the coalescence-time relatedness extras --
# printing the resulting files so you can eyeball that each stage worked.
#
# Usage: bash scripts/test_sim_mut.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$REPO_ROOT/scripts/test_output"
SEED=123

# delete and recreate the output directory
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

echo "=== Stage 1: simulate ==="
python "$REPO_ROOT/run_sim.py" \
    --Nh 3 --T 10 --K 2 --chi 0.18 --C 1 \
    --recombination_rate 7.4e-7 \
    --random_seed "$SEED" \
    --out_dir "$OUT_DIR" \
    --save_ts --save_ID --save_pairs --save_lineage

TS_FILE="$OUT_DIR/ts_${SEED}.trees"

echo
echo "=== Transmission map (only needs Stage 1's --save_lineage) ==="
python "$REPO_ROOT/run_transmission_map.py" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

echo
echo "=== Stage 2: add mutations ==="
python "$REPO_ROOT/run_mut.py" \
    --ts "$TS_FILE" \
    --mutation_rate 1e-6 \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

MUT_TS_FILE="$OUT_DIR/mutated_${SEED}.trees"

echo
echo "=== Stage 3: coalescence & mutation time tables ==="
python "$REPO_ROOT/run_coal_mut_t.py" \
    --ts "$MUT_TS_FILE" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

echo
echo "=== Pedigree relatedness (genome-wide, needs Stage 3's coal_T) ==="
python "$REPO_ROOT/run_relate_strain_genome.py" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

echo
echo "=== Stage 4: strain VCF ==="
python "$REPO_ROOT/run_vcf_strain.py" \
    --ts "$MUT_TS_FILE" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

STRAIN_VCF="$OUT_DIR/output_${SEED}_strain.vcf"

echo
echo "=== Stage 5: host VCF ==="
python "$REPO_ROOT/run_vcf_host.py" \
    --strain_vcf "$STRAIN_VCF" \
    --random_seed "$SEED" \
    --depth 10 \
    --output_dir "$OUT_DIR" \
    --save_proportions

echo
echo "=== Stage 6: pairwise IBD fractions ==="
python "$REPO_ROOT/extras/relate_t_strain.py" \
    --random_seed "$SEED" \
    --time_boundary 2 \
    --output_dir "$OUT_DIR"

echo
echo "=== Stage 7: host relatedness ==="
python "$REPO_ROOT/extras/relate_t_host.py" \
    --random_seed "$SEED" \
    --output_dir "$OUT_DIR"

echo
echo "=== Output files ==="
ls -la "$OUT_DIR"

echo
echo "=== sim_ind_node_map_${SEED}.csv (first 5 rows) ==="
head -5 "$OUT_DIR/sim_ind_node_map_${SEED}.csv"

echo
echo "=== node_pairs_${SEED}.csv ==="
cat "$OUT_DIR/node_pairs_${SEED}.csv"

echo
echo "=== lineage_${SEED}.csv (first 5 rows) ==="
head -5 "$OUT_DIR/lineage_${SEED}.csv"

echo
echo "=== node_id_map_${SEED}.csv (first 5 rows) ==="
head -5 "$OUT_DIR/node_id_map_${SEED}.csv"

echo
echo "=== transmission_map_${SEED}.html ==="
ls -lh "$OUT_DIR/transmission_map_${SEED}.html"

echo
echo "=== pedigree_relate_strain_genome_${SEED}.csv (first 10 rows) ==="
head -10 "$OUT_DIR/pedigree_relate_strain_genome_${SEED}.csv"

echo
echo "=== tree sequence summary (pre-mutation) ==="
python -c "
import tskit
ts = tskit.load('$TS_FILE')
print(f'samples={ts.num_samples}, trees={ts.num_trees}, sequence_length={ts.sequence_length:.0f}')
"

echo
echo "=== mutated tree sequence summary ==="
python -c "
import tskit
ts = tskit.load('$MUT_TS_FILE')
print(f'samples={ts.num_samples}, sites={ts.num_sites}, mutations={ts.num_mutations}')
"

echo
echo "=== mut_T_${SEED}.csv (first 5 rows) ==="
head -5 "$OUT_DIR/mut_T_${SEED}.csv"

echo
echo "=== coal_T_${SEED}.csv (first 5 rows) ==="
head -5 "$OUT_DIR/coal_T_${SEED}.csv"

echo
echo "=== output_${SEED}_strain.vcf (strain-level, first 15 lines) ==="
head -15 "$STRAIN_VCF"

echo
echo "=== output_${SEED}_host.vcf (host-level, first 15 lines) ==="
head -15 "$OUT_DIR/output_${SEED}_host.vcf"

echo
echo "=== host_proportions_${SEED}.csv ==="
cat "$OUT_DIR/host_proportions_${SEED}.csv"

echo
echo "=== relate_strain_${SEED}.csv ==="
cat "$OUT_DIR/relate_strain_${SEED}.csv"

echo
echo "=== relate_host_${SEED}.csv ==="
cat "$OUT_DIR/relate_host_${SEED}.csv"

echo
echo "Pipeline completed successfully. Output in: $OUT_DIR"

echo
echo "=== Benchmark (optional): hmmIBD vs ground truth ==="
HMMIBD_BIN="${HMMIBD:-hmmibd}"
if command -v "$HMMIBD_BIN" &>/dev/null || [ -x "$HMMIBD_BIN" ]; then
    python "$REPO_ROOT/extras/benchmark_hmmibd.py" \
        --strain_vcf    "$STRAIN_VCF" \
        --relate_strain "$OUT_DIR/relate_strain_${SEED}.csv" \
        --random_seed   "$SEED" \
        --output_dir    "$OUT_DIR" \
        --hmmibd        "$HMMIBD_BIN"

    echo
    echo "=== ibd_benchmark_${SEED}.csv (first 10 rows) ==="
    head -10 "$OUT_DIR/ibd_benchmark_${SEED}.csv"
    echo "(scatter plot saved to $OUT_DIR/ibd_benchmark_${SEED}.png)"
else
    echo "SKIP: hmmIBD not found ('$HMMIBD_BIN'). To run the benchmark:"
    echo "  git clone https://github.com/glipsnort/hmmIBD.git"
    echo "  cd hmmIBD && gcc -o hmmibd hmmIBD.c -lm"
    echo "  Then re-run with: HMMIBD=/path/to/hmmibd bash scripts/test_sim_mut.sh"
fi

echo
echo "=== Regression check ==="
compare_exit=0
python "$REPO_ROOT/scripts/compare_outputs.py" \
    --actual  "$OUT_DIR" \
    --expected "$REPO_ROOT/scripts/expected" || compare_exit=$?

if [ -t 1 ]; then
    echo
    read -r -p "Update expected results with these outputs? [y/N] " reply
    if [[ "$reply" =~ ^[Yy]$ ]]; then
        python "$REPO_ROOT/scripts/compare_outputs.py" \
            --actual   "$OUT_DIR" \
            --expected "$REPO_ROOT/scripts/expected" \
            --bless
    fi
fi

exit $compare_exit
