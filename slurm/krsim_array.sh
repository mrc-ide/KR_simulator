#!/bin/bash
# Run many krsim simulations as a SLURM array job: one simulation per array task.
#
# Each task reads one line of params.tsv (after the header), runs the krsim
# pipeline for it (simulate -> mutate -> vcf-strain -> vcf-host) and writes the
# results to results/Nh<Nh>_T<T>_K<K>_chi<chi>_C<C>/.
#
# Submit from this folder:
#   mkdir -p logs
#   sbatch --array=1-$(($(wc -l < params.tsv) - 1)) krsim_array.sh
#
# Each command uses one core. Request memory for the largest command and time
# for all of them together. Every task gets the --mem and --time below, so set
# them for the largest line in params.tsv, or put large runs in their own params
# file and override on submission. See "Running at scale" on the krsim website
# for what each run size needs, e.g. for Nh=1000, T=1000:
#   sbatch --mem=32G --time=02:00:00 --array=1-5 krsim_array.sh params_large.tsv

#SBATCH --job-name=krsim
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/krsim_%A_%a.out

set -euo pipefail

PARAMS=${1:-params.tsv}

# Which commands to run. The VCFs are the largest outputs (about 20 GB for the
# strain VCF and 10 GB for the host VCF at Nh=1000), so leave out any you don't need.
STEPS="simulate mutate vcf-strain vcf-host"
DEPTH=10   # sequencing depth for vcf-host

# activate the krsim environment (adjust for your cluster, e.g. `module load` first)
eval "$(conda shell.bash hook)"
conda activate krsim

# line 1 is the header, so array task N reads line N + 1
read -r NH T K CHI C SEED < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PARAMS")
OUT="results/Nh${NH}_T${T}_K${K}_chi${CHI}_C${C}"
mkdir -p "$OUT"

echo "task ${SLURM_ARRAY_TASK_ID}: Nh=${NH} T=${T} K=${K} chi=${CHI} C=${C} seed=${SEED} -> ${OUT}"

for step in $STEPS; do
    case $step in
        simulate)
            krsim simulate --Nh "$NH" --T "$T" --K "$K" --chi "$CHI" --C "$C" \
                --random_seed "$SEED" --output_dir "$OUT" ;;
        mutate)
            krsim mutate --ts "$OUT/ts_${SEED}.trees" --random_seed "$SEED" --output_dir "$OUT" ;;
        vcf-strain)
            krsim vcf-strain --ts "$OUT/ts_mut_${SEED}.trees" --random_seed "$SEED" --output_dir "$OUT" ;;
        vcf-host)
            krsim vcf-host --strain_vcf "$OUT/output_${SEED}_strain.vcf" \
                --ind_node_map_csv "$OUT/sim_ind_node_map_${SEED}.csv" --depth "$DEPTH" \
                --random_seed "$SEED" --output_dir "$OUT" ;;
        *)
            echo "unknown step: $step" >&2; exit 1 ;;
    esac
done
