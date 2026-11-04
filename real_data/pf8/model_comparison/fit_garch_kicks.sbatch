#!/usr/bin/env bash
# SLURM array: fit_garch_kicks.py (GARCH t kicks) on every sample fitted by fit_chosen_models.py, split over the array tasks.
# ~6 s per sample per core; 24,145 samples over 100 tasks x 8 cores is ~3 min per task.
#
# Usage, from this folder (the last argument must match the --array size):
#   sbatch --array=0-99%50 fit_garch_kicks.sbatch 100
#   then, when every task is done:  ../../../.venv/bin/python fit_garch_kicks.py merge   (quick: reads only the tables)
#SBATCH --job-name=pf8_garch
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth/logs/%x_%A_%a.out

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
../../../.venv/bin/python fit_garch_kicks.py fit $SLURM_ARRAY_TASK_ID $1
