#!/usr/bin/env bash
# SLURM array: fit_chosen_models.py on every QC-pass Pf8 gDNA and sWGA sample, split over the array tasks.
# ~10 s per sample per core; 24,145 samples over 100 tasks x 8 cores is ~5 min per task.
#
# Usage, from this folder (the last argument must match the --array size):
#   sbatch --array=0-99%50 fit_chosen_models.sbatch 100
#   then, when every task is done:  sbatch --cpus-per-task=16 --wrap "../../../.venv/bin/python fit_chosen_models.py merge"
#SBATCH --job-name=pf8_fit
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth/logs/%x_%A_%a.out

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
../../../.venv/bin/python fit_chosen_models.py fit $SLURM_ARRAY_TASK_ID $1
