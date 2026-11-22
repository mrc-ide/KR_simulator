#!/usr/bin/env bash
# SLURM array: one krsim simulation per task for coal_sweep.py (51 runs: 17 settings x 3 seeds; COAL_SET picks the parameter set).
# The largest runs (Nh 400) take ~10 min and ~11 GB.
#
# Usage, from this folder:
#   sbatch --array=0-50 coal_sweep.sbatch
#   then, when every task is done:  ../../.venv/bin/python coal_sweep.py plot
#SBATCH --job-name=coal_sweep
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/data/home/bty208/scratch/IMPERIAL/krsim_coal_sweep/logs/%x_%A_%a.out

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
../../.venv/bin/python coal_sweep.py run $SLURM_ARRAY_TASK_ID
