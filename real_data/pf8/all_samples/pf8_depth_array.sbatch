#!/usr/bin/env bash
# SLURM array: pf8_depth_task.sh for every sample in a list, one array task per sample.
#
# Usage, from this folder (set --array to 1-<number of samples>, %K = at most K at once):
#   sbatch --array=1-16%16 pf8_depth_array.sbatch /data/home/bty208/scratch/IMPERIAL/pf8_mosdepth/input/timing_subset.tsv
#SBATCH --job-name=pf8_depth
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth/logs/%x_%A_%a.out

module load miniforge/26.3.2
eval "$(conda shell.bash hook)"
conda activate /data/home/bty208/.conda/envs/mosdepth

bash pf8_depth_task.sh $1
