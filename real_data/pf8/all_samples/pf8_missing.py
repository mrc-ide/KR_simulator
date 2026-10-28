"""
Every Pf8 sample in Pf8_crams.tsv that has no window stats yet, as a list to run again.

Unlike $PF8/failed.tsv (written by a task that fails), this also catches tasks SLURM killed
(time limit, memory, node failure) or that never ran. The reason a sample failed, when known,
is the last line for it in failed.tsv.

Writes $PF8/input/Pf8_crams_missing.tsv (same columns as Pf8_crams.tsv; in parts of 10000 if longer)
and prints the sbatch command for each, e.g.
  sbatch --array=1-<number missing>%50 pf8_depth_array.sbatch $PF8/input/Pf8_crams_missing.tsv

Run only when no pf8_depth jobs are queued or running, or samples still in the queue are listed too.

Usage: python pf8_missing.py
"""

import os

import pandas as pd

# Folder for the data (change this on another computer)
PF8 = "/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth"

crams = pd.read_csv(f"{PF8}/input/Pf8_crams.tsv", sep="\t")
done = crams["sample"].map(lambda s: os.path.getsize(p) > 0 if os.path.exists(p := f"{PF8}/output/{s}.window_stats.tsv") else False)
missing = crams[~done]
missing.to_csv(f"{PF8}/input/Pf8_crams_missing.tsv", sep="\t", index=False)

print(f"{done.sum()} of {len(crams)} samples done, {len(missing)} missing "
      f"({missing.qc_pass.sum()} of them QC pass)")
if os.path.exists(f"{PF8}/failed.tsv"):
    failed = pd.read_csv(f"{PF8}/failed.tsv", sep="\t", header=None,
                         names=["sample", "type", "qc_pass", "cram_bytes", "url", "stage", "task", "date"])
    reasons = failed.drop_duplicates("sample", keep="last").set_index("sample").stage
    print("stage each missing sample last failed at:")
    print(missing["sample"].map(reasons).fillna("not in failed.tsv (killed by SLURM or not run yet)")
          .value_counts().to_string())
# one SLURM array holds at most 10000 tasks, so a long list goes in parts
if len(missing):
    print("\nrun them again (from all_samples/; at most 20000 tasks can be queued at once):")
for i, start in enumerate(range(0, len(missing), 10000), 1):
    part = missing.iloc[start:start + 10000]
    path = f"{PF8}/input/Pf8_crams_missing.tsv" if len(missing) <= 10000 else f"{PF8}/input/Pf8_crams_missing_part{i}.tsv"
    part.to_csv(path, sep="\t", index=False)
    print(f"  sbatch --array=1-{len(part)}%50 pf8_depth_array.sbatch {path}")
