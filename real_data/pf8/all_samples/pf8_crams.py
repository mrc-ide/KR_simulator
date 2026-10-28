"""
Every Pf8 sample's CRAM (URL and size, from ENA) in one table, and a small test subset of it.

Needs in $PF8/input: Pf8_samples.txt (Zenodo) and Pf8_ena_runs.tsv, the ENA list of the project's runs:
  curl -G "https://www.ebi.ac.uk/ena/portal/api/search" --data-urlencode 'result=read_run' \
    --data-urlencode 'query=study_accession=PRJEB96783' --data-urlencode 'limit=0' \
    --data-urlencode 'fields=run_accession,sample_alias,submitted_bytes,submitted_ftp' \
    --data-urlencode 'format=tsv' -o Pf8_ena_runs.tsv

Writes in $PF8/input:
  Pf8_crams.tsv      every sample with a CRAM: sample, type, QC pass, CRAM bytes, URL (QC-pass samples first)
  Pf8_crams_part<i>.tsv  the same in parts of 10000, one SLURM array each
  timing_subset.tsv  N QC-pass sWGA and N QC-pass gDNA samples, spread evenly over each type's CRAM sizes,
                     plus the 99th percentile and the largest CRAM of each type (the evenly spread ones stop at ~3.5 GB)

  random_subset.tsv  R random sWGA and R random gDNA samples (QC pass or not), none from timing_subset.tsv

and $PF8/output/Pf8_sample_info.tsv: every Pf8 sample's metadata (QC pass, exclusion reason, % callable, ...)
with its CRAM size, for filtering the window stats later.

Usage: python pf8_crams.py [N per sample type, default 6] [R per sample type, default 50]
"""

import os
import sys

import numpy as np
import pandas as pd

# Folder for the data (change this on another computer)
PF8 = "/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
PART_SIZE = 10000

# one row per CRAM file (a run lists its CRAM and .crai, separated by ";")
runs = pd.read_csv(f"{PF8}/input/Pf8_ena_runs.tsv", sep="\t").dropna(subset=["submitted_ftp"])
crams = pd.DataFrame(
    [(r.sample_alias, int(size), f"https://{path}")
     for r in runs.itertuples()
     for size, path in zip(str(r.submitted_bytes).split(";"), r.submitted_ftp.split(";"))
     if path.endswith(".cram")],
    columns=["sample", "cram_bytes", "url"])

meta = pd.read_csv(f"{PF8}/input/Pf8_samples.txt", sep="\t")
table = (meta[["Sample", "Sample type", "QC pass"]]
         .rename(columns={"Sample": "sample", "Sample type": "type", "QC pass": "qc_pass"})
         .merge(crams, on="sample"))
print(f"{len(table)} of {len(meta)} Pf8 samples have a CRAM in ENA ({table.cram_bytes.sum() / 1e12:.1f} TB)")

# the list for the full run, QC-pass samples first, split so each part fits in one SLURM array (at most 10000 tasks)
full = table.sort_values("qc_pass", ascending=False, kind="stable")
full.to_csv(f"{PF8}/input/Pf8_crams.tsv", sep="\t", index=False)
for i, start in enumerate(range(0, len(full), PART_SIZE), 1):
    full.iloc[start:start + PART_SIZE].to_csv(f"{PF8}/input/Pf8_crams_part{i}.tsv", sep="\t", index=False)

# every Pf8 sample's metadata (QC pass, exclusion reason, % callable, ...) and CRAM size, next to the window stats
os.makedirs(f"{PF8}/output", exist_ok=True)
meta.merge(crams[["sample", "cram_bytes"]], left_on="Sample", right_on="sample", how="left").drop(columns="sample") \
    .to_csv(f"{PF8}/output/Pf8_sample_info.tsv", sep="\t", index=False)
print(table.groupby(["type", "qc_pass"]).cram_bytes.agg(["count", "median", "sum"])
      .assign(median_GB=lambda d: (d["median"] / 1e9).round(2), total_TB=lambda d: (d["sum"] / 1e12).round(2))
      [["count", "median_GB", "total_TB"]])

# test subset: for each type, the samples nearest N evenly spaced size quantiles (2.5% ... 97.5%)
picks = []
for sample_type in ["sWGA", "gDNA"]:
    pool = table[(table.type == sample_type) & table.qc_pass].sort_values("cram_bytes").reset_index(drop=True)
    rows = np.round((np.arange(N) + 0.5) / N * (len(pool) - 1)).astype(int)
    rows = np.append(rows, [round(0.99 * (len(pool) - 1)), len(pool) - 1])   # and the large ones
    picks.append(pool.iloc[rows])
subset = pd.concat(picks)
subset.to_csv(f"{PF8}/input/timing_subset.tsv", sep="\t", index=False)
print(f"\ntiming_subset.tsv ({len(subset)} samples):")
print(subset.assign(GB=(subset.cram_bytes / 1e9).round(2))[["sample", "type", "GB"]].to_string(index=False))

# second test: R random samples of each type (QC pass or not, so sizes follow Pf8's real mix), not in the first
R = int(sys.argv[2]) if len(sys.argv) > 2 else 50
pool = table[~table["sample"].isin(subset["sample"])]
random_subset = pd.concat([pool[pool.type == t].sample(R, random_state=1) for t in ["sWGA", "gDNA"]])
random_subset.to_csv(f"{PF8}/input/random_subset.tsv", sep="\t", index=False)
print(f"\nrandom_subset.tsv: {R} sWGA + {R} gDNA, {random_subset.cram_bytes.sum() / 1e9:.0f} GB, "
      f"{(~random_subset.qc_pass).sum()} failed QC")
