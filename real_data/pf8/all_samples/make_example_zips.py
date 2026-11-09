"""
The two example zips for `krsim vcf-host-depth --depth_zip`: 100 sWGA and 100 gDNA Pf8 samples.

Samples are picked at random (SEED below) from the QC-pass samples that have window stats. Each is checked
with krsim's check_profile (3D7 chromosome names and lengths, 500 bp windows); a sample that fails is
skipped and the next random one taken. Writes, in $PF8/example_zips/:
  pf8_sWGA_100.zip, pf8_gDNA_100.zip   <sample>.window_stats.tsv files, ready for --depth_zip
  samples.tsv                          the samples in each zip: sample, type, country, year, study

The zips are meant for a GitHub release of krsim, not for git (about 47 MB each).

Usage: PYTHONPATH=<repo>/src python make_example_zips.py
"""

import os
import zipfile

import numpy as np
import pandas as pd

from krsim.output.pf8_profiles import check_profile

PF8 = "/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth"
OUT = f"{PF8}/example_zips"
SEED = 2026
N = 100
TYPES = ["sWGA", "gDNA"]


def main():
    os.makedirs(OUT, exist_ok=True)
    info = pd.read_csv(f"{PF8}/output/Pf8_sample_info.tsv", sep="\t")
    info = info[info["QC pass"] & info["Sample type"].isin(TYPES)]
    info = info[[os.path.exists(f"{PF8}/output/{s}.window_stats.tsv") for s in info["Sample"]]]
    rng = np.random.default_rng(SEED)
    rows = []
    for t in TYPES:
        pool = info[info["Sample type"] == t].sort_values("Sample")
        order = rng.permutation(len(pool))
        zip_path = f"{OUT}/pf8_{t}_{N}.zip"
        picked = 0
        with zipfile.ZipFile(zip_path + ".part", "w", compression=zipfile.ZIP_DEFLATED) as z:
            for i in order:
                r = pool.iloc[i]
                path = f"{PF8}/output/{r['Sample']}.window_stats.tsv"
                try:
                    check_profile(pd.read_csv(path, sep="\t"), path)
                except ValueError as exc:
                    print(f"  skipped {r['Sample']}: {exc}")
                    continue
                z.write(path, arcname=os.path.basename(path))
                rows.append({"zip": os.path.basename(zip_path), "sample": r["Sample"], "type": t,
                             "country": r["Country"], "year": int(r["Year"]), "study": r["Study"]})
                picked += 1
                if picked == N:
                    break
        os.replace(zip_path + ".part", zip_path)
        print(f"{zip_path}: {picked} samples, {os.path.getsize(zip_path) / 1e6:.1f} MB (from {len(pool):,} QC-pass)")
    pd.DataFrame(rows).to_csv(f"{OUT}/samples.tsv", sep="\t", index=False)
    print(f"wrote {OUT}/samples.tsv")


if __name__ == "__main__":
    main()
