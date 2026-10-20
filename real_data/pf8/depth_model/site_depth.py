"""Step 1: the depth at each SNP, from a real sample's depth windows.

For every site in a krsim strain VCF, find the window it falls in and record
that window's mean and SD of depth.

Usage:
    python site_depth.py <strain>.vcf <sample>.window_stats.tsv <output>.tsv
(window_stats.tsv is made by window_stats.py)
"""
import sys

import pandas as pd

strain_vcf, window_stats, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

# windows: chromosome, start (0-based), end, mean and SD of depth
windows = pd.read_csv(window_stats, sep="\t")
window_size = int((windows["end"] - windows["start"]).mode()[0])  # 500 bp
stats = {(row.chrom, row.start): (row.mean, row.sd) for row in windows.itertuples()}

# SNP positions from the strain VCF: krsim names chromosomes 1-14, Pf8 names them Pf3D7_01_v3 ...
sites = pd.read_csv(strain_vcf, sep="\t", comment="#", header=None, usecols=[0, 1], names=["chrom", "pos"])
sites["pf8_chrom"] = sites["chrom"].map(lambda c: f"Pf3D7_{int(c):02d}_v3")

# VCF positions start at 1, mosdepth windows at 0: position p is in the window starting at (p - 1) // 500 * 500
sites["window_start"] = (sites["pos"] - 1) // window_size * window_size
sites["window_end"] = sites["window_start"] + window_size
window = [stats[(c, s)] for c, s in zip(sites["pf8_chrom"], sites["window_start"])]
sites["mean"] = [m for m, _ in window]
sites["sd"] = [sd for _, sd in window]

sites.to_csv(out_tsv, sep="\t", index=False)
print(f"{len(sites)} sites -> {out_tsv}")
print(sites[["mean", "sd"]].describe().round(1).to_string())
