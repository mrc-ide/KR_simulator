"""Plot mosdepth window depth along each of the 14 P. falciparum chromosomes.

Usage:
    python plot_depth.py <sample>.regions.bed.gz <output>.png
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

regions_bed, out_png = sys.argv[1], sys.argv[2]
sample = regions_bed.split("/")[-1].replace(".regions.bed.gz", "")

df = pd.read_csv(regions_bed, sep="\t", header=None, names=["chrom", "start", "end", "depth"])
nuclear = df[~df["chrom"].str.contains("API|MIT")]
chroms = nuclear["chrom"].unique()
median = nuclear["depth"].median()

# colours: one series, so one hue; text and reference line in neutral inks
LINE, MEDIAN, TEXT, MUTED, GRID = "#2a78d6", "#52514e", "#0b0b0b", "#8a8984", "#e6e5e1"

fig, axes = plt.subplots(len(chroms), 1, figsize=(11, 1.15 * len(chroms)), sharex=True, sharey=True)
max_mb = nuclear.groupby("chrom")["end"].max().max() / 1e6

for ax, chrom in zip(axes, chroms):
    c = nuclear[nuclear["chrom"] == chrom]
    mid_mb = (c["start"] + c["end"]) / 2 / 1e6
    ax.plot(mid_mb, c["depth"].clip(lower=0.5), color=LINE, lw=0.6)  # 0 depth drawn at 0.5 (log axis)
    ax.axhline(median, color=MEDIAN, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_yticks([1, 10, 100, 1000, 10000], ["1", "10", "100", "1,000", "10,000"])
    ax.minorticks_off()
    ax.set_xlim(0, max_mb)
    ax.text(0.003, 0.88, f"{chrom}   mean {c['depth'].mean():.0f}×", transform=ax.transAxes,
            fontsize=8, color=TEXT, va="top")
    ax.grid(axis="y", color=GRID, lw=0.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7, labelcolor=TEXT)

axes[0].set_ylim(0.5, nuclear["depth"].max() * 1.5)
axes[-1].set_xlabel("Position on chromosome (Mb)", fontsize=9, color=TEXT)
fig.supylabel("Mean depth per 500 bp window (log scale)", fontsize=9, color=TEXT)
fig.suptitle(f"{sample}: read depth along each chromosome  (dashed line = genome median, {median:.0f}×)",
             fontsize=11, color=TEXT, x=0.01, ha="left")
fig.tight_layout(rect=(0.01, 0, 1, 0.985))
fig.savefig(out_png, dpi=150, facecolor="white")
print(f"wrote {out_png}")
