"""Plot window depth (mean, with a mean +/- SD band) along each of the 14 P. falciparum chromosomes.

Usage:
    python plot_depth.py <sample>.window_stats.tsv <output>.png                       # all 14 chromosomes
    python plot_depth.py <sample>.window_stats.tsv <output>.png Pf3D7_07_v3:350000-450000   # zoom into one region
(window_stats.tsv is made by step 8 of pf8_mosdepth.sh, or by `krsim vcf-host-depth`; the SD band is only readable when zoomed in)
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

window_stats, out_png = sys.argv[1], sys.argv[2]
region = sys.argv[3] if len(sys.argv) > 3 else None
sample = window_stats.split("/")[-1].replace(".window_stats.tsv", "")

df = pd.read_csv(window_stats, sep="\t").rename(columns={"mean": "depth"})
nuclear = df[~df["chrom"].str.contains("API|MIT")]
median = nuclear["depth"].median()  # genome median, also used for the zoomed view

if region:
    chrom, span = region.split(":")
    lo, hi = (int(x) for x in span.replace(",", "").split("-"))
    nuclear = nuclear[(nuclear["chrom"] == chrom) & (nuclear["end"] > lo) & (nuclear["start"] < hi)]
chroms = nuclear["chrom"].unique()

# colours: one series, so one hue (line, and a lighter band for the SD); text and reference line in neutral inks
LINE, MEDIAN, TEXT, MUTED, GRID = "#2a78d6", "#52514e", "#0b0b0b", "#8a8984", "#e6e5e1"

fig, axes = plt.subplots(len(chroms), 1, figsize=(11, 3.2 if region else 1.15 * len(chroms)),
                         sharex=True, sharey=True, squeeze=False)
axes = axes[:, 0]
x_lo, x_hi = (lo / 1e6, hi / 1e6) if region else (0, nuclear.groupby("chrom")["end"].max().max() / 1e6)

for ax, chrom in zip(axes, chroms):
    c = nuclear[nuclear["chrom"] == chrom]
    mid_mb = (c["start"] + c["end"]) / 2 / 1e6
    # depth 0 is drawn at 0.5, the bottom of the log axis
    ax.fill_between(mid_mb, (c["depth"] - c["sd"]).clip(lower=0.5), (c["depth"] + c["sd"]).clip(lower=0.5),
                    color=LINE, alpha=0.25, lw=0)
    ax.plot(mid_mb, c["depth"].clip(lower=0.5), color=LINE, lw=0.6)
    ax.axhline(median, color=MEDIAN, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_yticks([1, 10, 100, 1000, 10000], ["1", "10", "100", "1,000", "10,000"])
    ax.minorticks_off()
    ax.set_xlim(x_lo, x_hi)
    ax.text(0.003, 0.88, f"{chrom}   mean {c['depth'].mean():.0f}×", transform=ax.transAxes,
            fontsize=8, color=TEXT, va="top")
    ax.grid(axis="y", color=GRID, lw=0.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7, labelcolor=TEXT)

axes[0].set_ylim(0.5, (nuclear["depth"] + nuclear["sd"]).max() * 1.5)
axes[-1].set_xlabel("Position on chromosome (Mb)", fontsize=9, color=TEXT)
fig.supylabel("Mean depth per 500 bp window (log scale)", fontsize=9, color=TEXT)
fig.suptitle(f"{sample}: read depth {'in ' + region if region else 'along each chromosome'}  (line = window mean, band = ± SD, "
             f"dashed = genome median {median:.0f}×)",
             fontsize=11, color=TEXT, x=0.01, ha="left")
fig.tight_layout(rect=(0.01, 0, 1, 0.985))
fig.savefig(out_png, dpi=150, facecolor="white")
print(f"wrote {out_png}")
