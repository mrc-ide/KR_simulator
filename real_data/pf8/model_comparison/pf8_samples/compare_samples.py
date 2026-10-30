"""Compare depth spread across samples: one mean-vs-SD panel per sample, plus a summary table.

Each panel is like plot_mean_sd.py: every dot is a 500 bp window, and the dashed
line is the Poisson SD (sqrt(mean)). All panels share the same axes.

Usage:
    python plot_compare_samples.py <output>.png <summary>.tsv <a>.window_stats.tsv <b>.window_stats.tsv ...
"""
import math
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

out_png, out_tsv, stats_files = sys.argv[1], sys.argv[2], sys.argv[3:]

# sample details (country, year) from the Pf8 metadata, which sits next to the inputs
meta = pd.read_csv("/Users/bty208/Downloads/Pf8_samples.txt", sep="\t").set_index("Sample")

rows, data = [], []
for path in stats_files:
    sample = path.split("/")[-1].replace(".window_stats.tsv", "")
    df = pd.read_csv(path, sep="\t")
    nuc = df[~df["chrom"].str.contains("API|MIT")]
    pos = nuc[nuc["mean"] > 0]
    rows.append({
        "sample": sample,
        "type": meta.loc[sample, "Sample type"],
        "country": meta.loc[sample, "Country"],
        "year": int(meta.loc[sample, "Year"]),
        "median_depth": nuc["mean"].median(),
        "windows_under_5x": (nuc["mean"] < 5).mean(),
        "windows_above_poisson": (nuc["sd"] ** 2 > nuc["mean"]).mean(),
        "median_variance_over_mean": (pos["sd"] ** 2 / pos["mean"]).median(),
    })
    data.append(pos[pos["sd"] > 0])  # log axes cannot show 0

summary = pd.DataFrame(rows)
summary.round(3).to_csv(out_tsv, sep="\t", index=False)
print(summary.round(2).to_string(index=False))

# colours: one series per panel, so one hue; reference line and text in neutral inks
POINT, REF, TEXT, MUTED, GRID = "#2a78d6", "#52514e", "#0b0b0b", "#8a8984", "#e6e5e1"

n = len(data)
ncol = 4
nrow = math.ceil(n / ncol)
fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.9 * nrow), sharex=True, sharey=True, squeeze=False)
x_all = pd.concat(d["mean"] for d in data)
y_all = pd.concat(d["sd"] for d in data)
x = np.logspace(np.log10(x_all.min()), np.log10(x_all.max()), 200)
ticks = [0.1, 1, 10, 100, 1000, 10000]

for ax, d, r in zip(axes.flat, data, rows):
    ax.scatter(d["mean"], d["sd"], s=2, color=POINT, alpha=0.12, lw=0, rasterized=True)
    ax.plot(x, np.sqrt(x), color=REF, lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ticks, [f"{t:,g}" for t in ticks])
    ax.set_yticks(ticks, [f"{t:,g}" for t in ticks])
    ax.minorticks_off()
    ax.set_title(f"{r['sample']}  {meta.loc[r['sample'], 'Sample type']}  {r['country']} {r['year']}\n"
                 f"median {r['median_depth']:.0f}×,  {r['windows_above_poisson']:.0%} above Poisson",
                 fontsize=8, color=TEXT, loc="left")
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7, labelcolor=TEXT)

for ax in list(axes.flat)[n:]:
    ax.set_visible(False)
axes[0, 0].set_xlim(x_all.min() * 0.8, x_all.max() * 1.3)
axes[0, 0].set_ylim(y_all.min() * 0.8, y_all.max() * 1.3)
fig.supxlabel("Mean depth in window", fontsize=10, color=TEXT)
fig.supylabel("SD of depth in window", fontsize=10, color=TEXT)
counts = summary["type"].value_counts()
fig.suptitle(f"Pf8 samples ({counts.get('sWGA', 0)} sWGA, {counts.get('gDNA', 0)} gDNA): spread of depth within each 500 bp window (dashed line = Poisson, SD = √mean)",
             fontsize=11, color=TEXT, x=0.01, ha="left")
fig.tight_layout(rect=(0.01, 0.01, 1, 0.97))
fig.savefig(out_png, dpi=150, facecolor="white")
print(f"wrote {out_png} and {out_tsv}")
