"""Plot each window's depth SD against its mean, with the Poisson line (SD = sqrt(mean)).

Windows above the line vary more than a Poisson draw would (over-dispersed).

Usage:
    python plot_mean_sd.py <sample>.window_stats.tsv <output>.png
(window_stats.tsv is made by window_stats.py)
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

window_stats, out_png = sys.argv[1], sys.argv[2]
sample = window_stats.split("/")[-1].replace(".window_stats.tsv", "")

df = pd.read_csv(window_stats, sep="\t")
nuclear = df[~df["chrom"].str.contains("API|MIT")]
shown = nuclear[(nuclear["mean"] > 0) & (nuclear["sd"] > 0)]  # log axes cannot show 0
above = (shown["sd"] ** 2 > shown["mean"]).mean()

# colours: one series, so one hue; reference line and text in neutral inks
POINT, REF, TEXT, MUTED, GRID = "#2a78d6", "#52514e", "#0b0b0b", "#8a8984", "#e6e5e1"

fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(shown["mean"], shown["sd"], s=3, color=POINT, alpha=0.15, lw=0, rasterized=True)

x = np.logspace(np.log10(shown["mean"].min()), np.log10(shown["mean"].max()), 200)
ax.plot(x, np.sqrt(x), color=REF, lw=1.2, ls="--")
label_x = 3000  # a sparse part of the plot, so the label sits clear of the points
ax.text(label_x, np.sqrt(label_x) * 0.75, "Poisson: SD = √mean", color=TEXT, fontsize=9, ha="center", va="top")

ax.set_xscale("log")
ax.set_yscale("log")
ticks = [0.1, 1, 10, 100, 1000, 10000]
ax.set_xticks(ticks, [f"{t:,g}" for t in ticks])
ax.set_yticks(ticks, [f"{t:,g}" for t in ticks])
ax.set_xlim(shown["mean"].min() * 0.8, shown["mean"].max() * 1.3)
ax.set_ylim(shown["sd"].min() * 0.8, shown["sd"].max() * 1.3)
ax.minorticks_off()
ax.set_xlabel("Mean depth in window", fontsize=10, color=TEXT)
ax.set_ylabel("SD of depth in window", fontsize=10, color=TEXT)
ax.set_title(f"{sample}: spread of depth within each 500 bp window\n"
             f"{len(shown):,} nuclear windows; {above:.0%} lie above the Poisson line (more spread than Poisson)",
             fontsize=10, color=TEXT, loc="left")
ax.grid(color=GRID, lw=0.5)
ax.set_axisbelow(True)  # grid behind the points
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color(MUTED)
ax.tick_params(colors=MUTED, labelsize=8, labelcolor=TEXT)

fig.tight_layout()
fig.savefig(out_png, dpi=150, facecolor="white")
print(f"wrote {out_png}")
