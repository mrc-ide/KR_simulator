"""Compare read depth across samples, using their 500 bp window stats.

Writes into <output_dir>:
  depth_boxplot.png        spread of window depth in every sample (box = middle half of windows, line = median)
  depth_<type>.png         one per sample type (sWGA, gDNA): the samples of that type together,
                           mean depth per 500 bp window (line) and mean +/- SD (band), one colour per sample,
                           along chr7 and zoomed into crt

Usage:
    python compare_depth.py <output_dir> <a>.window_stats.tsv <b>.window_stats.tsv ...
"""
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
import pandas as pd

out_dir, stats_files = sys.argv[1], sys.argv[2:]
os.makedirs(out_dir, exist_ok=True)
meta = pd.read_csv("/Users/bty208/Downloads/Pf8_samples.txt", sep="\t").set_index("Sample")
meta["Country"] = meta["Country"].replace({"Democratic Republic of the Congo": "DR Congo"})

samples = {}
for path in stats_files:
    sample = path.split("/")[-1].replace(".window_stats.tsv", "")
    df = pd.read_csv(path, sep="\t")
    samples[sample] = df[~df["chrom"].str.contains("API|MIT")]

# colours: sample type in the box plot; one hue per sample in the per-type plots; text in neutral inks
TYPE_COLOUR = {"sWGA": "#2a78d6", "gDNA": "#eb6834"}
SAMPLE_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a4fd6",
                  "#d6336c", "#0b8a96", "#7a5c2e", "#5b6b7a", "#9bb31a"]
TEXT, MUTED, GRID = "#0b0b0b", "#8a8984", "#e6e5e1"


def style(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,g}"))
    ax.minorticks_off()
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8, labelcolor=TEXT)


# --- spread of window depth in every sample ---
# grouped by sample type, then from lowest to highest median depth
names = sorted(samples, key=lambda s: (list(TYPE_COLOUR).index(meta.loc[s, "Sample type"]), samples[s]["mean"].median()))
fig, ax = plt.subplots(figsize=(max(8, 0.65 * len(names) + 2), 5.5))
types = [meta.loc[s, "Sample type"] for s in names]
box = ax.boxplot([samples[s]["mean"].clip(lower=0.1) for s in names], whis=(5, 95), showfliers=False,
                 patch_artist=True, widths=0.6, medianprops={"color": TEXT, "lw": 1.2})
for patch, t in zip(box["boxes"], types):
    patch.set(facecolor=TYPE_COLOUR[t], alpha=0.6, edgecolor=MUTED)
style(ax)
ax.set_xticks(range(1, len(names) + 1),
              [f"{s} ({meta.loc[s, 'Country']})" for s in names], fontsize=8, rotation=45, ha="right")
ax.set_ylabel("Mean depth per 500 bp window", fontsize=9, color=TEXT)
ax.legend(handles=[Patch(color=c, alpha=0.6, label=t) for t, c in TYPE_COLOUR.items()],
          frameon=False, fontsize=9, loc="upper right")
fig.suptitle("Depth across the genome in each sample (box = middle 50% of windows, whiskers = 5-95%)",
             fontsize=11, color=TEXT, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(os.path.join(out_dir, "depth_boxplot.png"), dpi=150, facecolor="white")
plt.close(fig)

# --- per sample type: mean and +/- SD along chr7, and around crt ---
for sample_type in TYPE_COLOUR:
    group = [s for s in samples if meta.loc[s, "Sample type"] == sample_type]
    if not group:
        continue
    fig, axes = plt.subplots(2, 1, figsize=(12, 8.5))
    for ax, (lo, hi), where in [(axes[0], (0, 1.45e6), "chr7 (Pf3D7_07_v3)"),
                                (axes[1], (350_000, 450_000), "crt, Pf3D7_07_v3:350,000-450,000")]:
        for s, colour in zip(group, SAMPLE_COLOURS):
            c = samples[s]
            c = c[(c["chrom"] == "Pf3D7_07_v3") & (c["end"] > lo) & (c["start"] < hi)]
            mid_mb = (c["start"] + c["end"]) / 2 / 1e6
            # depth 0 is drawn at 0.5, the bottom of the log axis
            ax.fill_between(mid_mb, (c["mean"] - c["sd"]).clip(lower=0.5), (c["mean"] + c["sd"]).clip(lower=0.5),
                            color=colour, alpha=0.12, lw=0)
            ax.plot(mid_mb, c["mean"].clip(lower=0.5), color=colour, lw=0.7,
                    label=f"{s} ({meta.loc[s, 'Country']}, median {samples[s]['mean'].median():.0f}×)")
        style(ax)
        ax.set_ylim(0.5, None)
        ax.set_xlim(lo / 1e6, hi / 1e6)
        ax.set_ylabel("Depth per 500 bp window", fontsize=9, color=TEXT)
        ax.set_title(f"Along {where}", loc="left", fontsize=10, color=TEXT)
    axes[1].set_xlabel("Position on chromosome (Mb)", fontsize=9, color=TEXT)
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0.01, 0.955), ncol=4)
    fig.suptitle(f"{sample_type} samples: read depth (line = window mean, band = ± SD)",
                 fontsize=12, color=TEXT, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(os.path.join(out_dir, f"depth_{sample_type}.png"), dpi=150, facecolor="white")
    plt.close(fig)

print(f"wrote depth_boxplot.png and depth_<type>.png in {out_dir}")
