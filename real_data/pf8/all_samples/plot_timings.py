"""
How long the Pf8 depth pipeline takes per sample, and an extrapolation to every Pf8 sample.

Reads $PF8/timings/*.tsv (from pf8_depth_task.sh) and $PF8/input/Pf8_crams.tsv (from pf8_crams.py).
For each stage, fits seconds = a + b * CRAM GB on the timed samples (or only those from one array job,
since download speed depends on how many tasks ran at once), then applies it to every Pf8 CRAM.

Writes in $PF8/figures:
  timing_vs_cram_size.png   download, mosdepth, window stats and total time against CRAM size; peak memory
  pf8_cram_sizes.png        CRAM sizes of all Pf8 samples, with the timed samples marked
  pf8_extrapolation.png     days to do every Pf8 sample against how many array tasks run at once
  timing_fits.tsv           the fitted a and b for each stage
  pf8_extrapolation.tsv     predicted hours per stage, summed over every Pf8 sample

Usage: python plot_timings.py [array job ID to fit on, default all samples]
"""

import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Folder for the data (change this on another computer)
PF8 = "/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth"
FIGS = f"{PF8}/figures"
os.makedirs(FIGS, exist_ok=True)

COLOURS = {"sWGA": "#2a78d6", "gDNA": "#eb6834"}
FIT = "#52514e"
CONCURRENT = [10, 25, 50, 100, 200]   # array tasks at once, for the extrapolation
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e4e3df", "grid.linewidth": 0.6, "axes.axisbelow": True})

# one row per timed sample, one column per stage
long = pd.concat([d for f in glob.glob(f"{PF8}/timings/*.tsv")
                  if len(d := pd.read_csv(f, sep="\t"))])   # tasks still starting have only the header
FIT_JOB = int(sys.argv[1]) if len(sys.argv) > 1 else None
if "job" not in long:
    long["job"] = np.nan
long["fit"] = True if FIT_JOB is None else long.job.eq(FIT_JOB)   # the first test has no job column (NaN)
failed = long[long.exit != 0]
if len(failed):
    print("stages that failed (left out):\n", failed[["sample", "stage", "exit"]].to_string(index=False))
long = long[long.exit == 0]
long["stage"] = long.stage.replace({"download_crai": "download"})   # the .crai is part of the download
times = long.pivot_table(index=["sample", "type", "cram_bytes", "fit"], columns="stage", values="seconds", aggfunc="sum")
times = times.dropna(subset=["download", "mosdepth", "window_stats"]).reset_index()   # only samples where every stage worked
# later samples make the .crai with samtools (stage "index") instead of downloading it
times["total"] = times[[c for c in ["download", "index", "mosdepth", "window_stats", "cleanup"] if c in times]].sum(axis=1)
times["GB"] = times.cram_bytes / 1e9
rss = long.pivot_table(index="sample", columns="stage", values="max_rss_gb", aggfunc="max")
print(f"{len(times)} samples timed, {times.fit.sum()} used for the fits")
fitted = times[times.fit]

# seconds = a + b * GB, for each stage
STAGES = ["download", "mosdepth", "window_stats", "total"]
fits = pd.DataFrame({s: np.polyfit(fitted.GB, fitted[s], 1)[::-1] for s in STAGES}, index=["a_seconds", "b_seconds_per_GB"]).T
fits["r2"] = [np.corrcoef(fitted.GB, fitted[s])[0, 1] ** 2 for s in STAGES]
fits["MB_per_s"] = np.where(fits.index == "download", 1000 / fits.b_seconds_per_GB, np.nan)
fits.round(3).to_csv(f"{FIGS}/timing_fits.tsv", sep="\t")
print(fits.round(2))

# figure 1: time against CRAM size, one panel per stage, plus peak memory
titles = {"download": "Download (CRAM + index)", "mosdepth": "mosdepth", "window_stats": "Window mean and SD",
          "total": "Whole pipeline"}
fig, axes = plt.subplots(1, 5, figsize=(19, 3.8), constrained_layout=True)
grid = np.linspace(0, times.GB.max() * 1.05, 50)
for ax, s in zip(axes, STAGES):
    for t, c in COLOURS.items():
        d = times[times.type == t]
        ax.scatter(d.GB[d.fit], d[s][d.fit] / 60, s=40, color=c, edgecolor="white", linewidth=1.5, label=t, zorder=3)
        ax.scatter(d.GB[~d.fit], d[s][~d.fit] / 60, s=40, facecolor="none", edgecolor=c, linewidth=1.2, zorder=3,
                   label=f"{t}, not in fit" if (~times.fit).any() else None)
    a, b = fits.loc[s, ["a_seconds", "b_seconds_per_GB"]]
    ax.plot(grid, (a + b * grid) / 60, color=FIT, lw=1.5, ls="--", label="linear fit")
    ax.set_title(f"{titles[s]}\n{b / 60:.2f} min per GB, r² = {fits.loc[s, 'r2']:.2f}", fontsize=10, loc="left")
    ax.set_xlabel("CRAM size (GB)")
    ax.set_ylim(bottom=0)
axes[0].set_ylabel("minutes")
axes[0].legend(frameon=False, fontsize=9)
ax = axes[4]
for s, m in [("mosdepth", "o"), ("window_stats", "s")]:
    for t, c in COLOURS.items():
        d = times[times.type == t]
        ax.scatter(d.GB, rss.loc[d["sample"], s], s=40, marker=m, color=c, edgecolor="white", linewidth=1.5, zorder=3)
ax.scatter([], [], marker="o", color=FIT, label="mosdepth")
ax.scatter([], [], marker="s", color=FIT, label="window stats")
ax.legend(frameon=False, fontsize=9)
ax.set_title("Peak memory", fontsize=10, loc="left")
ax.set_xlabel("CRAM size (GB)")
ax.set_ylabel("GB")
ax.set_ylim(bottom=0)
fig.savefig(f"{FIGS}/timing_vs_cram_size.png", dpi=150)

# every Pf8 CRAM, and the predicted hours to do them all
crams = pd.read_csv(f"{PF8}/input/Pf8_crams.tsv", sep="\t")
crams["GB"] = crams.cram_bytes / 1e9
rows = []
for name, d in [("all", crams), ("QC pass", crams[crams.qc_pass])]:
    for t in ["all", "sWGA", "gDNA"]:
        dd = d if t == "all" else d[d.type == t]
        row = {"samples": name, "type": t, "n": len(dd), "TB": dd.GB.sum() / 1000}
        for s in STAGES:
            a, b = fits.loc[s, ["a_seconds", "b_seconds_per_GB"]]
            row[f"{s}_hours"] = (len(dd) * a + b * dd.GB.sum()) / 3600
        rows.append(row)
extra = pd.DataFrame(rows)
for k in CONCURRENT:
    extra[f"days_{k}_at_once"] = extra.total_hours / 24 / k
extra.round(2).to_csv(f"{FIGS}/pf8_extrapolation.tsv", sep="\t", index=False)
print(extra.round(1).to_string(index=False))

# figure 2: CRAM sizes of all Pf8 samples, timed ones marked
fig, ax = plt.subplots(figsize=(8, 3.8), constrained_layout=True)
bins = np.linspace(0, crams.GB.quantile(0.999), 60)
for t, c in COLOURS.items():
    ax.hist(crams.GB[crams.type == t], bins=bins, color=c, alpha=0.55, label=f"{t} (n = {(crams.type == t).sum():,})")
    ax.plot(times.GB[times.type == t], np.zeros((times.type == t).sum()), "|", color=c, ms=18, mew=2.5,
            clip_on=False, zorder=4)
ax.set_xlabel("CRAM size (GB)")
ax.set_ylabel("Pf8 samples")
ax.set_title("CRAM sizes of all Pf8 samples (ticks: the timed samples)", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=9)
fig.savefig(f"{FIGS}/pf8_cram_sizes.png", dpi=150)

# figure 3: days for every Pf8 sample, against tasks at once
fig, ax = plt.subplots(figsize=(7, 3.8), constrained_layout=True)
k = np.array(CONCURRENT)
for name, colour in [("all", "#1baf7a"), ("QC pass", "#6250d6")]:   # not blue/orange: those mean sWGA/gDNA
    r = extra[(extra.samples == name) & (extra.type == "all")].iloc[0]
    ax.plot(k, r.total_hours / 24 / k, color=colour, lw=2, marker="o", ms=6,
            label=f"{name} samples ({r.n:,}, {r.total_hours:,.0f} task-hours)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks(k, [str(x) for x in k])
ax.set_xlabel("array tasks running at once")
ax.set_ylabel("days")
ax.set_title("Predicted time for every Pf8 sample\n(assumes download speed per task holds at that many at once)",
             fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=9)
fig.savefig(f"{FIGS}/pf8_extrapolation.png", dpi=150)
print(f"figures in {FIGS}")
