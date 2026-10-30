"""Explore the structure of read depth along the genome, to decide which depth models to test.

Works on log depth per 500 bp window, y = log(mean + 1), for the samples of one type, and splits it into
    y[sample, window] = sample level + shared window effect + sample-specific residual
  sample level         = the sample's mean of y
  shared window effect = the mean across samples of (y - sample level), per window (means, so the parts add up)
  residual             = what is left; its own patchiness along the genome

Writes <type>_explore.png and <type>_explore.txt into this folder.

Usage:
    python explore_depth.py sWGA
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
REF = os.path.join(TOP, "reference", "PlasmoDB-54-Pfalciparum3D7-Genome.fasta")
WINDOW = 500
MAX_LAG = 200  # windows (100 kb)

sample_type = sys.argv[1]
meta = pd.read_csv(os.path.join(TOP, "samples.tsv"), sep="\t")
samples = list(meta.loc[meta["Sample type"] == sample_type, "Sample"])
samples = [s for s in samples if os.path.exists(f"{TOP}/{s}/depth_plots/{s}.window_stats.tsv")]

# --- one matrix: rows = samples, columns = nuclear 500 bp windows ---
frames = []
for s in samples:
    df = pd.read_csv(f"{TOP}/{s}/depth_plots/{s}.window_stats.tsv", sep="\t")
    frames.append(df[~df["chrom"].str.contains("API|MIT")].reset_index(drop=True))
windows = frames[0][["chrom", "start", "end"]]
mean = np.vstack([f["mean"].to_numpy() for f in frames])
sd = np.vstack([f["sd"].to_numpy() for f in frames])

y = np.log(mean + 1)
level = y.mean(axis=1, keepdims=True)
shared = (y - level).mean(axis=0)
resid = y - level - shared

# --- GC content of each window, from the reference ---
seqs, name = {}, None
for line in open(REF):
    if line.startswith(">"):
        name = line[1:].split()[0]
        seqs[name] = []
    else:
        seqs[name].append(line.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
gc = np.array([(lambda w: (w.count("G") + w.count("C")) / max(len(w), 1))(seqs[c][s:e])
               for c, s, e in windows.itertuples(index=False)])


def acf(x, chrom, max_lag):
    """Autocorrelation at lags 1..max_lag, using only pairs of windows on the same chromosome."""
    x = x - x.mean()
    out = []
    for lag in range(1, max_lag + 1):
        same = chrom[lag:] == chrom[:-lag]
        a, b = x[lag:][same], x[:-lag][same]
        out.append((a * b).mean() / x.var())
    return np.array(out)


chrom = windows["chrom"].to_numpy()
acfs = np.array([acf(r, chrom, MAX_LAG) for r in resid])
acf_shared = acf(shared, chrom, MAX_LAG)

# --- numbers ---
var_total = y.var(axis=1).mean()
var_shared = shared.var()
var_resid = resid.var(axis=1).mean()
kurt = [pd.Series(r).kurt() for r in resid]
pos = mean > 0
disp = np.median((sd[pos] ** 2) / mean[pos])
lines = [
    f"{sample_type}: {len(samples)} samples, {y.shape[1]:,} nuclear 500 bp windows",
    "",
    "median depth: " + ", ".join(f"{s} {np.median(m):.0f}x" for s, m in zip(samples, mean)),
    "",
    "variance of log depth within a sample, split into parts:",
    f"  total {var_total:.3f} = shared across samples {var_shared:.3f} ({var_shared / var_total:.0%})"
    f" + sample-specific {var_resid:.3f} ({var_resid / var_total:.0%})",
    f"  correlation of GC content with the shared effect: {np.corrcoef(gc, shared)[0, 1]:.2f}",
    "",
    "autocorrelation of the sample-specific residual (median over samples):",
    "  lag (windows):  " + "  ".join(f"{l:>5}" for l in (1, 2, 5, 10, 20, 50, 100, 200)),
    "  acf:            " + "  ".join(f"{np.median(acfs[:, l - 1]):5.2f}" for l in (1, 2, 5, 10, 20, 50, 100, 200)),
    "  AR(1) would give:" + "  ".join(f"{np.median(acfs[:, 0]) ** l:5.2f}" for l in (1, 2, 5, 10, 20, 50, 100, 200)),
    "autocorrelation of the shared effect:",
    "  acf:            " + "  ".join(f"{acf_shared[l - 1]:5.2f}" for l in (1, 2, 5, 10, 20, 50, 100, 200)),
    "",
    f"residual excess kurtosis per sample (0 = Gaussian): " + ", ".join(f"{k:.1f}" for k in kurt),
    f"windows with zero depth: {np.mean(mean == 0):.2%}",
    f"within-window variance / mean (median over windows): {disp:.1f}  (1 = Poisson)",
]
txt = "\n".join(lines)
open(os.path.join(HERE, f"{sample_type}_explore.txt"), "w").write(txt + "\n")
print(txt)

# --- plots ---
TEXT, MUTED, GRID, BLUE, ORANGE = "#0b0b0b", "#8a8984", "#e6e5e1", "#2a78d6", "#eb6834"
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

ax = axes[0, 0]
lags = np.arange(1, MAX_LAG + 1)
for a in acfs:
    ax.plot(lags * WINDOW / 1000, a, color=BLUE, lw=0.8, alpha=0.5)
ax.plot(lags * WINDOW / 1000, np.median(acfs[:, 0]) ** lags, color=TEXT, ls="--", lw=1.2, label="AR(1) with the lag-1 value")
ax.plot(lags * WINDOW / 1000, acf_shared, color=ORANGE, lw=1.5, label="shared effect")
ax.set_xscale("log")
ax.set_xlabel("Distance (kb)")
ax.set_ylabel("Autocorrelation")
ax.set_title("Autocorrelation of sample-specific log depth (blue, one per sample)", loc="left", fontsize=10)
ax.legend(frameon=False, fontsize=8)

ax = axes[0, 1]
i, j = 0, 1
ax.scatter(y[i] - level[i], y[j] - level[j], s=2, alpha=0.1, color=BLUE, lw=0, rasterized=True)
ax.set_xlabel(f"{samples[i]}: log depth - sample level")
ax.set_ylabel(f"{samples[j]}: log depth - sample level")
r = np.corrcoef(y[i], y[j])[0, 1]
ax.set_title(f"Two samples, window by window (r = {r:.2f}): how much is shared", loc="left", fontsize=10)

ax = axes[1, 0]
ax.scatter(gc * 100, shared, s=2, alpha=0.1, color=BLUE, lw=0, rasterized=True)
ax.set_xlabel("GC content of window (%)")
ax.set_ylabel("Shared window effect (log)")
ax.set_title(f"Shared effect against GC content (r = {np.corrcoef(gc, shared)[0, 1]:.2f})", loc="left", fontsize=10)

ax = axes[1, 1]
z = (resid - resid.mean(axis=1, keepdims=True)) / resid.std(axis=1, keepdims=True)
bins = np.linspace(-6, 6, 121)
ax.hist(z.ravel(), bins=bins, density=True, color=BLUE, alpha=0.6, label="residuals, all samples")
ax.plot(bins, np.exp(-bins ** 2 / 2) / np.sqrt(2 * np.pi), color=TEXT, ls="--", lw=1.2, label="Gaussian")
ax.set_yscale("log")
ax.set_ylim(1e-5, 1)
ax.set_xlabel("Standardised sample-specific residual")
ax.set_title("Tails of the residual compared with a Gaussian", loc="left", fontsize=10)
ax.legend(frameon=False, fontsize=8)

for ax in axes.flat:
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
fig.suptitle(f"{sample_type} ({len(samples)} samples): structure of log depth in 500 bp windows", x=0.01, ha="left", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(os.path.join(HERE, f"{sample_type}_explore.png"), dpi=130, facecolor="white")
print("wrote", os.path.join(HERE, f"{sample_type}_explore.png"))
