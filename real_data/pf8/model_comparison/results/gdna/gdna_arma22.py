"""gDNA: what the fitted ARMA(2,2) says about each sample.

ARMA(2,2) = a fast AR(1) level + a slow AR(1) level + independent noise on every window.
For each gDNA sample: fit ARMA(2,2), read the two fading speeds from its AR part, and split the
variance of log depth into fast, slow and noise by matching the fitted autocovariance.
Writes gdna_arma22.tsv and gdna_arma22.png next to this script.
"""
import os
import sys

import numpy as np
import pandas as pd

REPO = "/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/TSKIT_MALARIA_SIM/git/KR_simulator/real_data/pf8/model_comparison"
sys.path.insert(0, REPO)
from compare_processes import WINDOW, arma_gamma, fit, lag_matrix, load, models, pacf_to_coefs  # noqa: E402

PF8 = "/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/pf8_samples"
HERE = os.path.dirname(os.path.abspath(__file__))

info = pd.read_csv(f"{PF8}/samples.tsv", sep="\t").set_index("Sample")
fits = pd.read_csv(f"{PF8}/model/process_comparison/fits.tsv", sep="\t")
gdna = [s for s in info.index if info.loc[s, "Sample type"] == "gDNA"]
build, start, _ = models()["ARMA(2,2)"]

rows, acfs = [], {}
for s in gdna:
    ys = load(f"{PF8}/{s}/depth_plots/{s}.window_stats.tsv")
    Ls = [lag_matrix(y) for y in ys]
    x, _ = fit(ys, Ls, build, start(ys))
    # the two fading speeds: roots of the AR polynomial 1 - c1 z - c2 z^2  ->  phi = 1 / root
    c = pacf_to_coefs(np.tanh(x[1:3]))
    phis = np.sort(np.real_if_close(np.roots([1, -c[0], -c[1]])))     # phi_fast < phi_slow
    lags = np.arange(201)
    g = arma_gamma(x, 2, 2, lags)
    # gamma(k) = s_fast phi_fast^k + s_slow phi_slow^k for k >= 1; noise = what is left at lag 0
    A = np.column_stack([phis[0] ** lags[1:], phis[1] ** lags[1:]])
    (s_fast, s_slow), *_ = np.linalg.lstsq(A, g[1:], rcond=None)
    noise = g[0] - s_fast - s_slow
    length_kb = lambda phi: -WINDOW / 1e3 / np.log(phi) if 0 < phi < 1 else np.nan
    f = fits[fits["sample"] == s].set_index("model")
    real = np.concatenate(ys)
    rows.append({
        "sample": s, "study": info.loc[s, "Study"], "country": info.loc[s, "Country"],
        "median depth": round(float(pd.read_csv(f"{PF8}/{s}/depth_plots/{s}.window_stats.tsv", sep="\t")
                                     .query("~chrom.str.contains('API|MIT')", engine="python")["mean"].median())),
        "best model": f["dHeldout"].idxmin(),
        "ARMA(2,2) rank": int(f.loc["ARMA(2,2)", "rank Heldout"]),
        "AR(1) behind": round(f.loc["AR(1)", "dHeldout"]),
        "AR(20) behind": round(f.loc["AR(20)", "dHeldout"]),
        "var log depth": round(g[0], 3),
        "fast length (kb)": round(length_kb(phis[0]), 2),
        "slow length (kb)": round(length_kb(phis[1]), 1),
        "fast share": round(s_fast / g[0], 2),
        "slow share": round(s_slow / g[0], 2),
        "noise share": round(noise / g[0], 2),
    })
    acfs[s] = (g[:61] / g[0], np.mean([[1] + [np.corrcoef(y[:-k], y[k:])[0, 1] for k in range(1, 61)] for y in ys], axis=0))
    print(s, "done", flush=True)

tab = pd.DataFrame(rows).sort_values("median depth")
tab.to_csv(f"{HERE}/gdna_arma22.tsv", sep="\t", index=False)
with pd.option_context("display.width", 250, "display.max_columns", 30):
    print(tab.to_string(index=False))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
TEXT, MUTED, GRID = "#0b0b0b", "#8a8984", "#e6e5e1"
FAST, SLOW, NOISE = "#2a78d6", "#1baf7a", "#b9b8b3"
fig, axes = plt.subplots(1, 2, figsize=(14, 5.8), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axes[0]
ok = tab[~tab["study"].str.contains("TRAC2")]
yl = [f"{r['sample']} ({r['median depth']}×)" for _, r in ok.iterrows()]
y = np.arange(len(ok))
h = 0.26
ax.barh(y + h, ok["fast share"], height=h, color=FAST, label="fast part")
ax.barh(y, ok["slow share"], height=h, color=SLOW, label="slow part")
ax.barh(y - h, ok["noise share"], height=h, color=NOISE, label="per-window noise (negative = none; see text)")
for i, (_, r) in enumerate(ok.iterrows()):
    ax.text(1.32, i, f"fast {r['fast length (kb)']:.1f} kb · slow {r['slow length (kb)']:.0f} kb", va="center", fontsize=8, color=TEXT)
ax.axvline(0, color=TEXT, lw=0.8)
ax.set_yticks(y, yl, fontsize=8)
ax.set_xlim(-0.45, 1.9)
ax.set_xlabel("share of the variation in log depth (the three add up to 1)", fontsize=9)
ax.set_title("ARMA(2,2) split into parts, typical gDNA samples (TRAC2 samples left out: they do not split this way)\n"
             "right: how far the fast and slow ups and downs last", loc="left", fontsize=10, color=TEXT)
ax.legend(frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(1.0, -0.22), ncol=3)
ax = axes[1]
lag_kb = np.arange(61) * WINDOW / 1e3
for s, (model_acf, real_acf) in acfs.items():
    ax.plot(lag_kb, real_acf, color=TEXT, lw=0.8, alpha=0.5)
    ax.plot(lag_kb, model_acf, color="#eb6834", lw=0.8, alpha=0.7)
ax.plot([], [], color=TEXT, label="real (one line per sample)")
ax.plot([], [], color="#eb6834", label="fitted ARMA(2,2)")
ax.set_xlim(0, 15)
ax.set_xlabel("distance (kb)", fontsize=9)
ax.set_ylabel("autocorrelation", fontsize=9)
ax.set_title("How similar windows are, by distance: real vs fitted", loc="left", fontsize=10, color=TEXT)
ax.legend(frameon=False, fontsize=8)
for ax in axes:
    ax.grid(axis="x" if ax is axes[0] else "y", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
fig.tight_layout()
fig.savefig(f"{HERE}/gdna_arma22.png", dpi=150, facecolor="white", bbox_inches="tight")
print("wrote", f"{HERE}/gdna_arma22.png")
