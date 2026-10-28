"""
Plots that make the results of compare_processes.py easier to read.

    python plot_processes.py <pf8_dir> <output_dir>

<output_dir> is the folder compare_processes.py wrote fits.tsv into. Writes there:
  ranking.png    models from best to worst, sWGA and gDNA side by side (median held-out score)
  ranks.png      every sample x model: its rank under the held-out score (1 = best for that sample)
  examples.png   every sample: real autocorrelation over the first 5 kb next to AR(1) and the best model of
                 its type, ordered by how far AR(1) is behind
  ar_lookback.png  AR(p): median held-out score against how far back it looks (p windows = p x 0.5 kb)
  arma_gp.png      ARMA and GP models: median held-out score, sWGA and gDNA
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
from scipy.signal import lfilter

from compare_processes import K, WINDOW, fit, lag_matrix, load, models

TYPES = {"sWGA": "#2a78d6", "gDNA": "#eb6834"}
TEXT, MUTED, GRID, GREY = "#0b0b0b", "#8a8984", "#e6e5e1", "#b9b8b3"
SCORE = "dHeldout"


def style(ax, grid="y"):
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8, labelcolor=TEXT)


def ranking(res, out_png):
    """Models from best to worst for each sample type, by median held-out score."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5))
    for ax, (t, colour) in zip(axes, TYPES.items()):
        g = res[res["sample type"] == t]
        med = g.groupby("model")[SCORE].median().sort_values(ascending=False)
        wins = g[g[SCORE] == 0]["model"].value_counts()
        colours = [colour if i == len(med) - 1 else GREY for i in range(len(med))]
        ax.barh(range(len(med)), med + 1, color=colours)
        ax.set_xscale("log")
        ax.set_yticks(range(len(med)), [f"{m}  ({wins.get(m, 0)} wins)" if wins.get(m, 0) else m for m in med.index],
                      fontsize=8)
        ax.set_xlabel("Median distance from each sample's best model (−2 × held-out LL, + 1, log scale)", fontsize=8)
        ax.set_title(f"{t} ({g['sample'].nunique()} samples): lowest median = {med.index[-1]}",
                     loc="left", fontsize=11, color=TEXT)
        style(ax, grid="x")
    fig.suptitle("Which process fits read depth best? Shorter bar = better; '(n wins)' = best for n samples",
                 x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def ranks(res, out_png):
    """Heatmap: rank of every model in every sample (1 = best), samples grouped by type."""
    order = res.groupby("model")[f"rank Heldout"].median().sort_values().index
    samples = res.sort_values(["sample type", "sample"], ascending=[False, True])["sample"].unique()
    grid = res.pivot(index="sample", columns="model", values="rank Heldout").loc[samples, order]
    fig, ax = plt.subplots(figsize=(13, 0.38 * len(samples) + 2.5))
    im = ax.imshow(grid, cmap="Blues_r", aspect="auto", vmin=1, vmax=len(order))
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            r = int(grid.iat[i, j])
            ax.text(j, i, r, ha="center", va="center", fontsize=7, color="white" if r <= 6 else TEXT)
    t_of = res.drop_duplicates("sample").set_index("sample")["sample type"]
    ax.set_yticks(range(len(samples)), [f"{s} ({t_of[s]})" for s in samples], fontsize=8)
    ax.set_xticks(range(len(order)), order, rotation=45, ha="right", fontsize=8)
    n_swga = (t_of[samples] == "sWGA").sum()
    ax.axhline(n_swga - 0.5, color=TEXT, lw=1.5)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label="rank (1 = best)")
    ax.set_title("Rank of each model in each sample, by held-out likelihood (models ordered by median rank)",
                 loc="left", fontsize=11, color=TEXT)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def model_acf(build, x, max_lag, n_psi=5000):
    """Autocorrelation implied by a fitted model, from its best predictor on K previous windows."""
    phi, v = build(x)
    impulse = np.zeros(n_psi)
    impulse[0] = 1
    psi = lfilter([1], np.r_[1, -phi[K, 1:]], impulse)
    gamma = np.array([psi[:n_psi - k] @ psi[k:] for k in range(max_lag + 1)])
    return gamma / gamma[0]


def acf(ys, max_lag):
    out = []
    for y in ys:
        y = y - y.mean()
        out.append([1.0] + [(y[:-k] * y[k:]).sum() / (y * y).sum() for k in range(1, max_lag + 1)])
    return np.mean(out, axis=0)


def examples(res, pf8_dir, out_png, max_lag=10, ncol=4):
    """Every sample: real autocorrelation over the first 5 kb vs AR(1) and the type's best model."""
    all_models = models()
    lag_kb = np.arange(max_lag + 1) * WINDOW / 1e3
    best = {t: res[res["sample type"] == t].groupby("model")[SCORE].median().idxmin() for t in TYPES}
    rows_per_type = {t: int(np.ceil(res.loc[res["sample type"] == t, "sample"].nunique() / ncol)) for t in TYPES}
    fig, axes = plt.subplots(sum(rows_per_type.values()), ncol, figsize=(15, 3.1 * sum(rows_per_type.values())))
    row0 = 0
    for t, colour in TYPES.items():
        g = res[res["sample type"] == t]
        ar1_gap = g[g["model"] == "AR(1)"].set_index("sample")[SCORE].sort_values()   # AR(1) closest to best first
        for i, sample in enumerate(ar1_gap.index):
            ax = axes[row0 + i // ncol, i % ncol]
            ys = load(os.path.join(pf8_dir, sample, "depth_plots", f"{sample}.window_stats.tsv"))
            Ls = [lag_matrix(y) for y in ys]
            ax.plot(lag_kb, acf(ys, max_lag), color=TEXT, lw=2.4, label="real")
            for name, c, ls in (("AR(1)", GREY, "--"), (best[t], colour, "-")):
                build, start, _ = all_models[name]
                x, _ = fit(ys, Ls, build, start(ys))
                ax.plot(lag_kb, model_acf(build, x, max_lag), color=c, lw=1.6, ls=ls, label=name)
            ax.set_title(f"{sample} ({t})\nAR(1) behind its best model by {ar1_gap[sample]:,.0f}",
                         loc="left", fontsize=8.5, color=TEXT)
            ax.set_xticks(range(6))
            style(ax)
            if i == 0:
                ax.legend(frameon=False, fontsize=7.5)
        for j in range(len(ar1_gap), rows_per_type[t] * ncol):
            axes[row0 + j // ncol, j % ncol].set_visible(False)
        row0 += rows_per_type[t]
    fig.supxlabel("Distance (kb)", fontsize=10)
    fig.supylabel("How similar windows are (autocorrelation)", fontsize=10)
    fig.suptitle("Real autocorrelation (black) vs AR(1) (grey dashed) and each type's best model, first 5 kb of distance, every sample\n"
                 "Panels ordered from AR(1) closest to best to furthest; each panel has its own y-axis, so small differences show",
                 x=0.01, ha="left", fontsize=11, color=TEXT)
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.97))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)

def ar_lookback(res, out_png):
    """Median held-out score of AR(p) against how far back it looks, one line per sample type."""
    ar = res[res["model"].str.match(r"AR\(")].copy()
    ar["p"] = ar["model"].str.extract(r"(\d+)", expand=False).astype(int)
    fig, ax = plt.subplots(figsize=(9, 5))
    for t, colour in TYPES.items():
        g = ar[ar["sample type"] == t]
        med = g.groupby("p")[SCORE].median()
        for _, row in g.iterrows():
            ax.scatter(row["p"] * WINDOW / 1e3, row[SCORE] + 1, s=8, color=colour, alpha=0.25, lw=0)
        ax.plot(med.index * WINDOW / 1e3, med + 1, "o-", color=colour, lw=2, ms=5,
                label=f"{t} (median of {g['sample'].nunique()} samples; dots = samples)")
        best_p = med.idxmin()
        ax.annotate(f"lowest: AR({best_p})", (best_p * WINDOW / 1e3, med[best_p] + 1), xytext=(0, -16),
                    textcoords="offset points", ha="center", fontsize=8, color=colour)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,g}"))
    ax.set_xticks([p * WINDOW / 1e3 for p in sorted(ar["p"].unique())],
                  [f"{p * WINDOW / 1e3:g}\nAR({p})" for p in sorted(ar["p"].unique())], fontsize=7.5)
    ax.set_xlabel("How far back the model looks (kb)", fontsize=9)
    ax.set_ylabel("Held-out score − best model of the sample, + 1 (log scale; lower = better)", fontsize=9)
    ax.set_title("AR(p): how many previous depths do you need to look back at?\n"
                 "sWGA levels off at about 3–4 kb; gDNA stays behind its best model (ARMA(2,2)) however far back AR looks",
                 loc="left", fontsize=10, color=TEXT)
    ax.legend(frameon=False, fontsize=8)
    style(ax)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def arma_gp(res, out_png):
    """Median held-out score of the ARMA and GP models, sWGA and gDNA side by side."""
    groups = [("ARMA", ["ARMA(1,1)", "ARMA(2,1)", "ARMA(2,2)"]),
              ("Gaussian process", ["GP exponential", "GP Matern 3/2", "GP Matern 5/2", "GP squared exp", "GP two exponentials"])]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [3, 5]}, sharey=True)
    for ax, (family, names) in zip(axes, groups):
        x = np.arange(len(names))
        for i, (t, colour) in enumerate(TYPES.items()):
            g = res[res["sample type"] == t]
            med = g.groupby("model")[SCORE].median().reindex(names)
            bars = ax.bar(x + (i - 0.5) * 0.38, med + 1, width=0.36, color=colour, label=t)
            for b, v in zip(bars, med):
                ax.text(b.get_x() + b.get_width() / 2, (v + 1) * 1.12, f"{v:,.0f}", ha="center", fontsize=7.5, color=TEXT)
        ax.set_xticks(x, [n.replace("GP ", "") for n in names], fontsize=8.5)
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,g}"))
        ax.set_title(family, loc="left", fontsize=11, color=TEXT)
        style(ax)
    axes[0].set_ylabel("Median held-out score − best model of the sample, + 1\n(log scale; lower = better)", fontsize=9)
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("ARMA and Gaussian-process models: ARMA(2,2) is best for gDNA (0); the smoothest GP kernel is worst for both",
                 x=0.01, ha="left", fontsize=11, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pf8_dir")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    res = pd.read_csv(os.path.join(args.output_dir, "fits.tsv"), sep="\t")
    ranking(res, os.path.join(args.output_dir, "ranking.png"))
    ranks(res, os.path.join(args.output_dir, "ranks.png"))
    examples(res, args.pf8_dir, os.path.join(args.output_dir, "examples.png"))
    ar_lookback(res, os.path.join(args.output_dir, "ar_lookback.png"))
    arma_gp(res, os.path.join(args.output_dir, "arma_gp.png"))
    print(f"wrote ranking.png, ranks.png, examples.png, ar_lookback.png and arma_gp.png in {args.output_dir}")


if __name__ == "__main__":
    main()
