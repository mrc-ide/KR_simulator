"""
Figures for fit_chosen_models.py: the chosen depth models over every QC-pass Pf8 gDNA and sWGA sample.

    python plot_chosen_models.py          # reads the tables in OUT (fit_chosen_models.py), writes OUT/figures/

1_models_sWGA.png, 2_models_gDNA.png   the four models for one sample type: how far behind the best each one is
                                       in every sample, how often it wins, and whether its extreme kicks are as
                                       common as it says
3_selection_criteria.png               which model AIC, BIC and held-out likelihood pick, and how often they agree
4_all_samples.png                      the chosen model over all samples: tail weight (nu), tail calibration,
                                       one shared parameter set vs each sample's own fit, and by population
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

from best_simple_model import label
from fit_chosen_models import CHOSEN, OUT

TEXT, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
CHOSEN_COL, OTHER_T, GAUSS = "#2a78d6", "#7d7c78", "#c9c8c3"
TYPE_COL = {"gDNA": "#2a78d6", "sWGA": "#eb6834"}
CRIT_COL = {"AIC": "#2a78d6", "BIC": "#eb6834", "Held-out": "#1baf7a"}
MODELS = ["ARMA(2,2) t", "AR(8) t", "ARMA(2,2) Gaussian", "AR(8) Gaussian"]
DROP, SPIKE = "lowest 0.1% (expect 0.1)", "highest 0.1% (expect 0.1)"
comma = FuncFormatter(lambda v, _: f"{v:,g}")


def chosen_name(t):
    (p, q), kick = CHOSEN[t]
    return f"{label(p, q)} {kick}"


def model_order(t):
    c = chosen_name(t)
    return [c] + [m for m in MODELS if m != c]


def model_col(t, m):
    return CHOSEN_COL if m == chosen_name(t) else GAUSS if m.endswith("Gaussian") else OTHER_T


def tidy(ax, grid="x"):
    ax.grid(axis=grid, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#b9b8b3")
    ax.tick_params(colors=MUTED, labelsize=8.5)


def title(ax, text, sub=None):
    ax.set_title(text, loc="left", fontsize=10.5, color=TEXT, pad=18 if sub else 6)
    if sub:
        ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=8.5, color=MUTED, va="bottom")


def ranges(ax, ys, data, cols, log=True):
    """10-90% line, 25-75% bar and median dot for each row."""
    for y, v, c in zip(ys, data, cols):
        v = np.asarray(v)
        v = v[np.isfinite(v)]
        lo, q1, med, q3, hi = np.percentile(v, [10, 25, 50, 75, 90])
        ax.plot([lo, hi], [y, y], color=c, lw=2, solid_capstyle="round")
        ax.plot([q1, q3], [y, y], color=c, lw=7, solid_capstyle="round")
        ax.plot(med, y, "o", ms=8, color="white", mec=c, mew=2, zorder=3)
    if log:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(comma)


def plot_type(res, t, path):
    g = res[res["sample type"] == t]
    order = model_order(t)
    n = g["sample"].nunique()
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), gridspec_kw={"width_ratios": [1.5, 0.8, 1.1]})
    ys = np.arange(len(order))[::-1]
    cols = [model_col(t, m) for m in order]

    ax = axes[0]
    ranges(ax, ys, [g.loc[g["model"] == m, "dHeldout"] + 1 for m in order], cols)
    ax.set_yticks(ys, order, fontsize=9.5)
    for y, m in zip(ys, order):
        w = (g.loc[g["model"] == m, "dHeldout"] == 0).sum()
        ax.text(1.01, y, f"best in {w:,} ({100 * w / n:.0f}%)", transform=ax.get_yaxis_transform(),
                va="center", fontsize=8.5, color=TEXT)
    ax.set_xlabel("held-out score − sample's best, + 1  (log scale; 1 = this model is the best)", fontsize=8.5, color=MUTED)
    title(ax, f"{t}: how far behind the best model, per sample ({n:,} samples)",
          "dot = median, thick bar = middle 50% of samples, line = middle 80%")
    tidy(ax)

    ax = axes[1]
    w = np.array([(g.loc[g["model"] == m, "dHeldout"] == 0).mean() * 100 for m in order])
    close = np.array([g.loc[g["model"] == m, "close to best"].mean() * 100 for m in order])
    ax.barh(ys, close, height=0.55, color=cols, alpha=0.35)
    ax.barh(ys, w, height=0.55, color=cols)
    for y, a, b in zip(ys, w, close):
        ax.text(b + 1.5, y, f"{a:.0f}% best · {b:.0f}% within 10", va="center", fontsize=8, color=TEXT)
    ax.set_yticks(ys, [])
    ax.set_xlim(0, 135)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("% of samples", fontsize=8.5, color=MUTED)
    title(ax, "How often it wins", "solid = best · pale = within 10 of the best (held-out)")
    tidy(ax)

    ax = axes[2]
    off = 0.17
    for y, m in zip(ys, order):
        c = g[g["model"] == m]
        for dy, col, mk in ((off, DROP, "v"), (-off, SPIKE, "^")):
            q1, med, q3 = np.percentile(c[col], [25, 50, 75])
            ax.plot([q1, q3], [y + dy, y + dy], color=MUTED, lw=1.5)
            ax.plot(med, y + dy, mk, ms=8, color="#eb6834" if col == DROP else "#1baf7a")
    ax.axvline(0.1, color=TEXT, ls="--", lw=1)
    ax.text(0.1, ys[0] + 0.55, " what a fitting model gives", fontsize=7.5, color=TEXT)
    ax.set_yticks(ys, [])
    ax.set_ylim(ys[-1] - 0.6, ys[0] + 0.8)
    ax.set_xlabel("% of windows beyond the model's 0.1% tail\n(median, middle 50% of samples)", fontsize=8, color=MUTED)
    ax.plot([], [], "v", color="#eb6834", ms=8, label="biggest drops")
    ax.plot([], [], "^", color="#1baf7a", ms=8, label="biggest spikes")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    title(ax, "Extreme kicks as common as the model says?")
    tidy(ax)

    fig.suptitle(f"{t}: the four models on every QC-pass Pf8 {t} sample (chosen: {chosen_name(t)}, blue)",
                 x=0.01, ha="left", fontsize=12.5, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def plot_criteria(res, path):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), gridspec_kw={"width_ratios": [1, 1, 0.75]})
    crit = {"AIC": "dAIC", "BIC": "dBIC", "Held-out": "dHeldout"}
    agree = {}
    for ax, t in zip(axes, ["sWGA", "gDNA"]):
        g = res[res["sample type"] == t]
        order = model_order(t)
        ys = np.arange(len(order))[::-1]
        picks = {c: g.loc[g[d] == 0].drop_duplicates("sample").set_index("sample")["model"] for c, d in crit.items()}
        h = 0.26
        for i, (c, s) in enumerate(picks.items()):
            share = s.value_counts(normalize=True).reindex(order).fillna(0) * 100
            bars = ax.barh(ys + (1 - i) * h, share, height=h - 0.03, color=CRIT_COL[c], label=c)
            for b, v in zip(bars, share):
                if v >= 0.5:
                    ax.text(v + 1, b.get_y() + b.get_height() / 2, f"{v:.0f}%", va="center", fontsize=7.5, color=TEXT)
        ax.set_yticks(ys, order, fontsize=9.5)
        ax.set_xlim(0, 105)
        ax.set_xlabel("% of samples where this model is the best", fontsize=8.5, color=MUTED)
        title(ax, f"{t} ({g['sample'].nunique():,} samples): the model each criterion picks")
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        tidy(ax)
        agree[t] = {f"{a} = {b}": (picks[a] == picks[b]).mean() * 100
                    for a, b in (("BIC", "Held-out"), ("AIC", "Held-out"), ("AIC", "BIC"))}
    ax = axes[2]
    pairs = list(agree["sWGA"])
    ys = np.arange(len(pairs))[::-1]
    for i, t in enumerate(["sWGA", "gDNA"]):
        v = [agree[t][p] for p in pairs]
        bars = ax.barh(ys + (0.5 - i) * 0.36, v, height=0.33, color=TYPE_COL[t], label=t)
        for b, x in zip(bars, v):
            ax.text(x + 1, b.get_y() + b.get_height() / 2, f"{x:.0f}%", va="center", fontsize=8, color=TEXT)
    ax.set_yticks(ys, [p.replace(" = ", " and ") for p in pairs], fontsize=9.5)
    ax.set_xlim(0, 112)
    ax.set_xlabel("% of samples where both pick the same model", fontsize=8.5, color=MUTED)
    title(ax, "How often two criteria agree")
    ax.legend(frameon=False, fontsize=8.5, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.14))
    tidy(ax)
    fig.suptitle("Model choice by AIC, BIC and held-out likelihood (held-out decides, BIC checks)",
                 x=0.01, ha="left", fontsize=12.5, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def plot_all(res, shared, info, path):
    ch = res[res["chosen"]].set_index("sample")
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    types = ["sWGA", "gDNA"]

    ax = axes[0, 0]
    bins = np.geomspace(2, 400, 70)
    for t in types:
        v = ch.loc[ch["sample type"] == t, "nu"]
        ax.hist(v, bins=bins, histtype="step", lw=2, color=TYPE_COL[t],
                label=f"{t} ({chosen_name(t)}): median {v.median():.1f}")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(comma)
    ax.set_xlabel("ν of the t kicks (log scale; small = frequent big jumps, > 30 ≈ Gaussian)", fontsize=8.5, color=MUTED)
    ax.set_ylabel("samples", fontsize=8.5, color=MUTED)
    title(ax, "Tail weight of the chosen model in every sample")
    ax.legend(frameon=False, fontsize=8.5)
    tidy(ax, "y")

    ax = axes[0, 1]
    rows = [(t, src, col) for t in types for src in ("own fit", "shared") for col in (DROP, SPIKE)]
    ys = np.arange(len(rows))[::-1]
    data, cols, labels = [], [], []
    for t, src, col in rows:
        d = (ch[ch["sample type"] == t] if src == "own fit" else shared[shared["sample type"] == t])[col]
        data.append(d.to_numpy() + 0.01)
        cols.append(TYPE_COL[t])
        labels.append(f"{t}, {src}: {'drops' if col == DROP else 'spikes'}")
    ranges(ax, ys, data, cols)
    ax.axvline(0.11, color=TEXT, ls="--", lw=1)
    ax.set_yticks(ys, labels, fontsize=8.5)
    ax.set_xlabel("% of windows beyond the model's 0.1% tail, + 0.01 (log scale; dashed = what a fitting model gives)",
                  fontsize=8, color=MUTED)
    title(ax, "Extreme kicks: as often as the model says?",
          "each sample's own fit, and the one shared parameter set per type; dot = median, bars = middle 50% / 80%")
    tidy(ax)

    ax = axes[1, 0]
    gl = []
    for t in types:
        (p, q), _ = CHOSEN[t]
        g = res[(res["sample type"] == t) & (res["correlation"] == label(p, q))].pivot(index="sample", columns="kicks")
        gl.append(pd.DataFrame({"sample type": t, "loss": (g["LL", "t"] - g["LL", "Gaussian"]) / g["n", "t"]}))
    gl = pd.concat(gl)
    bins = np.linspace(-0.05, 0.6, 80)
    for t in types:
        ax.hist(shared.loc[shared["sample type"] == t, "loss per window"].clip(upper=0.6), bins=bins, histtype="step",
                lw=2, color=TYPE_COL[t], label=f"{t}: shared t parameters instead of its own t fit")
        ax.hist(gl.loc[gl["sample type"] == t, "loss"].clip(upper=0.6), bins=bins, histtype="step", lw=1.5, ls=":",
                color=TYPE_COL[t], label=f"{t}: its own Gaussian fit instead of its own t fit")
    ax.set_xlabel("log-likelihood lost per 500 bp window (0 = as good as the sample's own t fit; > 0.6 piled at 0.6)",
                  fontsize=8.5, color=MUTED)
    ax.set_ylabel("samples", fontsize=8.5, color=MUTED)
    med = {t: (shared.loc[shared["sample type"] == t, "loss per window"].median(),
               gl.loc[gl["sample type"] == t, "loss"].median()) for t in types}
    title(ax, "One shared parameter set per type (what krsim would ship) vs each sample's own fit",
          "median loss, shared t / own Gaussian: " + ", ".join(f"{t} {a:.3f} / {b:.3f}" for t, (a, b) in med.items()))
    ax.legend(frameon=False, fontsize=8)
    tidy(ax, "y")

    ax = axes[1, 1]
    ch = ch.join(info[["Population"]])
    by = (ch.groupby(["Population", "sample type"])["close to best"].agg(["mean", "size"]).reset_index())
    pops = sorted(by["Population"].unique())
    ys = np.arange(len(pops))[::-1]
    for i, t in enumerate(types):
        b = by[by["sample type"] == t].set_index("Population").reindex(pops)
        bars = ax.barh(ys + (0.5 - i) * 0.38, b["mean"] * 100, height=0.35, color=TYPE_COL[t], label=t)
        for bar, (_, r) in zip(bars, b.iterrows()):
            if np.isfinite(r["mean"]):
                ax.text(r["mean"] * 100 + 1, bar.get_y() + bar.get_height() / 2, f"{r['mean'] * 100:.0f}% (n={r['size']:,.0f})",
                        va="center", fontsize=7, color=TEXT)
    ax.set_yticks(ys, pops, fontsize=8.5)
    ax.set_xlim(0, 125)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("% of samples where the chosen model is best or within 10 of the best (held-out)", fontsize=8.5, color=MUTED)
    title(ax, "By Pf8 population")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right", ncol=2, bbox_to_anchor=(1, 1.0))
    tidy(ax)

    n = {t: (ch["sample type"] == t).sum() for t in types}
    fig.suptitle(f"The chosen models on all QC-pass Pf8 samples: {n['sWGA']:,} sWGA ({chosen_name('sWGA')}), "
                 f"{n['gDNA']:,} gDNA ({chosen_name('gDNA')})", x=0.01, ha="left", fontsize=12.5, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    res = pd.read_csv(f"{OUT}/fits.tsv", sep="\t")
    shared = pd.read_csv(f"{OUT}/shared_scores.tsv", sep="\t")
    info = pd.read_csv(f"{OUT}/samples.tsv", sep="\t", index_col=0)
    os.makedirs(f"{OUT}/figures", exist_ok=True)
    plot_type(res, "sWGA", f"{OUT}/figures/1_models_sWGA.png")
    plot_type(res, "gDNA", f"{OUT}/figures/2_models_gDNA.png")
    plot_criteria(res, f"{OUT}/figures/3_selection_criteria.png")
    plot_all(res, shared, info, f"{OUT}/figures/4_all_samples.png")
    print(f"wrote 4 figures in {OUT}/figures")


if __name__ == "__main__":
    main()
