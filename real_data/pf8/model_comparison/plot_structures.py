"""
Figures for fit_structures.py: which correlation structure is best, on every QC-pass Pf8 gDNA and sWGA sample.

    python plot_structures.py          # reads OUT/structures/structures.tsv, writes OUT/figures/

0a_structures.png           per sample type: every structure (t kicks) by how far behind the best it is; how often
                            each one wins with Gaussian and with t kicks; how far back an AR model needs to look
0b_structure_criteria.png   which structure AIC, BIC and held-out likelihood pick (t kicks), and how often they agree
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from best_simple_model import label
from fit_chosen_models import CHOSEN, OUT
from fit_structures import STRUCT_OUT, STRUCTURES
from plot_chosen_models import CHOSEN_COL, CRIT_COL, MUTED, OTHER_T, TEXT, TYPE_COL, comma, ranges, tidy, title

KCOL = {"Gaussian": "#8a8984", "t": "#2a78d6"}
AR_P = [p for p, q in STRUCTURES if q == 0]


def chosen_corr(t):
    return label(*CHOSEN[t][0])


def by_median(g):
    """Structures from best to worst by median held-out score."""
    return g.groupby("correlation")["dHeldout"].median().sort_values().index.tolist()


def structures(res, path):
    types = ["sWGA", "gDNA"]
    fig, axes = plt.subplots(2, 3, figsize=(20, 10.5), gridspec_kw={"width_ratios": [1.35, 1, 0.95]})
    for row, t in enumerate(types):
        g = res[res["sample type"] == t]
        n = g["sample"].nunique()
        gt = g[g["kicks"] == "t"]
        order = by_median(gt)
        ys = np.arange(len(order))[::-1]
        cols = [CHOSEN_COL if c == chosen_corr(t) else OTHER_T for c in order]

        ax = axes[row, 0]
        ranges(ax, ys, [gt.loc[gt["correlation"] == c, "dHeldout"] + 1 for c in order], cols)
        for y, c in zip(ys, order):
            w = (gt.loc[gt["correlation"] == c, "dHeldout"] == 0).sum()
            ax.text(1.01, y, f"best in {w:,} ({100 * w / n:.1f}%)", transform=ax.get_yaxis_transform(), va="center",
                    fontsize=8.5, color=TEXT)
        k = gt.groupby("correlation")["k"].first()
        ax.set_yticks(ys, [f"{c}  (k={k[c]})" for c in order], fontsize=9.5)
        ax.set_xlabel("held-out score − best structure, + 1  (t kicks; log scale; 1 = best)", fontsize=8.5, color=MUTED)
        title(ax, f"{t} ({n:,} samples): every structure, t kicks, best to worst",
              f"chosen {chosen_corr(t)} in blue; dot = median, bars = middle 50% / 80% of samples")
        tidy(ax)

        ax = axes[row, 1]
        h = 0.38
        for i, kick in enumerate(["Gaussian", "t"]):
            gk = g[g["kicks"] == kick]
            share = (gk[gk["dHeldout"] == 0].groupby("correlation")["sample"].nunique() / n * 100).reindex(order).fillna(0)
            bars = ax.barh(ys + (0.5 - i) * h, share, height=h - 0.04, color=KCOL[kick], label=f"{kick} kicks")
            for b, v in zip(bars, share):
                if v >= 0.5:
                    ax.text(v + 1, b.get_y() + b.get_height() / 2, f"{v:.0f}%", va="center", fontsize=7.5, color=TEXT)
        ax.set_yticks(ys, order, fontsize=9.5)
        ax.set_xlim(0, 105)
        ax.set_xlabel("% of samples where this structure is the best (held-out)", fontsize=8.5, color=MUTED)
        title(ax, "Does the winner depend on the kicks?", "best structure with Gaussian kicks and with t kicks")
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        tidy(ax)

        ax = axes[row, 2]
        for kick, ls in (("Gaussian", ":"), ("t", "-")):
            gk = g[g["kicks"] == kick]
            med = [gk.loc[gk["correlation"] == f"AR({p})", "dHeldout"].median() + 1 for p in AR_P]
            ax.plot(AR_P, med, ls, marker="o", ms=7, lw=2, color=KCOL[kick], label=f"AR(p), {kick} kicks")
            arma = gk.loc[gk["correlation"] == "ARMA(2,2)", "dHeldout"].median() + 1
            ax.axhline(arma, color=KCOL[kick], ls=ls, lw=1)
            ax.text(1.05, arma, f"ARMA(2,2), {kick}", fontsize=7.5, color=TEXT, va="bottom" if kick == "t" else "top")
        ax.set_xscale("log", base=2)
        ax.set_xticks(AR_P, [f"{p}\n({p * 0.5:g} kb)" for p in AR_P])
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(comma)
        ax.set_xlim(0.8, 40)
        ax.set_xlabel("AR order p (how far back it looks)", fontsize=8.5, color=MUTED)
        ax.set_ylabel("median held-out score − best, + 1 (log)", fontsize=8.5, color=MUTED)
        title(ax, "How far back does depth need to look?", "AR(p) against p; lines = ARMA(2,2)")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        tidy(ax, "y")
    fig.suptitle("Correlation structure of log depth along the genome, every QC-pass Pf8 sample "
                 "(all on the same windows; held-out = fit odd chromosomes, score even, then swap)",
                 x=0.01, ha="left", fontsize=13, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def criteria(res, path):
    crit = {"AIC": "dAIC", "BIC": "dBIC", "Held-out": "dHeldout"}
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6), gridspec_kw={"width_ratios": [1, 1, 0.75]})
    agree = {}
    for ax, t in zip(axes, ["sWGA", "gDNA"]):
        g = res[(res["sample type"] == t) & (res["kicks"] == "t")]
        order = by_median(g)
        ys = np.arange(len(order))[::-1]
        picks = {c: g.loc[g[d] == 0].drop_duplicates("sample").set_index("sample")["correlation"] for c, d in crit.items()}
        h = 0.27
        for i, (c, s) in enumerate(picks.items()):
            share = s.value_counts(normalize=True).reindex(order).fillna(0) * 100
            bars = ax.barh(ys + (1 - i) * h, share, height=h - 0.03, color=CRIT_COL[c], label=c)
            for b, v in zip(bars, share):
                if v >= 0.5:
                    ax.text(v + 1, b.get_y() + b.get_height() / 2, f"{v:.0f}%", va="center", fontsize=7, color=TEXT)
        k = g.groupby("correlation")["k"].first()
        ax.set_yticks(ys, [f"{c}  (k={k[c]})" for c in order], fontsize=9)
        ax.set_xlim(0, 105)
        ax.set_xlabel("% of samples where this structure is the best (t kicks)", fontsize=8.5, color=MUTED)
        title(ax, f"{t} ({g['sample'].nunique():,} samples): the structure each criterion picks")
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        tidy(ax)
        agree[t] = {f"{a} and {b}": np.equal(*picks[a].align(picks[b])).mean() * 100
                    for a, b in (("BIC", "Held-out"), ("AIC", "Held-out"), ("AIC", "BIC"))}
    ax = axes[2]
    pairs = list(agree["sWGA"])
    ys = np.arange(len(pairs))[::-1]
    for i, t in enumerate(["sWGA", "gDNA"]):
        v = [agree[t][p] for p in pairs]
        bars = ax.barh(ys + (0.5 - i) * 0.36, v, height=0.33, color=TYPE_COL[t], label=t)
        for b, x in zip(bars, v):
            ax.text(x + 1, b.get_y() + b.get_height() / 2, f"{x:.0f}%", va="center", fontsize=8, color=TEXT)
    ax.set_yticks(ys, pairs, fontsize=9.5)
    ax.set_xlim(0, 112)
    ax.set_xlabel("% of samples where both pick the same structure", fontsize=8.5, color=MUTED)
    title(ax, "How often two criteria agree")
    ax.legend(frameon=False, fontsize=8.5, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.14))
    tidy(ax)
    fig.suptitle("Structure choice by AIC, BIC and held-out likelihood, t kicks (k = number of parameters)",
                 x=0.01, ha="left", fontsize=13, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    res = pd.read_csv(f"{STRUCT_OUT}/structures.tsv", sep="\t")
    os.makedirs(f"{OUT}/figures", exist_ok=True)
    structures(res, f"{OUT}/figures/0a_structures.png")
    criteria(res, f"{OUT}/figures/0b_structure_criteria.png")
    print(f"wrote 0a_structures.png and 0b_structure_criteria.png in {OUT}/figures")


if __name__ == "__main__":
    main()
