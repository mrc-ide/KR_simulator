"""
How much Student t kicks improved on Gaussian kicks, from the tables written by best_simple_model.py.

    python plot_kick_gain.py            # reads ../simple_models/, writes the PNGs here

1_kick_shapes.png   what changes: a Gaussian and a t with the same SD (nu = 7 as gDNA, 3.6 as sWGA)
2_score_gain.png    held-out score gained by t kicks, per sample: every structure, and kicks vs structure
3_tail_check.png    chosen model per type: how often kicks land in the model's own tails, Gaussian vs t
4_all_models.png    every structure x Gaussian / t kicks: held-out score behind each sample's best model
4b_all_models_new20.png  the same for the 20 new samples (best_simple_model.py ~/pf8_data/pf8_samples ../simple_models_new20)
5_kicks_compared.png  simulated Gaussian vs t kicks (same SD) and the depth they make along a genome
6_how_scored.png    how the score is built: a penalty per window, added up; simulated sWGA-like kicks
7a/7b_real_*.png    5_kicks_compared on a real gDNA and sWGA sample (needs the depth data in PF8_SAMPLES)
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(HERE, "..", "simple_models")
PF8_SAMPLES = os.path.expanduser("~/pf8_data/pf8_samples")      # depth data for figures 7a/7b (change on another computer)
TEXT, MUTED, GRID = "#0b0b0b", "#6b6a66", "#e6e5e1"
GAUSS, T = "#a3a29d", "#2a78d6"
CHOSEN = {"gDNA": "ARMA(2,2)", "sWGA": "AR(8)"}
NU = {"gDNA": 7.0, "sWGA": 3.6}
CORR = ["AR(1)", "AR(2)", "AR(4)", "AR(8)", "ARMA(1,1)", "ARMA(2,1)", "ARMA(2,2)"]
TYPES = ["gDNA", "sWGA"]


def tidy(ax, grid="y"):
    ax.grid(axis=grid, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def kick_shapes(out_png):
    """Gaussian vs t with the same SD (t scale = sqrt((nu - 2) / nu)), so only the shape differs."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    z = np.linspace(-8, 8, 801)
    ax = axes[0]
    ax.plot(z, stats.norm.pdf(z), color=GAUSS, lw=2, label="Gaussian")
    for t, ls in zip(TYPES, ["--", "-"]):
        s = np.sqrt((NU[t] - 2) / NU[t])
        ax.plot(z, stats.t.pdf(z / s, NU[t]) / s, color=T, lw=2, ls=ls, label=f"t, ν = {NU[t]} ({t})")
    ax.set_yscale("log")
    ax.set_ylim(1e-9, 1)
    ax.set_xlabel("kick, in SDs", fontsize=9)
    ax.set_ylabel("probability density (log scale)", fontsize=9)
    ax.set_title("Same spread, different tails: the t keeps big jumps possible", loc="left", fontsize=10, color=TEXT)
    ax.legend(frameon=False, fontsize=8.5)
    tidy(ax)

    ax = axes[1]
    k = np.array([2, 3, 4, 5, 6])
    x = np.arange(len(k))
    tail = {"Gaussian": 2 * stats.norm.sf(k)}
    for t in TYPES:
        s = np.sqrt((NU[t] - 2) / NU[t])
        tail[f"t, ν = {NU[t]} ({t})"] = 2 * stats.t.sf(k / s, NU[t])
    cols = [GAUSS, "#8db8ea", T]
    for i, (name, v) in enumerate(tail.items()):
        bars = ax.bar(x + (i - 1) * 0.27, v, width=0.25, color=cols[i], label=name)
        for b, p in zip(bars, v):
            ax.text(b.get_x() + b.get_width() / 2, p * 1.4, f"1 in {1 / p:,.0f}", ha="center", va="bottom",
                    fontsize=6.5, color=TEXT, rotation=90)
    ax.set_yscale("log")
    ax.set_ylim(1e-10, 30)
    ax.set_xticks(x, [f"beyond {v} SD" for v in k], fontsize=9)
    ax.set_ylabel("share of kicks (log scale)", fontsize=9)
    ax.set_title("How often a kick is this big", loc="left", fontsize=10, color=TEXT)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    tidy(ax)
    fig.suptitle("What changes when the kicks are Student t instead of Gaussian", x=0.01, ha="left",
                 fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def score_gain(res, out_png):
    """Left: t gain over Gaussian for every structure and sample. Right: kicks gain vs structure gain."""
    s = res.pivot_table(index=["sample type", "sample", "correlation"], columns="kicks", values="Heldout").reset_index()
    s["gain"] = s["Gaussian"] - s["t"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), gridspec_kw={"width_ratios": [1.7, 1]})
    rng = np.random.default_rng(1)
    fmt = FuncFormatter(lambda v, _: f"{v:,g}")
    for row, t in enumerate(TYPES):
        g = s[s["sample type"] == t]
        n = g["sample"].nunique()
        ax = axes[row, 0]
        for i, c in enumerate(CORR):
            v = g.loc[g["correlation"] == c, "gain"].to_numpy()
            x = i + rng.uniform(-0.15, 0.15, len(v))
            ax.scatter(x, np.maximum(v, 3), s=22, color=T, alpha=0.75, lw=0)
            for xi, vi in zip(x[v < 3], v[v < 3]):        # t no better: drawn on the floor, with its value
                ax.scatter(xi, 3, s=40, facecolor="white", edgecolor=TEXT, zorder=3)
                ax.annotate(f"{vi:,.0f} (a tie)", (xi, 3), (xi + 0.1, 5), fontsize=7.5, color=TEXT)
            m = np.median(v)
            ax.plot([i - 0.3, i + 0.3], [m, m], color=TEXT, lw=1.5)
            ax.text(i + 0.33, m, f"{m:,.0f}", va="center", fontsize=7.5, color=TEXT)
        ax.set_yscale("log")
        ax.set_ylim(2, None)
        ax.yaxis.set_major_formatter(fmt)
        ax.set_xticks(range(len(CORR)), CORR, fontsize=9)
        ax.set_ylabel("held-out score: Gaussian − t\n(log scale; > 0 = t better)", fontsize=9)
        ax.set_title(f"{t} ({n} samples): t better in {(g['gain'] > 0).sum()} of {len(g)} fits; "
                     f"bar = median", loc="left", fontsize=10, color=TEXT)
        tidy(ax)

        # paired per sample: switching kicks (chosen structure) vs switching structure (AR(1) -> chosen, t kicks)
        ax = axes[row, 1]
        w = g.set_index(["sample", "correlation"])
        kicks = w.xs(CHOSEN[t], level="correlation")["gain"]
        struct = w.xs("AR(1)", level="correlation")["t"] - w.xs(CHOSEN[t], level="correlation")["t"]
        for smp in kicks.index:
            ax.plot([0, 1], [struct[smp], kicks[smp]], color=GRID, lw=1, zorder=1)
        ax.scatter(np.zeros(len(struct)), struct, s=30, color=MUTED, zorder=2)
        ax.scatter(np.ones(len(kicks)), kicks, s=30, color=T, zorder=2)
        for x, v in ((0, struct), (1, kicks)):
            ax.text(x + 0.08, v.median(), f"median {v.median():,.0f}", va="center", fontsize=8, color=TEXT)
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(fmt)
        ax.set_xlim(-0.3, 1.6)
        ax.set_xticks([0, 1], [f"structure\nAR(1) → {CHOSEN[t]}, t kicks", f"kicks\nGaussian → t, {CHOSEN[t]}"],
                      fontsize=8.5)
        ax.set_ylabel("held-out score gained (log scale)", fontsize=9)
        ax.set_title("What mattered more? (one line per sample)", loc="left", fontsize=10, color=TEXT)
        tidy(ax)
    fig.suptitle("Held-out score gained by switching to t kicks (lower score = better fit to unseen chromosomes)",
                 x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def tail_check(calib, out_png):
    """Chosen model per type: share of kicks in each of the model's tails, Gaussian vs t, one line per sample."""
    cols = [("lowest 0.1% (expect 0.1)", "biggest drops\n(model's lowest 0.1%)", 0.1),
            ("highest 0.1% (expect 0.1)", "biggest spikes\n(model's highest 0.1%)", 0.1),
            ("lowest 1% (expect 1)", "big drops\n(model's lowest 1%)", 1.0),
            ("highest 1% (expect 1)", "big spikes\n(model's highest 1%)", 1.0),
            ("KS distance", "overall mismatch\n(KS distance)", 0.0)]
    fig, axes = plt.subplots(2, len(cols), figsize=(15, 8.5))
    for row, t in enumerate(TYPES):
        g = calib[(calib["sample type"] == t) & calib["model"].str.startswith(CHOSEN[t] + " ")]
        g = g.assign(kicks=g["model"].str.split().str[-1]).pivot(index="sample", columns="kicks")
        for j, (col, title, expect) in enumerate(cols):
            ax = axes[row, j]
            a, b = g[(col, "Gaussian")], g[(col, "t")]
            for smp in g.index:
                ax.plot([0, 1], [a[smp], b[smp]], color=GRID, lw=1, zorder=1)
            ax.scatter(np.zeros(len(a)), a, s=26, color=GAUSS, zorder=2)
            ax.scatter(np.ones(len(b)), b, s=26, color=T, zorder=2)
            ax.axhline(expect, color=TEXT, ls="--", lw=1)
            for x, v in ((0, a), (1, b)):
                ax.text(x + 0.1, v.median(), f"{v.median():.2f}", va="center", fontsize=8, color=TEXT,
                        bbox=dict(facecolor="white", edgecolor="none", pad=1))
            ax.set_xlim(-0.35, 1.5)
            ax.set_ylim(0, max(a.max(), b.max(), expect) * 1.12)
            ax.set_xticks([0, 1], ["Gaussian", "t"], fontsize=9)
            ax.set_title(title, loc="left", fontsize=9.5, color=TEXT)
            if j == 0:
                ax.set_ylabel(f"{t}: {CHOSEN[t]}\n% of kicks", fontsize=10)
            elif j == len(cols) - 1:
                ax.set_ylabel("0 = perfect", fontsize=8.5)
            tidy(ax)
    fig.suptitle("Do extreme kicks happen as often as the model says? Dashed line = a perfect fit; "
                 "one line per sample", x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def all_models(res, out_png, which="the 21 comparison samples"):
    """Every structure x Gaussian / t: held-out score behind the sample's best model, one dot per sample."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    rng = np.random.default_rng(1)
    for ax, t in zip(axes, TYPES):
        g = res[res["sample type"] == t]
        n = g["sample"].nunique()
        for i, c in enumerate(CORR):
            for j, (kick, col) in enumerate([("Gaussian", GAUSS), ("t", T)]):
                v = g.loc[(g["correlation"] == c) & (g["kicks"] == kick), "dHeldout"].to_numpy()
                x = i + (j - 0.5) * 0.36
                ax.scatter(x + rng.uniform(-0.1, 0.1, len(v)), v, s=20, color=col, alpha=0.8, lw=0, zorder=2,
                           label=f"{kick} kicks" if i == 0 else None)
                m = np.median(v)
                ax.plot([x - 0.15, x + 0.15], [m, m], color=TEXT, lw=1.5, zorder=3)
                ax.text(x + (0.17 if j else -0.17), m, f"{m:,.0f}", ha="left" if j else "right", va="center",
                        fontsize=7.5, color=TEXT, bbox=dict(facecolor="white", edgecolor="none", pad=0.5), zorder=4)
                wins = int((v == 0).sum())
                if wins:
                    ax.text(x, -0.02, f"best in {wins} of {n}", ha="center", va="top", fontsize=8, color=TEXT,
                            fontweight="bold", transform=ax.get_xaxis_transform())
        ax.set_yscale("symlog", linthresh=10, linscale=0.5)
        top = 10 ** np.ceil(np.log10(g["dHeldout"].max() * 1.5))
        ax.set_ylim(-1.6, top)
        ax.set_yticks([v for v in (0, 10, 100, 1000, 10000, 100000) if v <= top])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.set_xticks(range(len(CORR)), CORR, fontsize=10)
        ax.tick_params(axis="x", pad=18)
        ax.set_ylabel("held-out score behind the sample's best model\n(0 = best; log scale above 10; lower = better)",
                      fontsize=9)
        ax.set_title(f"{t} ({n} samples): chosen = {CHOSEN[t]} + t kicks", loc="left", fontsize=11, color=TEXT)
        ax.legend(frameon=False, fontsize=9, loc="upper right", ncol=2)
        tidy(ax)
    fig.suptitle(f"Every model shape with Gaussian or t kicks, {which}: dots = samples, bar and number = median",
                 x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def t_kicks(rng, nu, n):
    """t kicks rescaled to SD 1, so they differ from Gaussian kicks only in shape."""
    return rng.standard_t(nu, n) * np.sqrt((nu - 2) / nu)


def kicks_compared(out_png):
    """Simulated kicks and the depth they make: Gaussian vs t (nu = 7, 3.6), same SD. Illustration, not a fit."""
    n, phi, sd = 1500, 0.8, 0.15           # an AR(1) with roughly realistic spread of log depth; illustrative only
    rng = np.random.default_rng(7)
    kinds = [("Gaussian kicks", rng.standard_normal(n), GAUSS),
             ("t kicks, ν = 7 (like gDNA)", t_kicks(rng, 7.0, n), "#8db8ea"),
             ("t kicks, ν = 3.6 (like sWGA)", t_kicks(rng, 3.6, n), T)]
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), sharey="row")
    for j, (name, e, col) in enumerate(kinds):
        ax = axes[0, j]
        big = np.abs(e) > 3
        ec = np.clip(e, -8, 8)              # one or two huge kicks would squash the rest; drawn at the edge
        ax.vlines(np.arange(n), 0, ec, color=col, lw=0.6)
        ax.scatter(np.flatnonzero(big), ec[big], s=14, color=TEXT, zorder=3)
        for i in np.flatnonzero(np.abs(e) > 8):
            ax.annotate(f"{e[i]:.0f} SD", (i, ec[i]), (i + 40, ec[i] * 0.95), fontsize=8, color=TEXT)
        for h in (-3, 3):
            ax.axhline(h, color=TEXT, ls="--", lw=0.8)
        ax.set_ylim(-8.8, 8.8)
        ax.set_title(f"{name}\n{big.sum()} of {n:,} kicks beyond 3 SD (black dots)", loc="left", fontsize=10,
                     color=TEXT)
        if j == 0:
            ax.set_ylabel("kick at each window\n(in SDs)", fontsize=9)
        tidy(ax)

        y = np.zeros(n)                     # log depth: new = 0.8 x previous + kick
        for i in range(1, n):
            y[i] = phi * y[i - 1] + sd * np.sqrt(1 - phi ** 2) * e[i]
        ax = axes[1, j]
        ax.plot(np.arange(n), 50 * np.exp(y), color=col, lw=0.8)
        ax.set_xlabel("window along the genome", fontsize=9)
        if j == 0:
            ax.set_ylabel("depth made by these kicks\n(mean 50 reads)", fontsize=9)
        tidy(ax)
    fig.suptitle("The same model driven by Gaussian or t kicks (all with SD 1): t gives mostly calm windows "
                 "plus sudden dips and spikes", x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def how_scored(out_png):
    """The score: each window's kick gets a penalty, -2 log(probability the model gave it); add them all up."""
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), gridspec_kw={"width_ratios": [1, 1.3]})
    z = np.linspace(-8, 8, 801)
    s = np.sqrt((3.6 - 2) / 3.6)
    pen = {"Gaussian": -2 * stats.norm.logpdf(z), "t, ν = 3.6": -2 * (stats.t.logpdf(z / s, 3.6) - np.log(s))}
    ax = axes[0]
    for (name, p), col in zip(pen.items(), [GAUSS, T]):
        ax.plot(z, p, color=col, lw=2, label=name)
    for k in (1, 6):
        g = -2 * stats.norm.logpdf(k)
        t = -2 * (stats.t.logpdf(k / s, 3.6) - np.log(s))
        ax.scatter([k, k], [g, t], color=[GAUSS, T], s=40, zorder=3, edgecolor=TEXT, lw=0.6)
        ax.annotate(f"kick of {k} SD:\nGaussian {g:.1f}, t {t:.1f}", (k, max(g, t)),
                    (k - 4.2, max(g, t) + 6) if k > 3 else (k + 1.5, max(g, t) + 12),
                    fontsize=8.5, color=TEXT, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    ax.set_ylim(0, 75)
    ax.set_xlabel("kick, in SDs", fontsize=9)
    ax.set_ylabel("penalty for this window\n= −2 × log(probability the model gave it)", fontsize=9)
    ax.set_title("Each window's penalty: small kicks cost about the same,\nbig kicks cost far more under a Gaussian",
                 loc="left", fontsize=10, color=TEXT)
    ax.legend(frameon=False, fontsize=9, loc="upper center")
    tidy(ax)

    # the same t-like kicks scored by both models; each model fitted to them (Gaussian: their SD; t: nu, scale)
    rng = np.random.default_rng(7)
    e = t_kicks(rng, 3.6, 3000)
    nu, loc, scale = stats.t.fit(e, floc=0)
    pg = -2 * stats.norm.logpdf(e, 0, e.std())
    pt = -2 * stats.t.logpdf(e, nu, loc, scale)
    ax = axes[1]
    gap = np.cumsum(pg - pt)
    ax.plot(gap, color=TEXT, lw=1.2)
    big = np.abs(e) > 4
    ax.scatter(np.flatnonzero(big), gap[big], s=14, color="#eb6834", zorder=3, label="windows with a kick beyond 4 SD")
    ax.set_xlabel("window along the genome", fontsize=9)
    ax.set_ylabel("running total: Gaussian score − t score\n(> 0 = t ahead)", fontsize=9)
    ax.set_title(f"Score = all the penalties added up. Same {len(e):,} windows scored by both models:\n"
                 f"t ends {gap[-1]:,.0f} better: big steps at big kicks, plus a steady climb on ordinary windows",
                 loc="left", fontsize=10, color=TEXT)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    tidy(ax)
    fig.suptitle("How the score is calculated (simulated sWGA-like kicks)", x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def real_compared(sample, sample_type, path, out_png, chrom=6):
    """Like 5_kicks_compared, on a real sample: its kicks and depth next to genomes simulated from its own
    Gaussian and t fits (chosen structure for its type), on one chromosome (index into the nuclear ones)."""
    import sys
    from scipy.signal import lfilter
    sys.path.insert(0, os.path.join(HERE, "..", ".."))
    from compare_arma_kicks import SKIP, arma_kicks, fit, unpack
    from compare_processes import PSEUDO, TRIM_KB, WINDOW, load, pacf_to_coefs

    p, q = {"gDNA": (2, 2), "sWGA": (8, 0)}[sample_type]
    ys = load(path)
    # undo load(): log depth -> reads, using the sample's median and the mean it removed
    df = pd.read_csv(path, sep="\t")
    df = df[~df["chrom"].str.contains("API|MIT")]
    med = df["mean"].median()
    trim = TRIM_KB * 1000 // WINDOW
    raw = [np.log(df.loc[df["chrom"] == c, "mean"].to_numpy() / med + PSEUDO)[trim:-trim] for c in df["chrom"].unique()]
    mu_y = np.concatenate(raw).mean()
    to_depth = lambda y: np.maximum(med * (np.exp(y + mu_y) - PSEUDO), 0)

    xg, _ = fit(ys, p, q, "Gaussian")
    xt, _ = fit(ys, p, q, "t", xg[:p + q])
    real = arma_kicks(ys, xt, p, q)
    sd = np.concatenate(real).std()
    n = len(real[chrom])
    rng = np.random.default_rng(3)

    def simulate(x, kick):
        d = unpack(kick, x[p + q:])
        burn = 500
        e = rng.standard_normal(n + burn) * d["sd"] if kick == "Gaussian" else \
            rng.standard_t(d["nu"], n + burn) * d["scale"]
        ar = np.r_[1, -pacf_to_coefs(np.tanh(x[:p]))]
        ma = np.r_[1, -pacf_to_coefs(np.tanh(x[p:p + q]))]
        return e[burn:], lfilter(ma, ar, e + d["mu"])[burn:]

    nu = unpack("t", xt[p + q:])["nu"]
    cols = [("Real data", real[chrom], ys[chrom][SKIP:], TEXT),
            ("Simulated from the Gaussian fit", *simulate(xg, "Gaussian"), GAUSS),
            (f"Simulated from the t fit (ν = {nu:.1f})", *simulate(xt, "t"), T)]
    lim = np.abs(real[chrom]).max() * 1.1
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), sharey="row")
    for j, (name, e, y, col) in enumerate(cols):
        ax = axes[0, j]
        big = np.abs(e) > 3 * sd
        ec = np.clip(e, -lim, lim)
        ax.vlines(np.arange(n), 0, ec, color=col, lw=0.5)
        ax.scatter(np.flatnonzero(big), ec[big], s=12, color="#eb6834", zorder=3)
        for h in (-3 * sd, 3 * sd):
            ax.axhline(h, color=TEXT, ls="--", lw=0.8)
        ax.set_ylim(-lim, lim)
        ax.set_title(f"{name}\n{big.sum()} of {n:,} kicks beyond 3 SD (orange dots)", loc="left", fontsize=10,
                     color=TEXT)
        if j == 0:
            ax.set_ylabel("kick at each window\n(change in log depth)", fontsize=9)
        tidy(ax)
        ax = axes[1, j]
        ax.plot(np.arange(n), to_depth(y) + 1, color=col, lw=0.6)
        ax.set_yscale("log")
        ax.set_yticks([1, 11, 101, 1001, 10001], ["0", "10", "100", "1,000", "10,000"])
        ax.minorticks_off()
        ax.set_xlabel(f"500 bp window along chromosome {chrom + 1}", fontsize=9)
        if j == 0:
            ax.set_ylabel("read depth (log scale)", fontsize=9)
        tidy(ax)
    allk = np.concatenate(real)
    exp_t = 2 * stats.t.sf(3 * sd / unpack("t", xt[p + q:])["scale"], nu)
    fig.suptitle(f"{sample} ({sample_type}, median depth {med:.0f}), model {CHOSEN[sample_type]}: real vs simulated. "
                 f"Whole genome: {100 * np.mean(np.abs(allk) > 3 * sd):.2f}% of real kicks beyond 3 SD; "
                 f"Gaussian expects 0.27%, the t fit {100 * exp_t:.2f}%",
                 x=0.01, ha="left", fontsize=11, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    res = pd.read_csv(os.path.join(TABLES, "simple_fits.tsv"), sep="\t")
    calib = pd.read_csv(os.path.join(TABLES, "simple_calibration.tsv"), sep="\t")
    kick_shapes(os.path.join(HERE, "1_kick_shapes.png"))
    score_gain(res, os.path.join(HERE, "2_score_gain.png"))
    tail_check(calib, os.path.join(HERE, "3_tail_check.png"))
    all_models(res, os.path.join(HERE, "4_all_models.png"))
    new20 = os.path.join(HERE, "..", "simple_models_new20", "simple_fits.tsv")
    if os.path.exists(new20):
        all_models(pd.read_csv(new20, sep="\t"), os.path.join(HERE, "4b_all_models_new20.png"),
                   "the 20 new samples")
    kicks_compared(os.path.join(HERE, "5_kicks_compared.png"))
    how_scored(os.path.join(HERE, "6_how_scored.png"))
    for s, t, tag in [("PD0048-C", "gDNA", "7a"), ("SPT43924", "sWGA", "7b")]:
        path = os.path.join(PF8_SAMPLES, s, "depth_plots", f"{s}.window_stats.tsv")
        if os.path.exists(path):
            real_compared(s, t, path, os.path.join(HERE, f"{tag}_real_{t}_{s}.png"))
    print(f"wrote figures in {HERE}")


if __name__ == "__main__":
    main()
