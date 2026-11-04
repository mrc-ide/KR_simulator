"""
Figures for fit_garch_kicks.py: Gaussian vs t vs GARCH t kicks on the chosen structure, every QC-pass sample.

    python plot_garch_kicks.py          # reads OUT/garch/kick_compare.tsv, writes OUT/figures/

5_kicks_evidence.png        for each sample type, over every sample: held-out score behind the best of the three,
                            what each step gains, how much big kicks bunch together, tail calibration
6a_kicks_explained_sWGA.png one real sample: each model's kick range along a real chromosome, depth simulated
6b_kicks_explained_gDNA.png from each model next to the real depth, bunching of big kicks, and the tails
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import lfilter

from best_simple_model import label
from compare_arma_kicks import SKIP, arma_kicks, garch_s2, unpack
from compare_processes import PSEUDO, TRIM_KB, WINDOW, pacf_to_coefs
from fit_chosen_models import CHOSEN, OUT, depth_path, load
from fit_garch_kicks import GARCH_OUT, KICKS, chosen_x, clustering, standardised
from plot_chosen_models import GRID, MUTED, TEXT, comma, ranges, tidy, title

KCOL = {"Gaussian": "#8a8984", "t": "#2a78d6", "GARCH t": "#eb6834"}
EXAMPLES = {"sWGA": "SPT43924", "gDNA": "PD0048-C"}      # the samples of results/kick_gain/ figures 7a/7b
DROP, SPIKE = "lowest 0.1% (expect 0.1)", "highest 0.1% (expect 0.1)"


def evidence(kc, path):
    types = ["sWGA", "gDNA"]
    fig, axes = plt.subplots(2, 4, figsize=(21, 9), gridspec_kw={"width_ratios": [1.25, 1.1, 0.9, 0.9]})
    order = ["GARCH t", "t", "Gaussian"]
    ys = np.arange(len(order))[::-1]
    for row, t in enumerate(types):
        g = kc[kc["sample type"] == t]
        n = g["sample"].nunique()
        (p, q), _ = CHOSEN[t]
        cols = [KCOL[k] for k in order]

        ax = axes[row, 0]
        ranges(ax, ys, [g.loc[g["kicks"] == k, "dHeldout"] + 1 for k in order], cols)
        for y, k in zip(ys, order):
            w = (g.loc[g["kicks"] == k, "dHeldout"] == 0).sum()
            ax.text(1.01, y, f"best in {w:,}\n({100 * w / n:.1f}%)", transform=ax.get_yaxis_transform(), va="center",
                    fontsize=8.5, color=TEXT)
        ax.set_yticks(ys, [f"{label(p, q)} {k}" for k in order], fontsize=9.5)
        ax.set_xlabel("held-out score − best of the three, + 1 (log scale)", fontsize=8.5, color=MUTED)
        title(ax, f"{t} ({n:,} samples): how far behind the best", "dot = median, bars = middle 50% / 80% of samples")
        tidy(ax)

        ax = axes[row, 1]
        wide = g.pivot(index="sample", columns="kicks", values="Heldout")
        steps = [("GARCH t over t", wide["t"] - wide["GARCH t"], KCOL["GARCH t"]),
                 ("t over Gaussian", wide["Gaussian"] - wide["t"], KCOL["t"])]
        sy = np.array([1, 0])
        ranges(ax, sy, [np.clip(v, 1, None) for _, v, _ in steps], [c for _, _, c in steps])
        for y, (name, v, _) in zip(sy, steps):
            ax.text(0.02, y + 0.3, f"median {v.median():,.0f} · better in {100 * (v > 0).mean():.1f}% of samples",
                    transform=ax.get_yaxis_transform(), fontsize=8.5, color=TEXT)
        ax.set_yticks(sy, [n for n, _, _ in steps], fontsize=9.5)
        ax.set_ylim(-0.6, 1.8)
        ax.set_xlabel("held-out gain per sample (−2 log-likelihood; log scale, < 1 shown at 1)", fontsize=8.5, color=MUTED)
        title(ax, "What each step gains", "on the same windows; ~44,000 windows per sample")
        tidy(ax)

        ax = axes[row, 2]
        ranges(ax, ys, [g.loc[g["kicks"] == k, "clustering"] for k in order], cols, log=False)
        ax.axvline(0, color=TEXT, ls="--", lw=1)
        ax.set_yticks(ys, order, fontsize=9.5)
        ax.set_xlabel("correlation of |kick| with the next |kick|\n(after dividing by the model's scale; 0 = no bunching)",
                      fontsize=8.5, color=MUTED)
        title(ax, "Do big kicks bunch together?", "what the model leaves unexplained")
        tidy(ax)

        ax = axes[row, 3]
        off = 0.17
        for y, k in zip(ys, order):
            c = g[g["kicks"] == k]
            for dy, col, mk, mc in ((off, DROP, "v", "#e34948"), (-off, SPIKE, "^", "#1baf7a")):
                q1, med, q3 = np.percentile(c[col], [25, 50, 75])
                ax.plot([q1, q3], [y + dy, y + dy], color=MUTED, lw=1.5)
                ax.plot(med, y + dy, mk, ms=8, color=mc)
        ax.axvline(0.1, color=TEXT, ls="--", lw=1)
        ax.set_yticks(ys, order, fontsize=9.5)
        ax.set_ylim(ys[-1] - 0.6, ys[0] + 0.6)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(comma)
        ax.plot([], [], "v", color="#e34948", ms=8, label="biggest drops")
        ax.plot([], [], "^", color="#1baf7a", ms=8, label="biggest spikes")
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        ax.set_xlabel("% of windows beyond the model's own 0.1% tail\n(median, middle 50%; dashed = a fitting model)",
                      fontsize=8.5, color=MUTED)
        title(ax, "Extreme kicks as often as expected?")
        tidy(ax)
    fig.suptitle("Gaussian vs t vs GARCH t kicks on every QC-pass Pf8 sample "
                 f"(sWGA: {label(*CHOSEN['sWGA'][0])}, gDNA: {label(*CHOSEN['gDNA'][0])})",
                 x=0.01, ha="left", fontsize=13, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def garch_x(sample):
    for f in sorted(os.listdir(f"{GARCH_OUT}/tasks")):
        if f.startswith("garch_"):
            df = pd.read_csv(f"{GARCH_OUT}/tasks/{f}", sep="\t", usecols=["sample", "x"])
            hit = df.loc[df["sample"] == sample, "x"]
            if len(hit):
                return np.array(hit.iloc[0].split(","), dtype=float)


def sim_kicks(rng, kick, d, n):
    """Kicks drawn from a fitted kick model (mu not added)."""
    if kick == "Gaussian":
        return rng.standard_normal(n) * d["sd"]
    z = rng.standard_t(d["nu"], n)
    if kick == "t":
        return z * d["scale"]
    e = np.empty(n)
    s2 = d["w"] / (1 - d["a"] - d["b"])
    for i in range(n):
        if i:
            s2 = d["w"] + d["a"] * e[i - 1] ** 2 + d["b"] * s2
        e[i] = np.sqrt(s2) * z[i]
    return e


def explained(sample, sample_type, xs, path, chrom=6):
    (p, q), _ = CHOSEN[sample_type]
    path_d = depth_path(sample)
    ys = load(path_d)
    df = pd.read_csv(path_d, sep="\t")
    df = df[~df["chrom"].str.contains("API|MIT")]
    med = df["mean"].median()
    trim = TRIM_KB * 1000 // WINDOW
    raw = [np.log(df.loc[df["chrom"] == c, "mean"].to_numpy() / med + PSEUDO)[trim:-trim] for c in df["chrom"].unique()]
    mu_y = np.concatenate(raw).mean()
    to_depth = lambda y: np.maximum(med * (np.exp(y + mu_y) - PSEUDO), 0)

    d = {k: unpack(k, xs[k][p + q:]) for k in KICKS}
    real = arma_kicks(ys, xs["t"], p, q)[chrom] - d["t"]["mu"]
    n = len(real)
    rng = np.random.default_rng(3)

    def simulate(kick, m, burn=500):
        x = xs[kick]
        ar = np.r_[1, -pacf_to_coefs(np.tanh(x[:p]))]
        ma = np.r_[1, -pacf_to_coefs(np.tanh(x[p:p + q]))]
        e = sim_kicks(rng, kick, d[kick], m + burn)
        return e[burn:], lfilter(ma, ar, e + d[kick]["mu"])[burn:]

    fig = plt.figure(figsize=(20, 11.5))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1.1], hspace=0.55, wspace=0.12)
    lim = np.abs(real).max() * 1.08
    x = np.arange(n)

    # row 0: real kicks with each model's 99.9% range at every window
    for j, k in enumerate(["real"] + KICKS):
        ax = fig.add_subplot(gs[0, j])
        ax.vlines(x, 0, real, color="#b9b8b3" if j else TEXT, lw=0.5)
        if k == "real":
            title(ax, "Real kicks along chromosome 7", "the surprise at each window after the correlation part")
            ax.set_ylabel("kick (change in log depth)", fontsize=9, color=MUTED)
        else:
            e = arma_kicks(ys, xs[k], p, q)[chrom] - d[k]["mu"]      # this model's own kicks
            if k == "Gaussian":
                half = np.full(n, stats.norm.ppf(0.999) * d[k]["sd"])
            elif k == "t":
                half = np.full(n, stats.t.ppf(0.999, d[k]["nu"]) * d[k]["scale"])
            else:
                half = stats.t.ppf(0.999, d[k]["nu"]) * np.sqrt(garch_s2(e, d[k]["w"], d[k]["a"], d[k]["b"]))
            ax.fill_between(x, -half, half, color=KCOL[k], alpha=0.22, lw=0)
            ax.plot(x, half, color=KCOL[k], lw=0.8)
            ax.plot(x, -half, color=KCOL[k], lw=0.8)
            out = np.abs(e) > half
            ax.scatter(x[out], e[out], s=10, color=TEXT, zorder=3)
            sub = {"Gaussian": "same range everywhere", "t": "same range everywhere, wider tails",
                   "GARCH t": "range widens after big kicks, then settles"}[k]
            title(ax, f"{k}: {out.sum()} of {n:,} kicks outside its 99.8% range", f"{sub}; expect ~{0.002 * n:.0f} outside")
            ax.set_yticklabels([])
        ax.set_ylim(-lim, lim)
        tidy(ax, "y")

    # row 1: real depth and depth simulated from each model, same length
    sims = {k: simulate(k, n) for k in KICKS}
    y_real = ys[chrom][SKIP:]
    top = max(to_depth(y_real).max(), *(to_depth(s[1]).max() for s in sims.values())) * 1.3 + 1
    for j, k in enumerate(["real"] + KICKS):
        ax = fig.add_subplot(gs[1, j])
        y = y_real if k == "real" else sims[k][1]
        ax.plot(x, to_depth(y) + 1, color=TEXT if k == "real" else KCOL[k], lw=0.6)
        ax.set_yscale("log")
        ax.set_ylim(0.8, top)
        ax.set_yticks([1, 11, 101, 1001, 10001], ["0", "10", "100", "1,000", "10,000"])
        ax.minorticks_off()
        if k == "real":
            ax.set_ylabel("read depth (log scale)", fontsize=9, color=MUTED)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("500 bp window along chromosome 7", fontsize=8.5, color=MUTED)
        title(ax, "Real depth" if k == "real" else f"Simulated from the {k} fit")
        tidy(ax, "y")

    # row 2: bunching of big kicks, and the tails
    long = {k: sim_kicks(rng, k, d[k], 200_000) for k in KICKS}
    lags = np.arange(1, 41)
    zs_real = standardised(ys, xs["t"], p, q, "t")
    ax = fig.add_subplot(gs[2, :2])
    ax.plot(lags, [clustering(zs_real, L) for L in lags], color=TEXT, lw=2.5, label="real kicks")
    for k in KICKS:
        ax.plot(lags, [clustering([long[k]], L) for L in lags], color=KCOL[k], lw=2, label=f"simulated, {k}")
    ax.axhline(0, color=TEXT, lw=0.8)
    ax.set_xlabel("distance between the two windows (number of 500 bp windows)", fontsize=8.5, color=MUTED)
    ax.set_ylabel("correlation of |kick|", fontsize=9, color=MUTED)
    title(ax, "Do big kicks bunch together? Real kicks vs kicks simulated from each fit",
          "Gaussian and t kicks never bunch (≈ 0); GARCH t bunches like the real kicks over the nearest windows")
    ax.legend(frameon=False, fontsize=8.5)
    tidy(ax, "y")

    ax = fig.add_subplot(gs[2, 2:])
    allreal = np.concatenate(arma_kicks(ys, xs["t"], p, q)) - d["t"]["mu"]
    lo, hi = np.percentile(allreal, [0.005, 99.995])
    bins = np.linspace(lo, hi, 120)
    mids = (bins[1:] + bins[:-1]) / 2
    dens = lambda v: np.histogram(v, bins=bins, density=True)[0]
    ax.plot(mids, dens(allreal), color=TEXT, lw=2.5, label="real kicks (whole genome)")
    for k in KICKS:
        ax.plot(mids, dens(long[k]), color=KCOL[k], lw=2, label=f"simulated, {k}")
    ax.set_yscale("log")
    ax.set_xlabel("kick (change in log depth)", fontsize=8.5, color=MUTED)
    ax.set_ylabel("density (log scale)", fontsize=9, color=MUTED)
    title(ax, "How often kicks of each size happen (log scale shows the tails)",
          "Gaussian: too few big kicks; t and GARCH t follow the real tails")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    tidy(ax, "y")

    gd = d["GARCH t"]
    fig.suptitle(f"{sample} ({sample_type}, median depth {med:.0f}), {label(p, q)}: Gaussian vs t vs GARCH t kicks.  "
                 f"Fitted: Gaussian sd {d['Gaussian']['sd']:.2f} · t scale {d['t']['scale']:.2f}, ν {d['t']['nu']:.1f} · "
                 f"GARCH t a {gd['a']:.2f}, b {gd['b']:.2f}, ν {gd['nu']:.1f}",
                 x=0.01, ha="left", fontsize=12.5, color=TEXT)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.9, bottom=0.06)
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    kc = pd.read_csv(f"{GARCH_OUT}/kick_compare.tsv", sep="\t")
    os.makedirs(f"{OUT}/figures", exist_ok=True)
    evidence(kc, f"{OUT}/figures/5_kicks_evidence.png")
    xs = chosen_x()
    for tag, t in (("6a", "sWGA"), ("6b", "gDNA")):
        sample = EXAMPLES[t]
        if sample not in xs:                  # not QC-pass: use the sample with the median GARCH t gain
            w = kc[kc["sample type"] == t].pivot(index="sample", columns="kicks", values="Heldout")
            gain = (w["t"] - w["GARCH t"]).sort_values()
            sample = gain.index[len(gain) // 2]
        x = {**xs[sample], "GARCH t": garch_x(sample)}
        explained(sample, t, x, f"{OUT}/figures/{tag}_kicks_explained_{t}.png")
    print(f"wrote 5_kicks_evidence.png, 6a and 6b in {OUT}/figures")


if __name__ == "__main__":
    main()
