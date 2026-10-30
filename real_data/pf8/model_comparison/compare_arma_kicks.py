"""
ARMA order and kick distribution for Pf8 depth (written for gDNA; works for any sample type).

Step 1 (correlation): every ARMA(p, q) with p, q = 0..3 (not both 0), Gaussian kicks.
Step 2 (kicks): on the step-1 winner, six ways of drawing the kicks:
  Gaussian          Normal(0, sigma^2)
  t                 Student t (heavier tails)
  skew t            Student t that can be lopsided (Fernandez-Steel); gamma < 1 = bigger drops than spikes
  GARCH t           t whose scale grows after a big kick and settles back
  GARCH skew t      both
  Gauss + dropout   mostly Normal, plus an occasional downward-only drop (a two-part mixture)

    python compare_arma_kicks.py <pf8_dir> <output_dir>          # --type gDNA (default) --workers N

All models use the same conditional likelihood: kicks are computed along each chromosome with a filter,
and the first SKIP windows of every chromosome are not scored, so every model is scored on the same windows.
Held-out = fit on odd chromosomes, score even ones, then swap. BIC = k ln(n) - 2 LL.
Writes order_fits.tsv, kick_fits.tsv, calibration.tsv and arma_kicks.png into <output_dir>.
"""

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from scipy.signal import lfilter
from scipy.special import expit, gammaln, logsumexp

from compare_kick_models import garch_s2
from compare_processes import load, pacf_to_coefs

SKIP = 10                                   # windows at the start of each chromosome used only as a run-in
ORDERS = [(p, q) for p in range(4) for q in range(4) if p + q > 0]
LOG2PI = np.log(2 * np.pi)


# --- kicks from ARMA coefficients (AR and MA parts kept stationary / invertible via partial autocorrelations) ---

def arma_kicks(ys, x, p, q):
    ar = np.r_[1, -pacf_to_coefs(np.tanh(x[:p]))]
    ma = np.r_[1, -pacf_to_coefs(np.tanh(x[p:p + q]))]
    return [lfilter(ar, ma, y)[SKIP:] for y in ys]


# --- kick log-densities; th = the kick model's own parameters ---

def t_logpdf(e, scale, nu):
    z = e / scale
    return (gammaln((nu + 1) / 2) - gammaln(nu / 2) - 0.5 * np.log(nu * np.pi) - np.log(scale)
            - (nu + 1) / 2 * np.log1p(z * z / nu))


def skew_t_logpdf(e, scale, nu, g):
    """Fernandez-Steel skewed t: the right side stretched by g, the left by 1/g."""
    z = np.where(e >= 0, e / g, e * g)
    return np.log(2) - np.log(g + 1 / g) + t_logpdf(z, scale, nu)


def skew_t_cdf(e, scale, nu, g):
    c = stats.t.cdf
    return np.where(e < 0, 2 / (1 + g * g) * c(e * g / scale, nu),
                    1 / (1 + g * g) + 2 * g * g / (1 + g * g) * (c(e / (g * scale), nu) - 0.5))


def unpack(name, th):
    """Kick-model parameters on their natural scale."""
    if name == "Gaussian":
        return {"mu": th[0], "sd": np.exp(th[1])}
    if name in ("t", "skew t"):
        d = {"mu": th[0], "scale": np.exp(th[1]), "nu": 2.05 + np.exp(th[2])}
        if name == "skew t":
            d["gamma"] = np.exp(th[3])
        return d
    if name in ("GARCH t", "GARCH skew t"):
        persist = 0.999 * expit(th[2])
        d = {"mu": th[0], "w": np.exp(th[1]), "a": persist * expit(th[3]), "b": persist * (1 - expit(th[3])),
             "nu": 2.05 + np.exp(th[4])}
        if name == "GARCH skew t":
            d["gamma"] = np.exp(th[5])
        return d
    if name == "Gauss + dropout":
        return {"mu": th[0], "sd": np.exp(th[1]), "drop share": expit(th[2]),
                "drop shift": -np.exp(th[3]), "drop sd": np.exp(th[4])}


def kick_logpdf(name, e, th):
    d = unpack(name, th)
    e = e - d["mu"]
    if name == "Gaussian":
        return -0.5 * (LOG2PI + 2 * np.log(d["sd"]) + (e / d["sd"]) ** 2)
    if name == "t":
        return t_logpdf(e, d["scale"], d["nu"])
    if name == "skew t":
        return skew_t_logpdf(e, d["scale"], d["nu"], d["gamma"])
    if name in ("GARCH t", "GARCH skew t"):
        scale = np.sqrt(garch_s2(e, d["w"], d["a"], d["b"]))
        return t_logpdf(e, scale, d["nu"]) if name == "GARCH t" else skew_t_logpdf(e, scale, d["nu"], d["gamma"])
    if name == "Gauss + dropout":
        main = -0.5 * (LOG2PI + 2 * np.log(d["sd"]) + (e / d["sd"]) ** 2)
        de = e - d["drop shift"]
        drop = -0.5 * (LOG2PI + 2 * np.log(d["drop sd"]) + (de / d["drop sd"]) ** 2)
        return logsumexp(np.column_stack([np.log1p(-d["drop share"]) + main, np.log(d["drop share"]) + drop]), axis=1)


def kick_cdf(name, e, th):
    """Where each kick falls in the model's own distribution (should be uniform if the model fits)."""
    d = unpack(name, th)
    e = e - d["mu"]
    if name == "Gaussian":
        return stats.norm.cdf(e / d["sd"])
    if name == "t":
        return stats.t.cdf(e / d["scale"], d["nu"])
    if name == "skew t":
        return skew_t_cdf(e, d["scale"], d["nu"], d["gamma"])
    if name in ("GARCH t", "GARCH skew t"):
        scale = np.sqrt(garch_s2(e, d["w"], d["a"], d["b"]))
        return stats.t.cdf(e / scale, d["nu"]) if name == "GARCH t" else skew_t_cdf(e, scale, d["nu"], d["gamma"])
    if name == "Gauss + dropout":
        return ((1 - d["drop share"]) * stats.norm.cdf(e / d["sd"])
                + d["drop share"] * stats.norm.cdf((e - d["drop shift"]) / d["drop sd"]))


def kick_starts(name, ls):
    """Starting values for the kick parameters; ls = log of the Gaussian kick SD."""
    if name == "Gaussian":
        return [[0, ls]]
    if name == "t":
        return [[0, ls - 0.2, np.log(4)], [0, ls - 0.05, np.log(15)]]
    if name == "skew t":
        return [[0.05, ls - 0.2, np.log(4), -0.3], [0.1, ls - 0.1, np.log(10), -0.6]]
    if name == "GARCH t":
        lw = 2 * ls + np.log(0.2)
        return [[0, lw, 2.0, -1.0, np.log(5)], [0, lw, 1.0, 0.0, np.log(10)]]
    if name == "GARCH skew t":
        lw = 2 * ls + np.log(0.2)
        return [[0.05, lw, 2.0, -1.0, np.log(5), -0.3], [0.1, lw, 1.0, 0.0, np.log(10), -0.6]]
    if name == "Gauss + dropout":
        return [[0.05, ls - 0.2, -3.5, ls + 1.0, ls + 0.3], [0.02, ls - 0.1, -4.5, ls + 1.5, ls]]


KICKS = ["Gaussian", "t", "skew t", "GARCH t", "GARCH skew t", "Gauss + dropout"]


def fit(ys, p, q, name, arma_start=None):
    """Maximum likelihood for ARMA(p, q) + kick model; returns parameters and LL."""
    y = np.concatenate(ys)
    ls = np.log(y.std())

    def neg(x):
        es = arma_kicks(ys, x, p, q)
        ll = sum(kick_logpdf(name, e, x[p + q:]).sum() for e in es)
        return -ll if np.isfinite(ll) else 1e20

    if arma_start is not None:
        arma_starts = [arma_start]
    else:
        arma_starts = [np.r_[np.full(p, 0.5), np.full(q, 0.3)], np.r_[np.arctanh(np.r_[0.9, np.zeros(p - 1)]) if p else [],
                                                                       np.full(q, -0.3)]]
    best = None
    for a0 in arma_starts:
        for k0 in kick_starts(name, ls - (0.5 if p + q else 0)):
            r = minimize(neg, np.r_[a0, k0], method="L-BFGS-B")
            if best is None or r.fun < best.fun:
                best = r
    return best.x, -best.fun


def score(ys, p, q, name, arma_start=None):
    """Full fit (LL, k, n) and held-out LL (odd vs even chromosomes)."""
    x, ll = fit(ys, p, q, name, arma_start)
    halves = [list(range(0, len(ys), 2)), list(range(1, len(ys), 2))]
    held = 0.0
    for train, test in (halves, halves[::-1]):
        xt, _ = fit([ys[i] for i in train], p, q, name, x[:p + q])
        held += sum(kick_logpdf(name, e, xt[p + q:]).sum() for e in arma_kicks([ys[i] for i in test], xt, p, q))
    n = sum(len(y) - SKIP for y in ys)
    return x, {"LL": ll, "k": len(x) + 1, "n": n, "held-out LL": held}   # + 1 for the mean removed in load()


def step1(job):
    sample, path = job
    ys = load(path)
    rows, best = [], {}
    for p, q in ORDERS:
        x, r = score(ys, p, q, "Gaussian")
        rows.append({"sample": sample, "model": f"ARMA({p},{q})", **r})
        best[(p, q)] = x[:p + q]
    return rows, best


def step2(job):
    sample, path, p, q, arma_start = job
    t0 = time.time()
    ys = load(path)
    rows, calib = [], []
    for name in KICKS:
        x, r = score(ys, p, q, name, arma_start)
        th = x[p + q:]
        u = np.concatenate([kick_cdf(name, e, th) for e in arma_kicks(ys, x, p, q)])
        rows.append({"sample": sample, "model": name, **r, **{k: v for k, v in unpack(name, th).items()}})
        calib.append({"sample": sample, "model": name,
                      "lowest 1% (expect 1)": 100 * np.mean(u < 0.01), "highest 1% (expect 1)": 100 * np.mean(u > 0.99),
                      "lowest 0.1% (expect 0.1)": 100 * np.mean(u < 0.001), "highest 0.1% (expect 0.1)": 100 * np.mean(u > 0.999),
                      "middle 50% (expect 50)": 100 * np.mean((u > 0.25) & (u < 0.75)),
                      "KS distance": stats.kstest(u, "uniform").statistic})
    print(f"  {sample}: {time.time() - t0:.0f} s", flush=True)
    return rows, calib


def add_scores(res):
    res["BIC"] = res["k"] * np.log(res["n"]) - 2 * res["LL"]
    res["Heldout"] = -2 * res["held-out LL"]
    for s in ("BIC", "Heldout"):
        res[f"d{s}"] = res[s] - res.groupby("sample")[s].transform("min")
    return res


def summary(res):
    return res.groupby("model", sort=False).agg(
        best_heldout=("dHeldout", lambda v: int((v == 0).sum())), median_dHeldout=("dHeldout", "median"),
        best_BIC=("dBIC", lambda v: int((v == 0).sum())), median_dBIC=("dBIC", "median"), k=("k", "first")).round(1)


def plot(order_res, kick_res, calib, out_png, sample_type, winner):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    TEXT, MUTED, GRID, COL = "#0b0b0b", "#8a8984", "#e6e5e1", "#eb6834"
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw={"width_ratios": [1.15, 1, 1]})

    # 1: ARMA order grid, median held-out score
    ax = axes[0]
    grid = np.full((4, 4), np.nan)
    med = order_res.groupby("model")["dHeldout"].median()
    for (p, q) in ORDERS:
        grid[p, q] = med[f"ARMA({p},{q})"]
    im = ax.imshow(np.log10(grid + 1), cmap="Oranges_r", origin="lower")
    for p in range(4):
        for q in range(4):
            if not np.isnan(grid[p, q]):
                ax.text(q, p, f"{grid[p, q]:,.0f}", ha="center", va="center", fontsize=9,
                        color="white" if grid[p, q] < 30 else TEXT, weight="bold" if (p, q) == winner else None)
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xlabel("q: previous kicks remembered (MA part)", fontsize=9)
    ax.set_ylabel("p: previous depths looked back at (AR part)", fontsize=9)
    ax.set_title(f"Step 1: ARMA(p, q), Gaussian kicks\nmedian held-out score − best (lower = better); best = ARMA{winner}",
                 loc="left", fontsize=10, color=TEXT)

    # 2: kick models, held-out score per sample
    ax = axes[1]
    for j, name in enumerate(KICKS):
        v = kick_res.loc[kick_res["model"] == name, "dHeldout"]
        ax.scatter(np.full(len(v), j) + np.linspace(-0.15, 0.15, len(v)), v + 1, s=14, color=COL, alpha=0.5, lw=0)
        ax.scatter(j, v.median() + 1, s=260, marker="_", color=COL, lw=3)
        ax.text(j + 0.22, v.median() + 1, f"{v.median():,.0f}", va="center", fontsize=8, color=TEXT)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,g}"))
    ax.set_xticks(range(len(KICKS)), KICKS, rotation=30, ha="right", fontsize=8.5)
    ax.set_title(f"Step 2: ARMA{winner} with six kick models\nheld-out score − best + 1 (dots = samples, bar = median)",
                 loc="left", fontsize=10, color=TEXT)

    # 3: calibration: share of kicks in the lowest / highest 1%
    ax = axes[2]
    c = calib.groupby("model", sort=False)[["lowest 1% (expect 1)", "highest 1% (expect 1)"]].median().reindex(KICKS)
    x = np.arange(len(KICKS))
    ax.bar(x - 0.2, c.iloc[:, 0], width=0.38, color="#2a78d6", label="kicks in the model's lowest 1% (big drops)")
    ax.bar(x + 0.2, c.iloc[:, 1], width=0.38, color="#1baf7a", label="kicks in the model's highest 1% (big spikes)")
    ax.axhline(1, color=TEXT, ls="--", lw=1)
    ax.text(len(KICKS) - 0.5, 1.03, "a fitting model: 1%", ha="right", va="bottom", fontsize=8, color=TEXT)
    ax.set_xticks(x, KICKS, rotation=30, ha="right", fontsize=8.5)
    ax.set_ylabel("% of kicks (median over samples)", fontsize=9)
    ax.set_title("Does each model get the big drops and spikes right?\n(both bars at 1% = right)", loc="left", fontsize=10, color=TEXT)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    for ax in axes[1:]:
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle(f"{sample_type}: ARMA order and kick distribution ({kick_res['sample'].nunique()} samples)",
                 x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pf8_dir")
    parser.add_argument("output_dir")
    parser.add_argument("--type", default="gDNA")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    info = pd.read_csv(os.path.join(args.pf8_dir, "samples.tsv"), sep="\t")
    jobs = [(s, os.path.join(args.pf8_dir, s, "depth_plots", f"{s}.window_stats.tsv"))
            for s, t in zip(info["Sample"], info["Sample type"]) if t == args.type]
    jobs = [j for j in jobs if os.path.exists(j[1])]

    print(f"step 1: {len(jobs)} {args.type} samples, {len(ORDERS)} ARMA orders")
    with ProcessPoolExecutor(args.workers) as pool:
        out = list(pool.map(step1, jobs))
    order_res = add_scores(pd.DataFrame([r for rows, _ in out for r in rows]))
    order_res.round(3).to_csv(os.path.join(args.output_dir, "order_fits.tsv"), sep="\t", index=False)
    s1 = summary(order_res).sort_values("median_dHeldout")
    print(s1.to_string())
    winner = tuple(int(c) for c in s1.index[0][5:-1].split(","))
    print(f"step 1 winner: ARMA{winner}")

    print(f"step 2: ARMA{winner} with {len(KICKS)} kick models")
    starts = {job[0]: best[winner] for job, (_, best) in zip(jobs, out)}
    with ProcessPoolExecutor(args.workers) as pool:
        out2 = list(pool.map(step2, [(s, path, *winner, starts[s]) for s, path in jobs]))
    kick_res = add_scores(pd.DataFrame([r for rows, _ in out2 for r in rows]))
    calib = pd.DataFrame([c for _, cs in out2 for c in cs])
    kick_res.round(4).to_csv(os.path.join(args.output_dir, "kick_fits.tsv"), sep="\t", index=False)
    calib.round(3).to_csv(os.path.join(args.output_dir, "calibration.tsv"), sep="\t", index=False)
    print(summary(kick_res).to_string())
    print(calib.groupby("model", sort=False).median(numeric_only=True).round(2).to_string())
    params = [c for c in ("nu", "gamma", "a", "b", "drop share", "drop shift", "drop sd", "sd") if c in kick_res]
    print(kick_res.groupby("model", sort=False)[params].median().round(3).to_string())
    plot(order_res, kick_res, calib, os.path.join(args.output_dir, "arma_kicks.png"), args.type, winner)
    print(f"wrote order_fits.tsv, kick_fits.tsv, calibration.tsv and arma_kicks.png in {args.output_dir}")


if __name__ == "__main__":
    main()
