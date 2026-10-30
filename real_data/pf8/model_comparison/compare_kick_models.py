"""
Better kicks for sWGA depth: keep AR(p) for how depth relates to nearby windows, change how the kicks are drawn.

compare_processes.py found AR(~8) captures the correlation of sWGA log depth, but its leftover kicks are
not Gaussian: far too many extreme jumps (heavy tails), and big kicks cluster (calm and jumpy stretches).
This fits AR(p) with five kick models and compares them by held-out likelihood and BIC:

  Gaussian      kick ~ Normal(0, sigma^2), the same sigma everywhere (the model used so far)
  t             kick ~ Student t (heavier tails), the same scale everywhere
  GARCH         kick ~ Normal(0, sigma_t^2), sigma_t^2 = w + a * kick_{t-1}^2 + b * sigma_{t-1}^2:
                a big kick makes the next ones bigger, then it calms down
  GARCH t       GARCH with Student t kicks
  2-state sigma the genome switches between a calm and a jumpy sigma (hidden Markov model on sigma)

    python compare_kick_models.py <pf8_dir> <output_dir>      # --type sWGA (default) --p 8 --workers N

All models use the same conditional likelihood: each chromosome's first p windows are only used as
predictors, so every model is scored on exactly the same windows. Held-out = fit on odd chromosomes, score
even ones, then swap. Writes kick_fits.tsv, kick_summary.tsv and kick_models.png into <output_dir>.
"""

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from numba import njit
from scipy import stats
from scipy.optimize import minimize
from scipy.signal import lfilter
from scipy.special import expit, gammaln

from compare_processes import load

LOG2PI = np.log(2 * np.pi)


# --- data: per chromosome, the targets y_t and the p previous windows ---

def design(ys, p):
    out = []
    for y in ys:
        X = np.column_stack([y[p - i - 1:len(y) - i - 1] for i in range(p)])   # column i = i+1 windows back
        out.append((y[p:], X))
    return out


def kicks(chroms, phi):
    return [y - X @ phi for y, X in chroms]


# --- kick densities ---

def gauss_logpdf(e, s2):
    return -0.5 * (LOG2PI + np.log(s2) + e * e / s2)


def t_logpdf(e, s2, nu):
    """Student t with variance s2 (so sigma means the same as for the Gaussian), nu > 2."""
    scale2 = s2 * (nu - 2) / nu
    return (gammaln((nu + 1) / 2) - gammaln(nu / 2) - 0.5 * np.log(nu * np.pi * scale2)
            - (nu + 1) / 2 * np.log1p(e * e / (nu * scale2)))


def garch_s2(e, w, a, b):
    """sigma_t^2 = w + a e_{t-1}^2 + b sigma_{t-1}^2, starting at the long-run level w / (1 - a - b)."""
    u = w / (1 - a - b)
    x = np.empty_like(e)
    x[0] = w + a * u
    x[1:] = w + a * e[:-1] ** 2
    return lfilter([1.0], [1.0, -b], x, zi=[b * u])[0]


@njit(cache=True)
def hmm_loglik(logd, p11, p22):
    """Forward algorithm for a 2-state switch; logd[t, k] = log density of kick t in state k."""
    stay = np.array([p11, p22])
    pi1 = (1 - p22) / (2 - p11 - p22)
    a0, a1 = pi1, 1 - pi1
    ll = 0.0
    for t in range(logd.shape[0]):
        m = max(logd[t, 0], logd[t, 1])
        d0, d1 = np.exp(logd[t, 0] - m), np.exp(logd[t, 1] - m)
        pr0 = a0 * stay[0] + a1 * (1 - stay[1])
        pr1 = a0 * (1 - stay[0]) + a1 * stay[1]
        f0, f1 = pr0 * d0, pr1 * d1
        c = f0 + f1
        ll += np.log(c) + m
        a0, a1 = f0 / c, f1 / c
    return ll


# --- the five models: (number of extra parameters beyond phi, log-likelihood of kicks, starting values) ---

def nu_of(x):
    return 2.05 + np.exp(x)


def ll_gauss(es, x):
    s2 = np.exp(x[0])
    return sum(gauss_logpdf(e, s2).sum() for e in es)


def ll_t(es, x):
    s2, nu = np.exp(x[0]), nu_of(x[1])
    return sum(t_logpdf(e, s2, nu).sum() for e in es)


def garch_params(x):
    persist = 0.999 * expit(x[1])
    share = expit(x[2])
    return np.exp(x[0]), persist * share, persist * (1 - share)


def ll_garch(es, x):
    w, a, b = garch_params(x)
    return sum(gauss_logpdf(e, garch_s2(e, w, a, b)).sum() for e in es)


def ll_garch_t(es, x):
    w, a, b = garch_params(x)
    nu = nu_of(x[3])
    return sum(t_logpdf(e, garch_s2(e, w, a, b), nu).sum() for e in es)


def ll_switch(es, x):
    s2a = np.exp(x[0])
    s2b = s2a * (1 + np.exp(x[1]))           # the jumpy state has the larger sigma
    p11, p22 = expit(x[2]), expit(x[3])
    total = 0.0
    for e in es:
        logd = np.column_stack([gauss_logpdf(e, s2a), gauss_logpdf(e, s2b)])
        total += hmm_loglik(logd, p11, p22)
    return total


def starts(name, v):
    lv = np.log(v)
    if name == "Gaussian":
        return [[lv]]
    if name == "t":
        return [[lv, np.log(4)], [lv, np.log(1)]]
    if name == "GARCH":
        return [[lv + np.log(1 - 0.95), 3.0, -1.5], [lv + np.log(1 - 0.8), 1.4, -0.5]]
    if name == "GARCH t":
        return [[lv + np.log(1 - 0.95), 3.0, -1.5, np.log(6)], [lv + np.log(1 - 0.8), 1.4, -0.5, np.log(3)]]
    if name == "2-state sigma":
        return [[lv - 0.7, np.log(3), 4.0, 3.0], [lv - 0.3, np.log(1.5), 3.0, 3.0]]


KICK_MODELS = {"Gaussian": ll_gauss, "t": ll_t, "GARCH": ll_garch, "GARCH t": ll_garch_t, "2-state sigma": ll_switch}


def fit(chroms, name, p):
    """Maximum likelihood for AR(p) coefficients + kick-model parameters, from several starting values."""
    y = np.concatenate([c[0] for c in chroms])
    X = np.vstack([c[1] for c in chroms])
    phi0 = np.linalg.lstsq(X, y, rcond=None)[0]
    v = np.mean((y - X @ phi0) ** 2)
    ll_kicks = KICK_MODELS[name]

    def neg(z):
        ll = ll_kicks(kicks(chroms, z[:p]), z[p:])
        return -ll if np.isfinite(ll) else 1e20

    best = None
    for s in starts(name, v):
        r = minimize(neg, np.r_[phi0, s], method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    return best.x, -best.fun


def standardised_kicks(chroms, name, z, p):
    """Kicks divided by the model's sigma at each window (should then look like the model's kick distribution)."""
    es = kicks(chroms, z[:p])
    x = z[p:]
    out = []
    for e in es:
        if name in ("Gaussian", "t"):
            s2 = np.exp(x[0])
        elif name in ("GARCH", "GARCH t"):
            s2 = garch_s2(e, *garch_params(x))
        else:                                    # 2-state: sigma of the more likely state at each window (approximate)
            s2a, s2b = np.exp(x[0]), np.exp(x[0]) * (1 + np.exp(x[1]))
            s2 = np.where(e * e > np.sqrt(s2a * s2b), s2b, s2a)
        out.append(e / np.sqrt(s2))
    return np.concatenate(out), out


def fit_sample(args):
    sample, path, p = args
    t0 = time.time()
    ys = load(path)
    chroms = design(ys, p)
    n = sum(len(c[0]) for c in chroms)
    halves = [list(range(0, len(chroms), 2)), list(range(1, len(chroms), 2))]
    rows = []
    for name, ll_kicks in KICK_MODELS.items():
        z, ll = fit(chroms, name, p)
        k = len(z) + 1                               # + the mean removed in load()
        held = 0.0
        for train, test in (halves, halves[::-1]):
            zt, _ = fit([chroms[i] for i in train], name, p)
            held += ll_kicks(kicks([chroms[i] for i in test], zt[:p]), zt[p:])
        zc, zs = standardised_kicks(chroms, name, z, p)
        extra = {}
        if name in ("t", "GARCH t"):
            extra["nu (t tail weight)"] = nu_of(z[-1])
        if name in ("GARCH", "GARCH t"):
            w, a, b = garch_params(z[p:])
            extra["GARCH a"], extra["GARCH b"] = a, b
        if name == "2-state sigma":
            s2a = np.exp(z[p]); s2b = s2a * (1 + np.exp(z[p + 1]))
            p11, p22 = expit(z[p + 2]), expit(z[p + 3])
            extra.update({"sigma calm": np.sqrt(s2a), "sigma jumpy": np.sqrt(s2b),
                          "share jumpy": (1 - p11) / (2 - p11 - p22),
                          "jumpy stretch (kb)": 0.5 / (1 - p22), "calm stretch (kb)": 0.5 / (1 - p11)})
        rows.append({"sample": sample, "model": name, "k": k, "n": n, "LL": ll, "held-out LL": held,
                     "std kick kurtosis": stats.kurtosis(zc),
                     "std |kick| corr lag 1": np.mean([np.corrcoef(np.abs(e[:-1]), np.abs(e[1:]))[0, 1] for e in zs]),
                     **extra})
    print(f"  {sample}: {time.time() - t0:.0f} s", flush=True)
    return rows


def plot(res, out_png, sample_type, p):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    TEXT, MUTED, GRID = "#0b0b0b", "#8a8984", "#e6e5e1"
    COL = {"Gaussian": "#b9b8b3", "t": "#2a78d6", "GARCH": "#1baf7a", "GARCH t": "#eb6834", "2-state sigma": "#8a4fd6"}
    order = list(KICK_MODELS)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    for ax, col, label in ((axes[0], "dHeldout", "Held-out score − best (lower = better)"),
                           (axes[1], "dBIC", "BIC − best (lower = better)")):
        for j, name in enumerate(order):
            v = res.loc[res["model"] == name, col]
            ax.scatter(np.full(len(v), j) + np.linspace(-0.15, 0.15, len(v)), v + 1, s=14, color=COL[name], alpha=0.6, lw=0)
            ax.scatter(j, v.median() + 1, s=260, marker="_", color=COL[name], lw=3)
            ax.text(j + 0.22, v.median() + 1, f"{v.median():,.0f}", va="center", fontsize=8, color=TEXT)
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,g}"))
        ax.set_xticks(range(len(order)), order, fontsize=8.5)
        ax.set_title(f"{label} + 1, log scale\n(dots = samples, bar = median)", loc="left", fontsize=10, color=TEXT)
    ax = axes[2]
    for name in order:
        g = res[res["model"] == name]
        ax.scatter(g["std kick kurtosis"], g["std |kick| corr lag 1"], s=30, color=COL[name], label=name, alpha=0.85, lw=0)
    ax.axvline(0, color=MUTED, lw=0.8, ls=":")
    ax.axhline(0, color=MUTED, lw=0.8, ls=":")
    ax.set_xlabel("Excess kurtosis of standardised kicks (0 = Gaussian; t models expect > 0)", fontsize=8.5)
    ax.set_ylabel("Do big kicks cluster? (0 = no)", fontsize=8.5)
    ax.set_title("What is left over after each model\n(one dot per sample; nearer (0, 0) = less unexplained)", loc="left", fontsize=10, color=TEXT)
    ax.legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle(f"{sample_type}: AR({p}) with five ways of drawing the kicks ({res['sample'].nunique()} samples)",
                 x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pf8_dir")
    parser.add_argument("output_dir")
    parser.add_argument("--type", default="sWGA")
    parser.add_argument("--p", type=int, default=8)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    info = pd.read_csv(os.path.join(args.pf8_dir, "samples.tsv"), sep="\t")
    jobs = [(s, os.path.join(args.pf8_dir, s, "depth_plots", f"{s}.window_stats.tsv"), args.p)
            for s, t in zip(info["Sample"], info["Sample type"]) if t == args.type]
    jobs = [j for j in jobs if os.path.exists(j[1])]
    print(f"{len(jobs)} {args.type} samples, AR({args.p}) with {len(KICK_MODELS)} kick models, {args.workers} at a time")
    with ProcessPoolExecutor(args.workers) as pool:
        res = pd.DataFrame([r for rows in pool.map(fit_sample, jobs) for r in rows])

    res["BIC"] = res["k"] * np.log(res["n"]) - 2 * res["LL"]
    res["Heldout"] = -2 * res["held-out LL"]
    for score in ("BIC", "Heldout"):
        res[f"d{score}"] = res[score] - res.groupby("sample")[score].transform("min")
    res.round(4).to_csv(os.path.join(args.output_dir, "kick_fits.tsv"), sep="\t", index=False)

    summary = res.groupby("model", sort=False).agg(
        best_heldout=("dHeldout", lambda v: int((v == 0).sum())),
        median_dHeldout=("dHeldout", "median"),
        best_BIC=("dBIC", lambda v: int((v == 0).sum())),
        median_dBIC=("dBIC", "median"),
        median_kurtosis=("std kick kurtosis", "median"),
        median_kick_clustering=("std |kick| corr lag 1", "median"),
    ).round(2)
    summary.to_csv(os.path.join(args.output_dir, "kick_summary.tsv"), sep="\t")
    plot(res, os.path.join(args.output_dir, "kick_models.png"), args.type, args.p)
    print(summary.to_string())
    extra = [c for c in res.columns if c in ("nu (t tail weight)", "GARCH a", "GARCH b", "sigma calm", "sigma jumpy",
                                             "share jumpy", "jumpy stretch (kb)", "calm stretch (kb)")]
    print(res.groupby("model", sort=False)[extra].median().round(3).to_string())
    print(f"wrote kick_fits.tsv, kick_summary.tsv and kick_models.png in {args.output_dir}")


if __name__ == "__main__":
    main()
