"""
Which random process best describes read depth along the genome, for real Pf8 samples?

Fits every candidate process to every sample's log depth (500 bp windows), scores each fit three ways,
and compares sWGA with gDNA. The same models and likelihood as depth_model/depth_process_models.ipynb.

    python compare_processes.py <pf8_dir> <output_dir>

<pf8_dir> holds samples.tsv (Sample, Sample type, ...) and <sample>/depth_plots/<sample>.window_stats.tsv.

Models (each sample gets its own parameters):
  AR(p)        each window = weighted sum of the previous p windows + noise
  ARMA(p,q)    AR plus a weighted sum of past noise terms
  GP <kernel>  correlation is a function of distance: exponential (= AR(1)), Matern 3/2, Matern 5/2,
               squared exponential, two exponentials (a fast and a slow length-scale)
Not included: iid (no correlation; depth is clearly correlated) and hidden Markov models
(their correlation has the same exponential shape as AR(1)).

Scores (all on the -2 log-likelihood scale, so they can be compared; lower = better, 0 = best model of a sample):
  dAIC       AIC = 2k - 2 LL, fitted on all chromosomes
  dBIC       BIC = k ln(n) - 2 LL, fitted on all chromosomes; a larger penalty per parameter than AIC
  dHeldout   -2 x LL on chromosomes the model was not fitted to (odd vs even chromosomes, then swapped);
             no penalty needed, as extra parameters only help if they generalise

Writes into <output_dir>:
  fits.tsv      every sample x model: LL, k, n, AIC, BIC, held-out LL and the three scores
  summary.tsv   per sample type x model: wins and median rank under each score, median of each score
  scores.png    median score per model, sWGA against gDNA
"""

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.signal import lfilter
from statsmodels.tsa.stattools import pacf

WINDOW = 500      # bp per window, as made by krsim.output.pf8_profiles
PSEUDO = 0.05     # added before taking the log, so windows with zero depth do not give -inf
TRIM_KB = 50      # drop this many kb from each end of every chromosome (subtelomeric dropouts)
K = 200           # each window is predicted from at most this many previous windows (100 kb)
AR_ORDERS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]
ARMA_ORDERS = [(1, 1), (2, 1), (2, 2)]
N_PSI = 20000     # psi weights kept for the ARMA autocovariance

KERNELS = {
    "exponential": lambda d, l: np.exp(-d / l),
    "Matern 3/2": lambda d, l: (1 + np.sqrt(3) * d / l) * np.exp(-np.sqrt(3) * d / l),
    "Matern 5/2": lambda d, l: (1 + np.sqrt(5) * d / l + 5 * d ** 2 / (3 * l ** 2)) * np.exp(-np.sqrt(5) * d / l),
    "squared exp": lambda d, l: np.exp(-d ** 2 / (2 * l ** 2)),
}
LAG_D = np.arange(K + 1)


# --- data ---

def load(path):
    """log(depth / sample median + PSEUDO) per nuclear chromosome, ends trimmed, the sample's mean removed."""
    df = pd.read_csv(path, sep="\t")
    df = df[~df["chrom"].str.contains("API|MIT")]
    y_all = np.log(df["mean"] / df["mean"].median() + PSEUDO)
    trim = TRIM_KB * 1000 // WINDOW
    ys = [y_all[df["chrom"] == c].to_numpy()[trim:(df["chrom"] == c).sum() - trim] for c in df["chrom"].unique()]
    mu = np.concatenate(ys).mean()
    return [y - mu for y in ys]


# --- likelihood (Durbin-Levinson, each window predicted from up to K previous ones) ---

def levinson(gamma=None, pacf_values=None, v0=None):
    """Prediction coefficients phi[m] and error variances v[m] from the autocovariance or the PACF."""
    if gamma is not None:
        n_lags, v0 = len(gamma) - 1, gamma[0]
    else:
        n_lags = len(pacf_values)
    phi = np.zeros((n_lags + 1, n_lags + 1))
    v = np.empty(n_lags + 1)
    v[0] = v0
    for m in range(1, n_lags + 1):
        if gamma is not None:
            a = (gamma[m] - phi[m - 1, 1:m] @ gamma[m - 1:0:-1]) / v[m - 1]
        else:
            a = pacf_values[m - 1]
        phi[m, m] = a
        phi[m, 1:m] = phi[m - 1, 1:m] - a * phi[m - 1, m - 1:0:-1]
        v[m] = v[m - 1] * (1 - a * a)
    return phi, v


def lag_matrix(y):
    """Row t = the windows before window t, nearest first (only the first K rows are needed)."""
    L = np.zeros((K, K))
    for t in range(1, K):
        L[t, :t] = y[t - 1::-1][:t]
    return L


def loglik(ys, Ls, phi, v):
    if not (np.all(np.isfinite(v)) and np.all(v > 0)):
        return -np.inf
    ll = 0.0
    for y, L in zip(ys, Ls):
        e = lfilter(np.r_[1, -phi[K, 1:]], [1], y)      # prediction errors once K windows are available
        e[:K] = y[:K] - (phi[:K, 1:] * L).sum(1)        # the first K windows have fewer before them
        var = np.r_[v[:K], np.full(len(y) - K, v[K])]
        ll += -0.5 * np.sum(np.log(2 * np.pi * var) + e * e / var)
    return ll


def fit(ys, Ls, build, starts):
    """Maximum likelihood from several starting values; returns the best parameters and LL."""
    def neg_ll(x):
        ll = loglik(ys, Ls, *build(x))
        return -ll if np.isfinite(ll) else 1e20
    best = None
    for x0 in starts:
        r = minimize(neg_ll, x0, method="L-BFGS-B", bounds=[(-12, 12)] * len(x0))
        if best is None or r.fun < best.fun:
            best = r
    return best.x, -best.fun


# --- models: build(x) -> (phi, v), start(ys) -> starting values ---

def ar_model(p):
    def build(x):
        a = np.zeros(K)
        a[:p] = np.tanh(x[1:])
        return levinson(pacf_values=a, v0=np.exp(x[0]))
    def start(ys):
        y = np.concatenate(ys)
        return [np.r_[np.log(y.var()), np.arctanh(np.clip(pacf(y, nlags=p)[1:], -0.99, 0.99))]]
    return build, start


def pacf_to_coefs(r):
    c = np.zeros(0)
    for a in r:
        c = np.r_[c - a * c[::-1], a]
    return c


def arma_gamma(x, p, q, lags):
    sigma2 = np.exp(x[0])
    ar = np.r_[1, -pacf_to_coefs(np.tanh(x[1:1 + p]))]
    ma = np.r_[1, -pacf_to_coefs(np.tanh(x[1 + p:]))]
    impulse = np.zeros(N_PSI)
    impulse[0] = 1
    psi = lfilter(ma, ar, impulse)
    return sigma2 * np.array([psi[:N_PSI - k] @ psi[k:] for k in lags])


def arma_model(p, q):
    build = lambda x: levinson(gamma=arma_gamma(x, p, q, LAG_D))
    def start(ys):
        y = np.concatenate(ys)
        lv = np.log(y.var())
        slow = np.r_[np.arctanh(0.98), np.zeros(p - 1)]
        return [np.r_[lv - 2, slow, np.full(q, np.arctanh(-0.5))],
                np.r_[lv - 2, slow, np.full(q, np.arctanh(0.5))],
                np.r_[lv, np.arctanh(np.clip(pacf(y, nlags=p)[1:], -0.99, 0.99)), np.zeros(q)]]
    return build, start


def gp_model(kernel):
    def build(x):
        s2, l, t2 = np.exp(x)                 # variance, length-scale (windows), nugget
        g = s2 * KERNELS[kernel](LAG_D, l)
        g[0] += t2
        return levinson(gamma=g)
    def start(ys):
        lv = np.log(np.concatenate(ys).var())
        return [np.r_[lv, np.log(l), lv - 3] for l in (2, 10, 40)]
    return build, start


def gp2_model():
    def build(x):
        s1, l1, s2, l2, t2 = np.exp(x)        # fast part, slow part, nugget
        g = s1 * np.exp(-LAG_D / l1) + s2 * np.exp(-LAG_D / l2)
        g[0] += t2
        return levinson(gamma=g)
    def start(ys):
        lv = np.log(np.concatenate(ys).var())
        return [np.r_[lv - 0.7, np.log(l1), lv - 0.7, np.log(l2), lv - 3] for l1, l2 in ((1, 20), (3, 60), (5, 150))]
    return build, start


def models():
    """name -> (build, start, number of parameters, counting the mean removed in load())"""
    out = {}
    for p in AR_ORDERS:
        out[f"AR({p})"] = (*ar_model(p), p + 2)
    for p, q in ARMA_ORDERS:
        out[f"ARMA({p},{q})"] = (*arma_model(p, q), p + q + 2)
    for kernel in KERNELS:
        out[f"GP {kernel}"] = (*gp_model(kernel), 4)
    out["GP two exponentials"] = (*gp2_model(), 6)
    return out


# --- one sample: full fit and held-out fit for every model ---

def fit_sample(args):
    sample, sample_type, path = args
    t0 = time.time()
    ys = load(path)
    Ls = [lag_matrix(y) for y in ys]
    n = sum(len(y) for y in ys)
    halves = [list(range(0, len(ys), 2)), list(range(1, len(ys), 2))]   # odd vs even chromosomes
    rows = []
    for name, (build, start, k) in models().items():
        _, ll = fit(ys, Ls, build, start(ys))
        held = 0.0
        for train, test in (halves, halves[::-1]):
            y_tr = [ys[i] for i in train]
            x, _ = fit(y_tr, [Ls[i] for i in train], build, start(y_tr))
            held += loglik([ys[i] for i in test], [Ls[i] for i in test], *build(x))
        rows.append({"sample": sample, "sample type": sample_type, "model": name,
                     "LL": ll, "k": k, "n": n, "held-out LL": held})
    print(f"  {sample} ({sample_type}): {time.time() - t0:.0f} s", flush=True)
    return rows


# --- summary and plot ---

def summarise(res):
    rows = []
    for (t, name), g in res.groupby(["sample type", "model"], sort=False):
        row = {"sample type": t, "model": name, "samples": len(g)}
        for score in ("dAIC", "dBIC", "dHeldout"):
            row[f"best by {score[1:]}"] = int((g[score] == 0).sum())
        for score in ("dAIC", "dBIC", "dHeldout"):
            row[f"median rank ({score[1:]})"] = g[f"rank {score[1:]}"].median()
        for score in ("dAIC", "dBIC", "dHeldout"):
            row[f"median {score}"] = round(g[score].median(), 1)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["sample type", "median rank (BIC)"])


def plot(res, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    TYPE_COLOUR = {"sWGA": "#2a78d6", "gDNA": "#eb6834"}
    TEXT, MUTED, GRID = "#0b0b0b", "#8a8984", "#e6e5e1"
    order = list(dict.fromkeys(res["model"]))
    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
    for ax, score, label in zip(axes, ("dAIC", "dBIC", "dHeldout"),
                                ("ΔAIC", "ΔBIC", "Δ held-out (−2 × LL on unseen chromosomes)")):
        for i, (t, colour) in enumerate(TYPE_COLOUR.items()):
            g = res[res["sample type"] == t]
            if g.empty:
                continue
            x = np.arange(len(order)) + (i - 0.5) * 0.25
            for j, name in enumerate(order):
                v = g.loc[g["model"] == name, score]
                ax.scatter(np.full(len(v), x[j]), v + 1, s=10, color=colour, alpha=0.35, lw=0)
            med = g.groupby("model")[score].median().reindex(order)
            ax.scatter(x, med + 1, s=45, color=colour, marker="_", lw=2.5,
                       label=f"{t} ({g['sample'].nunique()} samples; dots = samples, bar = median)")
        ax.set_yscale("log")
        ax.set_ylabel(f"{label} + 1", fontsize=9, color=TEXT)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelsize=8, labelcolor=TEXT)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    axes[-1].set_xticks(range(len(order)), order, rotation=45, ha="right")
    fig.suptitle("Which process fits read depth best? Distance from each sample's best model (lower = better, log scale)",
                 x=0.01, ha="left", fontsize=11, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_png, dpi=150, facecolor="white")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pf8_dir")
    parser.add_argument("output_dir")
    parser.add_argument("--types", nargs="+", default=["sWGA", "gDNA"])
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    info = pd.read_csv(os.path.join(args.pf8_dir, "samples.tsv"), sep="\t")
    jobs = []
    for s, t in zip(info["Sample"], info["Sample type"]):
        path = os.path.join(args.pf8_dir, s, "depth_plots", f"{s}.window_stats.tsv")
        if t in args.types and os.path.exists(path):
            jobs.append((s, t, path))
    print(f"{len(jobs)} samples, {len(models())} models each, {args.workers} at a time")

    with ProcessPoolExecutor(args.workers) as pool:
        res = pd.DataFrame([r for rows in pool.map(fit_sample, jobs) for r in rows])

    res["AIC"] = 2 * res["k"] - 2 * res["LL"]
    res["BIC"] = res["k"] * np.log(res["n"]) - 2 * res["LL"]
    res["Heldout"] = -2 * res["held-out LL"]
    for score in ("AIC", "BIC", "Heldout"):
        res[f"d{score}"] = res[score] - res.groupby("sample")[score].transform("min")
        res[f"rank {score}"] = res.groupby("sample")[f"d{score}"].rank(method="min")

    res.round(2).to_csv(os.path.join(args.output_dir, "fits.tsv"), sep="\t", index=False)
    summary = summarise(res)
    summary.to_csv(os.path.join(args.output_dir, "summary.tsv"), sep="\t", index=False)
    plot(res, os.path.join(args.output_dir, "scores.png"))
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summary.to_string(index=False))
    print(f"wrote fits.tsv, summary.tsv and scores.png in {args.output_dir}")


if __name__ == "__main__":
    main()
