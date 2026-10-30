"""
The simplest good model of read depth along the genome, for each sample type.

Candidates (the same for gDNA and sWGA): correlation AR(1), AR(2), AR(4), AR(8), ARMA(1,1), ARMA(2,1), ARMA(2,2),
each with Gaussian or Student t kicks (t = one extra parameter, nu, for heavier tails). Scored on the same
windows by held-out likelihood (fit on odd chromosomes, score even, then swap) and BIC, and checked for
calibration (do the kicks land where the model's own distribution says they should?).

    python best_simple_model.py <pf8_dir> <output_dir>          # --workers N

Writes simple_fits.tsv, simple_summary.tsv, simple_calibration.tsv and simple_models.png into <output_dir>.
"""

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

from compare_arma_kicks import arma_kicks, kick_cdf, load, score, unpack

CORR = [(1, 0), (2, 0), (4, 0), (8, 0), (1, 1), (2, 1), (2, 2)]
KICKS = ["Gaussian", "t"]
CLOSE = 10            # held-out score within this of the sample's best counts as "as good"


def label(p, q):
    return f"AR({p})" if q == 0 else f"ARMA({p},{q})"


def fit_sample(job):
    sample, sample_type, path = job
    t0 = time.time()
    ys = load(path)
    rows, calib = [], []
    for p, q in CORR:
        start = None
        for kick in KICKS:
            x, r = score(ys, p, q, kick, start)
            start = x[:p + q]                    # the t fit starts from the Gaussian fit's correlation part
            th = x[p + q:]
            u = np.concatenate([kick_cdf(kick, e, th) for e in arma_kicks(ys, x, p, q)])
            name = f"{label(p, q)} {kick}"
            d = unpack(kick, th)
            rows.append({"sample": sample, "sample type": sample_type, "model": name, "correlation": label(p, q),
                         "kicks": kick, **r, "nu": d.get("nu", np.nan)})
            calib.append({"sample": sample, "sample type": sample_type, "model": name,
                          "lowest 0.1% (expect 0.1)": 100 * np.mean(u < 0.001),
                          "highest 0.1% (expect 0.1)": 100 * np.mean(u > 0.999),
                          "lowest 1% (expect 1)": 100 * np.mean(u < 0.01),
                          "highest 1% (expect 1)": 100 * np.mean(u > 0.99),
                          "KS distance": stats.kstest(u, "uniform").statistic})
    print(f"  {sample} ({sample_type}): {time.time() - t0:.0f} s", flush=True)
    return rows, calib


def plot(res, calib, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    TEXT, GRID = "#0b0b0b", "#e6e5e1"
    KCOL = {"Gaussian": "#b9b8b3", "t": "#2a78d6"}
    corr = [label(p, q) for p, q in CORR]
    fig, axes = plt.subplots(2, 2, figsize=(15, 9.5), gridspec_kw={"width_ratios": [1.6, 1]})
    for row, t in enumerate(["gDNA", "sWGA"]):
        g = res[res["sample type"] == t]
        ax = axes[row, 0]
        x = np.arange(len(corr))
        for i, kick in enumerate(KICKS):
            med = g[g["kicks"] == kick].groupby("correlation")["dHeldout"].median().reindex(corr)
            bars = ax.bar(x + (i - 0.5) * 0.38, med + 1, width=0.36, color=KCOL[kick], label=f"{kick} kicks")
            for b, v in zip(bars, med):
                ax.text(b.get_x() + b.get_width() / 2, (v + 1) * 1.15, f"{v:,.0f}", ha="center", fontsize=7.5, color=TEXT)
        ax.set_yscale("log")
        ax.set_ylim(0.7, ax.get_ylim()[1] * 8)          # headroom for the legend
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,g}"))
        ax.set_xticks(x, corr, fontsize=9)
        ax.set_ylabel("median held-out score − best, + 1\n(log scale; lower = better)", fontsize=9)
        ax.set_title(f"{t} ({g['sample'].nunique()} samples): every correlation structure with Gaussian or t kicks",
                     loc="left", fontsize=10, color=TEXT)
        ax.legend(frameon=False, fontsize=8.5, loc="upper right", ncol=2)
        ax = axes[row, 1]
        c = calib[calib["sample type"] == t].groupby("model")[["lowest 0.1% (expect 0.1)", "highest 0.1% (expect 0.1)"]].median()
        names = [f"{m} {k}" for m in corr for k in KICKS]
        c = c.reindex(names)
        y = np.arange(len(names))
        ax.barh(y + 0.2, c.iloc[:, 0], height=0.38, color="#eb6834", label="biggest drops (model's lowest 0.1%)")
        ax.barh(y - 0.2, c.iloc[:, 1], height=0.38, color="#1baf7a", label="biggest spikes (model's highest 0.1%)")
        ax.axvline(0.1, color=TEXT, ls="--", lw=1)
        ax.set_yticks(y, names, fontsize=7.5)
        ax.invert_yaxis()
        ax.set_xlabel("% of kicks (median over samples); a fitting model gives 0.1", fontsize=8.5)
        ax.set_title("Are the extreme kicks as common as the model says?", loc="left", fontsize=10, color=TEXT)
        ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    for ax in axes.flat:
        ax.grid(axis="y" if ax in axes[:, 0] else "x", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle("The simplest good model: correlation structure x Gaussian or t kicks", x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)


def replot(output_dir):
    """Redraw simple_models.png from the saved tables."""
    res = pd.read_csv(os.path.join(output_dir, "simple_fits.tsv"), sep="\t")
    calib = pd.read_csv(os.path.join(output_dir, "simple_calibration.tsv"), sep="\t")
    plot(res, calib, os.path.join(output_dir, "simple_models.png"))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pf8_dir")
    parser.add_argument("output_dir")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    info = pd.read_csv(os.path.join(args.pf8_dir, "samples.tsv"), sep="\t")
    jobs = [(s, t, os.path.join(args.pf8_dir, s, "depth_plots", f"{s}.window_stats.tsv"))
            for s, t in zip(info["Sample"], info["Sample type"])]
    jobs = [j for j in jobs if os.path.exists(j[2])]
    print(f"{len(jobs)} samples, {len(CORR) * len(KICKS)} models each")
    with ProcessPoolExecutor(args.workers) as pool:
        out = list(pool.map(fit_sample, jobs))
    res = pd.DataFrame([r for rows, _ in out for r in rows])
    calib = pd.DataFrame([c for _, cs in out for c in cs])
    res["BIC"] = res["k"] * np.log(res["n"]) - 2 * res["LL"]
    res["Heldout"] = -2 * res["held-out LL"]
    for s in ("BIC", "Heldout"):
        res[f"d{s}"] = res[s] - res.groupby("sample")[s].transform("min")
    res["close to best"] = res["dHeldout"] <= CLOSE
    res.round(3).to_csv(os.path.join(args.output_dir, "simple_fits.tsv"), sep="\t", index=False)
    calib.round(3).to_csv(os.path.join(args.output_dir, "simple_calibration.tsv"), sep="\t", index=False)

    summ = res.groupby(["sample type", "model"], sort=False).agg(
        k=("k", "first"), best_heldout=("dHeldout", lambda v: int((v == 0).sum())),
        close_to_best=("close to best", "sum"), median_dHeldout=("dHeldout", "median"),
        best_BIC=("dBIC", lambda v: int((v == 0).sum())), median_dBIC=("dBIC", "median"), median_nu=("nu", "median"),
    ).round(1).reset_index()
    cal = calib.groupby(["sample type", "model"], sort=False).median(numeric_only=True).round(2).reset_index()
    summ = summ.merge(cal, on=["sample type", "model"])
    summ.to_csv(os.path.join(args.output_dir, "simple_summary.tsv"), sep="\t", index=False)
    with pd.option_context("display.width", 260, "display.max_columns", 30):
        for t, g in summ.groupby("sample type"):
            print(f"\n== {t}")
            print(g.drop(columns="sample type").sort_values("median_dHeldout").to_string(index=False))
    plot(res, calib, os.path.join(args.output_dir, "simple_models.png"))
    print(f"wrote simple_fits.tsv, simple_summary.tsv, simple_calibration.tsv and simple_models.png in {args.output_dir}")


if __name__ == "__main__":
    main()
