"""
Coalescence times between strains in different hosts, varying one parameter at a time.

Two parameter sets (SETS, chosen by COAL_SET): "paper" (baseline Nh 100, K 5, chi 0.5, C 1; the default) and
"real_scale" (baseline Nh 100, K 10, chi 0.2, C 3, as in docs/real-scale.qmd); T 1000 in both. Each of Nh, K, chi
and C takes 5 values, the others held at the baseline; every setting is run with 3 seeds. Results in results/<set>/.

For each run: `krsim simulate`, then, from the tree sequence, the fraction of (pair x genome) coalescing in each
generation, for pairs of strains in different hosts (and, for reference, in the same host), with tskit's
pair_coalescence_counts (as in the real-scale notebook). Cumulated, this is the fraction of the genome coalesced
by each generation. The tree sequence is deleted after use.

    python coal_sweep.py list               # the runs, one per line (task id = line number - 1)
    sbatch --array=0-50 coal_sweep.sbatch   # one run per task -> results/<set>/runs/<param>_<value>_seed<seed>.npz
    python coal_sweep.py plot               # -> results/<set>/coal_sweep*.{png,pdf}, coal_sweep_summary.tsv
    (prefix with COAL_SET=real_scale for the other set; sbatch passes the variable on)
"""

import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
SCRATCH = "/data/home/bty208/scratch/IMPERIAL/krsim_coal_sweep"
KRSIM = os.path.join(HERE, "..", "..", ".venv", "bin", "krsim")

# parameter sets; choose with the environment variable COAL_SET (default "paper")
SETS = {
    # the real-scale tutorial baseline (docs/real-scale.qmd)
    "real_scale": ({"Nh": 100, "T": 1000, "K": 10, "chi": 0.2, "C": 3},
                   {"Nh": [25, 50, 100, 200, 400], "K": [1, 2, 5, 10, 20], "chi": [0, 0.1, 0.2, 0.5, 0.9],
                    "C": [1, 2, 3, 5, 10]}),
    # the values from the user's earlier write-up (Q there = K here)
    "paper": ({"Nh": 100, "T": 1000, "K": 5, "chi": 0.5, "C": 1},
              {"Nh": [50, 100, 150, 200, 300], "K": [3, 5, 10, 25, 50], "chi": [0.05, 0.1, 0.2, 0.5, 0.9],
               "C": [1, 2, 3, 5, 10]}),
}
SET = os.environ.get("COAL_SET", "paper")
BASELINE, SWEEP = SETS[SET]
RESULTS = os.path.join(RESULTS, SET)
SCRATCH = os.path.join(SCRATCH, SET)
SEEDS = [1, 2, 3]
LABEL = {"Nh": "Nh (hosts)", "K": "K (parasites per bite)", "chi": "chi (superinfection)", "C": "C (oocysts)"}


def runs():
    """Every distinct (setting, seed); the baseline is run once and shared by all four sweeps."""
    out, seen = [], set()
    for param, values in SWEEP.items():
        for value in values:
            params = {**BASELINE, param: value}
            key = tuple(params.items())
            if key in seen:
                continue
            seen.add(key)
            name = "baseline" if params == BASELINE else f"{param}_{value}"
            out += [(name, params, seed) for seed in SEEDS]
    return out


def base_text():
    return ", ".join(f"{k} {v:g}" for k, v in BASELINE.items() if k != "T") + f"; T {BASELINE['T']}"


def result_path(name, seed):
    return os.path.join(RESULTS, "runs", f"{name}_seed{seed}.npz")


def run_one(task):
    import pandas as pd
    import tskit

    name, params, seed = runs()[task]
    if os.path.exists(result_path(name, seed)):
        print("already done:", name, seed)
        return
    out = os.path.join(SCRATCH, "runs", f"{name}_seed{seed}")
    os.makedirs(out, exist_ok=True)
    cmd = [KRSIM, "simulate", "--random_seed", str(seed), "--output_dir", out]
    for k, v in params.items():
        cmd += [f"--{k}", str(v)]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

    ts = tskit.load(f"{out}/ts_{seed}.trees")
    node_map = pd.read_csv(f"{out}/sim_ind_node_map_{seed}.csv")
    hosts = [g["node_id"].to_numpy() for _, g in node_map[node_map["time"] == 0].groupby("ind_id")]
    T = params["T"]
    tw = np.append(np.arange(T + 1), np.inf).astype(float)   # one window per generation, then "not coalesced by T"
    total = ts.pair_coalescence_counts(time_windows=tw)
    within = ts.pair_coalescence_counts(sample_sets=hosts, indexes=[(i, i) for i in range(len(hosts))],
                                        time_windows=tw).sum(axis=0)
    sizes = np.array([len(h) for h in hosts])
    n_same = int((sizes * (sizes - 1) // 2).sum())
    n_all = ts.num_samples * (ts.num_samples - 1) // 2
    os.makedirs(os.path.dirname(result_path(name, seed)), exist_ok=True)
    np.savez(result_path(name, seed), different=(total - within) / (n_all - n_same),
             same=within / max(n_same, 1), n_strains=ts.num_samples, n_hosts=len(hosts), **params)
    shutil.rmtree(out)
    print(f"done: {name} seed {seed}, {ts.num_samples} strains in {len(hosts)} hosts", flush=True)


def generation_at(cdf, q):
    """First generation by which a share q of the genome has coalesced (nan if not by T)."""
    i = np.searchsorted(cdf, q)
    return float(i) if i < len(cdf) - 1 else np.nan


def plot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    # sequential ramp, one hue, light -> dark = low -> high value
    ramp = ["#b7d3f5", "#7aaeea", "#2a78d6", "#1a54a3", "#0d2f63"]
    ink, muted, grid = "#0f2333", "#52514e", "#e6e5e0"
    plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2), sharey=True)
    rows = []
    for ax, (param, values) in zip(axes, SWEEP.items()):
        for value, colour in zip(values, ramp):
            params = {**BASELINE, param: value}
            name = "baseline" if params == BASELINE else f"{param}_{value}"
            cdfs = np.array([np.cumsum(np.load(result_path(name, s))["different"]) for s in SEEDS])
            g = np.arange(cdfs.shape[1] - 1)          # generations 0..T (the last window is "not by T")
            lo, mid, hi = cdfs.min(0)[:-1], np.median(cdfs, 0)[:-1], cdfs.max(0)[:-1]
            ax.fill_between(g, lo, hi, step="post", color=colour, alpha=0.25, lw=0)
            label = f"{value:g}" + (" (baseline)" if value == BASELINE[param] else "")
            ax.step(g, mid, where="post", color=colour, lw=2, label=label)
            for s, c in zip(SEEDS, cdfs):
                rows.append({"parameter": param, "value": value, "seed": s,
                             "median coalescence (generations)": generation_at(c, 0.5),
                             "coalesced by 10": c[10], "coalesced by 100": c[100], "coalesced by T": c[-2]})
        ax.set_xscale("log")
        ax.set_xlim(1, BASELINE["T"])
        ax.set_ylim(0, 1)
        ax.grid(color=grid, lw=0.8)
        ax.set_axisbelow(True)
        ax.set_title(LABEL[param], loc="left", fontweight="bold")
        ax.set_xlabel("Generations ago (log scale)")
        ax.legend(frameon=False, fontsize=10, loc="upper left", title=param, title_fontsize=10)
    axes[0].set_ylabel("Fraction of the genome coalesced\n(pairs of strains in different hosts)")
    fig.suptitle("Coalescence between strains in different hosts, varying one parameter at a time "
                 f"(others at {base_text()}). Line = median of {len(SEEDS)} seeds, "
                 "band = range", x=0.01, ha="left", fontsize=13)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(RESULTS, f"coal_sweep.{ext}"), dpi=200)
    summary = (pd.DataFrame(rows).groupby(["parameter", "value"], sort=False)
               .agg(["median", "min", "max"]).drop(columns="seed"))
    summary.columns = [f"{a} ({b})" for a, b in summary.columns]
    summary.to_csv(os.path.join(RESULTS, "coal_sweep_summary.tsv"), sep="\t")
    pd.set_option("display.width", 250)
    print(summary[[c for c in summary.columns if c.startswith(("median coalescence", "coalesced by T"))]].round(3))
    print("wrote", os.path.join(RESULTS, "coal_sweep.{png,pdf}"))
    plt.close(fig)
    plot_density()                            # 0-1,000 generations, 20-generation periods
    plot_density(bin_width=1, max_gen=100)    # 0-100 generations, every generation
    plot_density(bin_width=1, max_gen=20)     # 0-20 generations, every generation (the early rise in close-up)


def plot_density(bin_width=20, max_gen=None):
    """The distribution itself: share of the genome coalescing in each bin of bin_width generations (linear axis),
    for pairs of strains in the same host (top row) and in different hosts (bottom row)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ramp = ["#b7d3f5", "#7aaeea", "#2a78d6", "#1a54a3", "#0d2f63"]
    grid = "#e6e5e0"
    T = BASELINE["T"]
    shown = max_gen or T                               # generations shown (the runs go to T)
    edges = np.arange(0, shown + 1, bin_width)
    # one-generation periods: plot at the generation itself (coalescence at g is in [g, g + 1))
    mids = edges[:-1].astype(float) if bin_width == 1 else (edges[:-1] + edges[1:]) / 2
    rows = [("same", "Same host"), ("different", "Different hosts")]
    YMIN = 1e-4
    suffix = "" if shown == T else f"_{shown}gen"
    plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 4, figsize=(19, 9), sharex=True, sharey=True)
    for r, (key, row_label) in enumerate(rows):
        for ax, (param, values) in zip(axes[r], SWEEP.items()):
            for value, colour in zip(values, ramp):
                params = {**BASELINE, param: value}
                name = "baseline" if params == BASELINE else f"{param}_{value}"
                per_gen = np.array([np.load(result_path(name, s))[key][:T] for s in SEEDS])  # generations 0..T-1
                cdf = np.c_[np.zeros(len(SEEDS)), np.cumsum(per_gen, axis=1)]               # cdf[:, g] = by g
                binned = cdf[:, edges[1:]] - cdf[:, edges[:-1]]
                lo, mid, hi = binned.min(0), np.median(binned, 0), binned.max(0)
                # log axis: a period where nothing coalesced (0) is left as a gap
                ax.fill_between(mids, np.maximum(lo, YMIN), hi, where=hi > 0, color=colour, alpha=0.25, lw=0)
                not_by_T = 1 - np.median(cdf[:, shown])
                label = f"{value:g}" + (" (baseline)" if value == BASELINE[param] else "")
                if not_by_T > 0.005:
                    label += f"  [{not_by_T:.0%} not by {shown}]"
                ax.plot(mids, np.where(mid > 0, mid, np.nan), color=colour, lw=2, label=label,
                        marker="o" if bin_width == 1 else None, ms=4 if shown <= 20 else 2)
            ax.set_xlim(0, shown)
            if shown <= 20:
                ax.set_xticks(range(0, shown + 1, 2))
            ax.set_yscale("log")
            ax.set_ylim(YMIN, 1)
            ax.grid(color=grid, lw=0.8)
            ax.set_axisbelow(True)
            ax.legend(frameon=False, fontsize=9.5, loc="upper right", title=param, title_fontsize=10)
            if r == 0:
                ax.set_title(LABEL[param], loc="left", fontweight="bold")
            else:
                ax.set_xlabel("Coalescence time (generations ago)")
        axes[r, 0].set_ylabel(f"{row_label}\nproportion of the genome coalescing\nin each {bin_width}-generation period (log)")
    fig.suptitle("Distribution of coalescence times, pairs of strains in the same host (top) and in different hosts "
                 f"(bottom), varying one parameter at a time\n(others at {base_text()}; first {shown} generations shown). "
                 f"Line = median of {len(SEEDS)} seeds, band = range", x=0.01, ha="left", fontsize=13)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(RESULTS, f"coal_sweep_density{suffix}.{ext}"), dpi=200)
    print("wrote", os.path.join(RESULTS, f"coal_sweep_density{suffix}.{{png,pdf}}"))
    plt.close(fig)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "list":
        for i, (name, params, seed) in enumerate(runs()):
            print(i, name, seed, params)
    elif cmd == "run":
        run_one(int(sys.argv[2]))
    elif cmd == "plot":
        plot()
