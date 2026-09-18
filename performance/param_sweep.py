"""
Time every stage of the pipeline (as run by scripts/test_sim_mut.sh) across a
one-parameter-at-a-time sweep around the standard test set
(--Nh 3 --T 10 --K 2 --chi 0.18 --C 1), and plot how each stage's runtime
scales with each parameter.

Usage (from repo root):
    python scripts/benchmark_params.py
    python scripts/benchmark_params.py --output_dir scripts/benchmark_output --seed 42
"""

import os
import sys
import time
import argparse
import subprocess

import pandas as pd
import matplotlib.pyplot as plt

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BASE = {"Nh": 3, "T": 10, "K": 2, "chi": 0.18, "C": 1}
RECOMBINATION_RATE = "7.4e-7"
MUTATION_RATE = "1e-6"
DEPTH = "10"
TIME_BOUNDARY = "2"

# one-at-a-time sweeps, holding the other parameters at BASE
SWEEPS = {
    "Nh":  [3, 10, 30, 60, 100, 200],
    "T":   [10, 100, 300, 600, 1000],
    "K":   [1, 2, 4, 8, 16],
    "chi": [0.05, 0.18, 0.4, 0.7, 0.95],
    "C":   [1, 2, 4, 8, 16],
}

REPLICATES = 3

STAGE_ORDER = [
    "1_sim", "2_transmission_map", "3_pairwise_coalescence", "4_mutations",
    "5_coal_mut_t", "6_vcf_strain", "7_vcf_host", "8_relate_strain", "9_relate_host",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output_dir", type=str,
                        default=os.path.expanduser("~/IMPERIAL/MRC_IDE/kr_sim_runs/benchmark_output"),
                        help="Directory for the timing CSV, plots, and per-config run artifacts. "
                             "Kept outside the repo since output here is machine-specific and can run into GBs.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed used for every run (only speed is measured).")
    parser.add_argument("--hmmibd", type=str, default=os.environ.get("HMMIBD", "/data/home/bty208/hmmIBD_org/hmmIBD/hmmibd"),
                        help="Path to the hmmIBD binary; the benchmark stage is skipped if this isn't runnable.")
    parser.add_argument("--replicates", type=int, default=REPLICATES,
                        help="Repeats per unique parameter combination, plotted as mean with min/max error bars.")
    return parser.parse_args()


def _cfg_key(cfg):
    return tuple(sorted(cfg.items()))


def _run(cmd, cwd):
    start = time.perf_counter()
    subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)
    return time.perf_counter() - start


def run_pipeline(cfg, workdir, seed, hmmibd_bin):
    """Run one full pipeline (Stage 1-5 + relatedness extras) for this
    parameter combination and return {stage_name: elapsed_seconds}."""
    os.makedirs(workdir, exist_ok=True)
    py = sys.executable
    timings = {}

    ts_file = os.path.join(workdir, f"ts_{seed}.trees")
    mut_ts_file = os.path.join(workdir, f"mutated_{seed}.trees")
    strain_vcf = os.path.join(workdir, f"output_{seed}_strain.vcf")
    relate_strain_csv = os.path.join(workdir, f"relate_strain_{seed}.csv")

    timings["1_sim"] = _run([
        py, "run_sim.py",
        "--Nh", str(cfg["Nh"]), "--T", str(cfg["T"]), "--K", str(cfg["K"]),
        "--chi", str(cfg["chi"]), "--C", str(cfg["C"]),
        "--recombination_rate", RECOMBINATION_RATE,
        "--random_seed", str(seed), "--out_dir", workdir,
        "--save_ts", "--save_ID", "--save_pairs", "--save_lineage",
    ], REPO_ROOT)

    timings["2_transmission_map"] = _run([
        py, "visualization/run_transmission_map.py",
        "--random_seed", str(seed), "--output_dir", workdir,
    ], REPO_ROOT)

    timings["3_pairwise_coalescence"] = _run([
        py, "visualization/run_pairwise_coalescence.py",
        "--random_seed", str(seed), "--output_dir", workdir,
    ], REPO_ROOT)

    timings["4_mutations"] = _run([
        py, "run_mut.py",
        "--ts", ts_file, "--mutation_rate", MUTATION_RATE,
        "--random_seed", str(seed), "--output_dir", workdir,
    ], REPO_ROOT)

    timings["5_coal_mut_t"] = _run([
        py, "run_coal_mut_t.py",
        "--ts", mut_ts_file, "--random_seed", str(seed), "--output_dir", workdir,
    ], REPO_ROOT)

    timings["6_vcf_strain"] = _run([
        py, "run_vcf_strain.py",
        "--ts", mut_ts_file, "--random_seed", str(seed), "--output_dir", workdir,
    ], REPO_ROOT)

    timings["7_vcf_host"] = _run([
        py, "run_vcf_host.py",
        "--strain_vcf", strain_vcf, "--random_seed", str(seed),
        "--depth", DEPTH, "--output_dir", workdir, "--save_proportions",
    ], REPO_ROOT)

    timings["8_relate_strain"] = _run([
        py, "extras/relate_t_strain.py",
        "--random_seed", str(seed), "--time_boundary", TIME_BOUNDARY, "--output_dir", workdir,
    ], REPO_ROOT)

    timings["9_relate_host"] = _run([
        py, "extras/relate_t_host.py",
        "--random_seed", str(seed), "--output_dir", workdir,
    ], REPO_ROOT)

    if hmmibd_bin and (os.path.isfile(hmmibd_bin) and os.access(hmmibd_bin, os.X_OK)):
        timings["10_hmmibd_benchmark"] = _run([
            py, "extras/benchmark_hmmibd.py",
            "--strain_vcf", strain_vcf, "--relate_strain", relate_strain_csv,
            "--random_seed", str(seed), "--output_dir", workdir, "--hmmibd", hmmibd_bin,
        ], REPO_ROOT)

    return timings


def collect_timings(output_dir, seed, hmmibd_bin, replicates):
    runs_dir = os.path.join(output_dir, "runs")
    unique_configs = {}
    for pname, values in SWEEPS.items():
        for v in values:
            cfg = dict(BASE)
            cfg[pname] = v
            unique_configs.setdefault(_cfg_key(cfg), cfg)

    total = len(unique_configs) * replicates
    done = 0
    results = {}  # key -> list of {stage: seconds} dicts, one per replicate
    for key, cfg in unique_configs.items():
        label = "_".join(f"{k}{v}" for k, v in key)
        results[key] = []
        for rep in range(replicates):
            done += 1
            print(f"[{done}/{total}] rep {rep + 1}/{replicates}  {cfg}")
            try:
                results[key].append(run_pipeline(cfg, os.path.join(runs_dir, f"{label}_r{rep}"), seed, hmmibd_bin))
            except subprocess.CalledProcessError as e:
                print(f"  FAILED: {e.cmd}\n  stderr: {e.stderr[-2000:]}")
                results[key].append({})

    rows = []
    for pname, values in SWEEPS.items():
        for v in values:
            cfg = dict(BASE)
            cfg[pname] = v
            for rep, timings in enumerate(results.get(_cfg_key(cfg), [])):
                for stage, seconds in timings.items():
                    rows.append({"sweep": pname, "value": v, "stage": stage, "replicate": rep, "seconds": seconds})
    return pd.DataFrame(rows)


# a fixed, ordered, colorblind-distinguishable qualitative palette (tab10),
# assigned to stages by position so a stage always gets the same color/marker
STAGE_COLORS = plt.get_cmap("tab10").colors
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<"]


def _agg(df, group_cols):
    """mean/min/max of `seconds` across replicates, plus asymmetric yerr."""
    g = df.groupby(group_cols, as_index=False)["seconds"].agg(["mean", "min", "max"])
    g["yerr_lo"] = g["mean"] - g["min"]
    g["yerr_hi"] = g["max"] - g["mean"]
    return g


def plot_sweep(df, pname, output_dir):
    sub = df[df["sweep"] == pname]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, stage in enumerate(STAGE_ORDER):
        stage_df = _agg(sub[sub["stage"] == stage], ["value"]).sort_values("value")
        if stage_df.empty:
            continue
        ax.errorbar(stage_df["value"], stage_df["mean"],
                    yerr=[stage_df["yerr_lo"], stage_df["yerr_hi"]],
                    marker=MARKERS[i % len(MARKERS)], color=STAGE_COLORS[i % len(STAGE_COLORS)],
                    label=stage, linewidth=1.6, markersize=5, capsize=3, elinewidth=1)
    ax.set_xlabel(pname)
    ax.set_ylabel("wall time (s), mean ± min/max over replicates")
    ax.set_yscale("log")
    held = ", ".join(f"{k}={v}" for k, v in BASE.items() if k != pname)
    ax.set_title(f"Runtime vs {pname}  (others held at {held}, n={REPLICATES})")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.grid(True, which="both", axis="y", alpha=0.25)
    fig.tight_layout()
    out_path = os.path.join(output_dir, f"timing_vs_{pname}.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def plot_overview(df, output_dir):
    fig, axes = plt.subplots(1, len(SWEEPS), figsize=(4 * len(SWEEPS), 4), sharey=True)
    per_replicate_totals = df.groupby(["sweep", "value", "replicate"], as_index=False)["seconds"].sum()
    totals = _agg(per_replicate_totals, ["sweep", "value"])
    for ax, pname in zip(axes, SWEEPS):
        sub = totals[totals["sweep"] == pname].sort_values("value")
        ax.errorbar(sub["value"], sub["mean"], yerr=[sub["yerr_lo"], sub["yerr_hi"]],
                    marker="o", color=STAGE_COLORS[0], linewidth=1.8, capsize=3, elinewidth=1)
        ax.set_xlabel(pname)
        ax.set_title(pname)
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("total pipeline wall time (s), mean ± min/max")
    fig.suptitle(f"Total pipeline runtime (sum of all stages) vs each parameter (n={REPLICATES})")
    fig.tight_layout()
    out_path = os.path.join(output_dir, "timing_overview_total.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    global REPLICATES
    args = parse_args()
    REPLICATES = args.replicates
    os.makedirs(args.output_dir, exist_ok=True)

    df = collect_timings(args.output_dir, args.seed, args.hmmibd, args.replicates)
    csv_path = os.path.join(args.output_dir, "timings.csv")
    df.to_csv(csv_path, index=False)
    print(f"wrote {csv_path}")

    for pname in SWEEPS:
        plot_sweep(df, pname, args.output_dir)
    plot_overview(df, args.output_dir)


if __name__ == "__main__":
    main()
