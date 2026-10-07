"""
Benchmarks behind the website's "Running at scale" page (docs/scaling.qmd).

    python scaling.py pipeline --Nh 100 --T 1000     # times every krsim command for one run, appends to pipeline.csv
    python scaling.py run --Nh 100 --T 1000           # times the parts of krsim simulate, appends to benchmarks.csv
    python scaling.py coal-rows ts_mut_42.trees       # rows krsim coal-mut would write, and its estimated cost
    python scaling.py plot                            # redraws docs/images/scaling_*.png

Parameters not given are at the baseline (Nh 100, T 1000, K 10, chi 0.2, C 3).
Run one benchmark at a time so the timings don't compete for the CPU.

pipeline.csv: each command run as its own process, with its time, peak RAM and output size.
benchmarks.csv: krsim simulate across parameters, split into its parts.
coal_mut.csv: krsim coal-mut runs (rows from `coal-rows`), plus row counts for
larger runs it cannot finish.
"""

import argparse
import os
import platform
import resource
import shutil
import subprocess
import tempfile
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "benchmarks.csv")
PIPELINE_CSV = os.path.join(HERE, "pipeline.csv")
COAL_CSV = os.path.join(HERE, "coal_mut.csv")
IMAGES = os.path.join(HERE, "..", "docs", "images")
BASELINE = {"Nh": 100, "T": 1000, "K": 10, "chi": 0.2, "C": 3}

# ru_maxrss is in bytes on macOS and in kilobytes on Linux
RSS_TO_GB = (1 if platform.system() == "Darwin" else 1024) / 2**30


def peak_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * RSS_TO_GB


def machine():
    """Short description of this computer, e.g. 'Apple M5 Pro, 18 cores, 64 GB'."""
    cpu, ram_gb = platform.processor() or platform.machine(), None
    try:
        if platform.system() == "Darwin":
            cpu = subprocess.check_output(["sysctl", "-n", "machdep.cpu.brand_string"], text=True).strip()
            ram_gb = int(subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True)) / 2**30
        else:
            ram_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30
    except (OSError, ValueError):
        pass
    return f"{cpu}, {os.cpu_count()} cores" + (f", {ram_gb:.0f} GB" if ram_gb else "")


def timed(cmd):
    """Run a command; return its wall time (s) and its own peak RAM (GB)."""
    start = time.perf_counter()
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    _, status, usage = os.wait4(proc.pid, 0)
    if status != 0:
        raise RuntimeError(f"failed: {' '.join(cmd)}")
    return time.perf_counter() - start, usage.ru_maxrss * RSS_TO_GB


def pipeline(args):
    import tskit

    p = {k: getattr(args, k) for k in BASELINE}
    seed = str(args.random_seed)
    out = tempfile.mkdtemp(prefix="krsim_bench_")
    f = lambda name: os.path.join(out, name)
    steps = [
        ("simulate", f(f"ts_{seed}.trees"),
         ["krsim", "simulate", "--Nh", str(p["Nh"]), "--T", str(p["T"]), "--K", str(p["K"]),
          "--chi", str(p["chi"]), "--C", str(p["C"]), "--random_seed", seed, "--output_dir", out]),
        ("mutate", f(f"ts_mut_{seed}.trees"),
         ["krsim", "mutate", "--ts", f(f"ts_{seed}.trees"), "--random_seed", seed, "--output_dir", out]),
        ("vcf-strain", f(f"output_{seed}_strain.vcf"),
         ["krsim", "vcf-strain", "--ts", f(f"ts_mut_{seed}.trees"), "--random_seed", seed, "--output_dir", out]),
        ("vcf-host", f(f"output_{seed}_host.vcf"),
         ["krsim", "vcf-host", "--strain_vcf", f(f"output_{seed}_strain.vcf"),
          "--ind_node_map_csv", f(f"sim_ind_node_map_{seed}.csv"), "--depth", "10",
          "--random_seed", seed, "--output_dir", out]),
    ]
    rows = []
    try:
        for name, output, cmd in steps:
            secs, gb = timed(cmd)
            rows.append({**p, "command": name, "time_s": round(secs, 2), "peak_ram_GB": round(gb, 3),
                         "output_MB": round(os.path.getsize(output) / 2**20, 1)})
            print(f"{name}: {secs:.1f} s, {gb:.2f} GB, {rows[-1]['output_MB']} MB")
        ts = tskit.load(f(f"ts_mut_{seed}.trees"))
        hosts = pd.read_csv(f(f"sim_ind_node_map_{seed}.csv"))
        sizes = {"strains": ts.num_samples, "sites": ts.num_sites,
                 "hosts": hosts.loc[hosts["time"] == 0, "ind_id"].nunique()}
    finally:
        shutil.rmtree(out)  # the VCFs can be tens of GB
    df = pd.DataFrame([{**r, **sizes} for r in rows])
    df.to_csv(PIPELINE_CSV, mode="a", header=not os.path.exists(PIPELINE_CSV), index=False)


def run(args):
    import msprime
    from krsim.sim.params import SimParams
    from krsim.sim.simulator import KRSimulator

    params = SimParams(Nh=args.Nh, T=args.T, K=args.K, chi=args.chi, C=args.C,
                       recombination_rate=7.4e-7, random_seed=args.random_seed)
    sim = KRSimulator(params)
    row = {"Nh": args.Nh, "T": args.T, "K": args.K, "chi": args.chi, "C": args.C}

    # the same steps as KRSimulator.run(), timed separately
    start = time.perf_counter()
    population = sim._initialise_population()
    for t in range(params.T - 1, -1, -1):
        population = sim._run_generation(t, population)
    row["generations_s"] = time.perf_counter() - start
    row["edges_before_simplify"] = sim.tables.edges.num_rows

    start = time.perf_counter()
    sim.tables.sort()
    row["sort_s"] = time.perf_counter() - start

    start = time.perf_counter()
    sim.tables.simplify(np.where(sim.tables.nodes.time == 0)[0], keep_input_roots=True)
    row["simplify_s"] = time.perf_counter() - start
    row["simulate_s"] = row["generations_s"] + row["sort_s"] + row["simplify_s"]

    ts = sim.tables.tree_sequence()
    row.update(strains=ts.num_samples, trees=ts.num_trees)

    start = time.perf_counter()
    msprime.sim_mutations(ts, rate=1e-6, random_seed=args.random_seed)
    row["mutate_s"] = time.perf_counter() - start

    start = time.perf_counter()
    ts.pair_coalescence_counts(time_windows=np.append(np.arange(params.T + 1), np.inf).astype(float))
    row["coalescence_s"] = time.perf_counter() - start

    row["peak_ram_GB"] = peak_rss_gb()
    row["machine"] = machine()

    df = pd.DataFrame([row])
    df.to_csv(CSV, mode="a", header=not os.path.exists(CSV), index=False)
    print(df.round(2).to_string(index=False))


def coal_rows(ts):
    """
    Number of rows `krsim coal-mut` writes for a tree sequence: one per strain
    pair per tree in which the pair has a common ancestor, i.e. for each tree,
    the pairs of samples below each of its roots.
    """
    total = 0
    for tree in ts.trees():
        n = np.array([tree.num_samples(root) for root in tree.roots])
        total += int((n * (n - 1) // 2).sum())
    return total


def coal_cost_per_row():
    """Time (s), peak RAM (GB) and CSV size (MB) per row, from the measured coal-mut runs."""
    m = pd.read_csv(COAL_CSV).dropna(subset=["time_s"])
    # fit through the larger runs, where start-up time no longer matters
    m = m[m["rows"] >= 1e7]
    return tuple((m[col] / m["rows"]).mean() for col in ["time_s", "peak_ram_GB", "csv_MB"])


def coal_rows_cmd(args):
    import tskit
    per_s, per_gb, per_mb = coal_cost_per_row()
    for path in args.trees:
        ts = tskit.load(path)
        rows = coal_rows(ts)
        print(f"{path}: {ts.num_samples:,} strains, {ts.num_trees:,} trees -> {rows:,} rows; "
              f"estimated {rows * per_s / 60:,.0f} min, {rows * per_gb:,.1f} GB RAM, "
              f"{rows * per_mb / 1024:,.1f} GB CSV (on the machine in benchmarks.csv)")


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
    plt.rcParams.update({
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.8,
        "lines.linewidth": 2, "lines.markersize": 6, "axes.axisbelow": True,
    })
    where = pd.read_csv(CSV)["machine"].iloc[0]

    # --- every command against Nh: time, peak RAM, output size ---
    pipe = pd.read_csv(PIPELINE_CSV)
    pipe = pipe[(pipe["T"] == BASELINE["T"]) & (pipe["K"] == BASELINE["K"])].sort_values("Nh")
    commands = [("simulate", BLUE), ("mutate", ORANGE), ("vcf-strain", AQUA), ("vcf-host", YELLOW)]
    # one column, three rows, so each panel is large on the website
    fig, axes = plt.subplots(3, 1, figsize=(9, 13), sharex=True)
    for ax, col, label, scale in [(axes[0], "time_s", "Time (minutes)", 60),
                                  (axes[1], "peak_ram_GB", "Peak RAM (GB)", 1),
                                  (axes[2], "output_MB", "Output size (GB)", 1024)]:
        ends = []
        for name, colour in commands:
            d = pipe[pipe["command"] == name]
            ax.plot(d["Nh"], d[col] / scale, "o-", color=colour, label=name)
            ends.append([d[col].iloc[-1] / scale, name])
        # label each line at its end, nudging labels apart (log axis) so they don't overlap
        ends.sort()
        for i in range(1, len(ends)):
            ends[i][0] = max(ends[i][0], ends[i - 1][0] * 1.35)
        for y, name in ends:
            ax.annotate(name, (pipe["Nh"].max(), y), xytext=(8, 0), textcoords="offset points",
                        va="center", color="#52514e", fontsize=9, annotation_clip=False)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(pipe["Nh"].unique())
        ax.set_xticklabels(pipe["Nh"].unique())
        ax.minorticks_off()
        ax.set_xlim(None, pipe["Nh"].max() * 2.2)
        ax.set_title(label, loc="left")
    axes[-1].set_xlabel("Nh (number of hosts)")
    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle(f"Every krsim command against the number of hosts (T = {BASELINE['T']}, K = {BASELINE['K']})\n"
                 f"{where}", x=0.01, ha="left")
    plt.tight_layout()
    plt.savefig(os.path.join(IMAGES, "scaling_pipeline.png"), dpi=110)
    plt.close()

    # --- coal-mut: measured runs and extrapolation by number of rows ---
    coal = pd.read_csv(COAL_CSV)
    measured = coal.dropna(subset=["time_s"])
    predicted = coal[coal["time_s"].isna()]
    per_s, per_gb, _ = coal_cost_per_row()
    line = np.geomspace(coal["rows"].min(), coal["rows"].max() * 1.5, 50)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, col, per, scale in [(axes[0], "time_s", per_s, 3600), (axes[1], "peak_ram_GB", per_gb, 1)]:
        ax.plot(line, line * per / scale, color=BLUE, linewidth=1.5, linestyle="--", label="extrapolated")
        ax.plot(measured["rows"], measured[col] / scale, "o", color=BLUE, label="measured")
        ax.plot(predicted["rows"], predicted["rows"] * per / scale, "o", color=BLUE,
                markerfacecolor="white", markeredgewidth=2)
        for _, r in coal.iterrows():
            if r["T"] == 1000:
                y = r[col] / scale if pd.notna(r[col]) else r["rows"] * per / scale
                ax.annotate(f"Nh = {int(r['Nh'])}", (r["rows"], y), xytext=(8, -10), textcoords="offset points",
                            ha="left", color="#52514e", fontsize=9)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Rows in coal_T table")
    for gb, label in [(16, "16 GB laptop"), (64, "64 GB")]:
        axes[1].axhline(gb, color=ORANGE, linewidth=1)
        axes[1].annotate(label, (line[-1], gb), xytext=(0, 3), textcoords="offset points",
                         ha="right", color=ORANGE, fontsize=9)
    axes[0].set_ylabel("Time (hours)")
    axes[1].set_ylabel("Peak RAM (GB)")
    axes[0].set_title("krsim coal-mut run time", loc="left")
    axes[1].set_title("krsim coal-mut peak RAM", loc="left")
    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle(f"krsim coal-mut cost grows with the number of rows (K = 10, T = 1000) — {where}",
                 x=0.01, ha="left")
    plt.tight_layout()
    plt.savefig(os.path.join(IMAGES, "scaling_coal_mut.png"), dpi=110)
    plt.close()
    print("wrote docs/images/scaling_pipeline.png, docs/images/scaling_coal_mut.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, help_text in [("pipeline", "time every krsim command for one parameter set (appends to pipeline.csv)"),
                            ("run", "time the parts of krsim simulate for one parameter set (appends to benchmarks.csv)")]:
        s = sub.add_parser(name, help=help_text)
        for param, typ in [("Nh", int), ("T", int), ("K", int), ("chi", float), ("C", int)]:
            s.add_argument(f"--{param}", type=typ, default=BASELINE[param])
        s.add_argument("--random_seed", type=int, default=42)
    c = sub.add_parser("coal-rows", help="rows krsim coal-mut would write for these .trees files, and its estimated cost")
    c.add_argument("trees", nargs="+")
    sub.add_parser("plot", help="draw the figures from the CSVs")
    args = parser.parse_args()
    {"pipeline": pipeline, "run": run, "coal-rows": coal_rows_cmd, "plot": plot}[args.command](args)
