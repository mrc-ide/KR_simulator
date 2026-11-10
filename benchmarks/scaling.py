"""
Benchmarks behind the website's "Running at scale" page (docs/scaling.qmd).

    python scaling.py pipeline --Nh 100 --T 1000     # times every krsim command for one run, appends to pipeline.csv
    python scaling.py run --Nh 100 --T 1000           # times the parts of krsim simulate, appends to benchmarks.csv
    python scaling.py plot                            # redraws docs/images/scaling_pipeline.png

Parameters not given are at the baseline (Nh 100, T 1000, K 10, chi 0.2, C 3).
Run one benchmark at a time so the timings don't compete for the CPU.

pipeline.csv: each command run as its own process, with its time, peak RAM and output size.
benchmarks.csv: krsim simulate across parameters, split into its parts.
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
IMAGES = os.path.join(HERE, "..", "docs", "images")
BASELINE = {"Nh": 100, "T": 1000, "K": 10, "chi": 0.2, "C": 3}

# typical laptop / cluster core is about this many times slower than the benchmark laptop
TYPICAL_SLOWDOWN = 2

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
    out = tempfile.mkdtemp(prefix="krsim_bench_")  # set TMPDIR to a disk with room: the VCFs reach tens of GB
    f = lambda name: os.path.join(out, name)
    steps = [
        ("simulate", f(f"ts_{seed}.trees"),
         ["krsim", "simulate", "--Nh", str(p["Nh"]), "--T", str(p["T"]), "--K", str(p["K"]),
          "--chi", str(p["chi"]), "--C", str(p["C"]), "--random_seed", seed, "--output_dir", out]),
        ("mutate", f(f"ts_mut_{seed}.trees"),
         ["krsim", "mutate", "--ts", f(f"ts_{seed}.trees"), "--random_seed", seed, "--output_dir", out]),
        ("vcf-strain", f(f"output_{seed}_strain.vcf"),
         ["krsim", "vcf-strain", "--ts", f(f"ts_mut_{seed}.trees"), "--random_seed", seed, "--output_dir", out]),
        ("vcf-host-depth", f(f"output_{seed}_host_depth.vcf"),
         ["krsim", "vcf-host-depth", "--strain_vcf", f(f"output_{seed}_strain.vcf"),
          "--ind_node_map_csv", f(f"sim_ind_node_map_{seed}.csv"), "--depth_model", "sWGA",
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
    # replace any earlier measurements of this parameter set rather than adding duplicates
    csv_path = args.csv or PIPELINE_CSV
    if os.path.exists(csv_path):
        old = pd.read_csv(csv_path)
        same = np.ones(len(old), dtype=bool)
        for k, v in p.items():
            same &= np.isclose(old[k].astype(float), float(v))
        df = pd.concat([old[~same], df], ignore_index=True)
    df.to_csv(csv_path, index=False)


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


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker

    BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
    plt.rcParams.update({
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.8,
        "lines.linewidth": 2.5, "lines.markersize": 8, "axes.axisbelow": True,
        "font.size": 15, "axes.titlesize": 19, "axes.titleweight": "bold",
        "axes.labelsize": 15, "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14,
    })
    # The benchmarks ran on a fast laptop; times are shown for a typical laptop
    # or cluster core, estimated as TYPICAL_SLOWDOWN x the measured time.
    # Memory and disk do not depend on the machine.
    where = pd.read_csv(CSV)["machine"].iloc[0]

    # --- every command against Nh: time, peak RAM, output size ---
    pipe = pd.read_csv(PIPELINE_CSV)
    pipe = pipe[(pipe["T"] == BASELINE["T"]) & (pipe["K"] == BASELINE["K"])].sort_values("Nh")
    commands = [("simulate", BLUE), ("mutate", ORANGE), ("vcf-strain", AQUA), ("vcf-host-depth", YELLOW)]
    # one column, three rows, so each panel is large on the website
    fig, axes = plt.subplots(3, 1, figsize=(10, 17))
    for ax, col, label, scale in [(axes[0], "time_s", "Time (minutes, typical laptop)", 60 / TYPICAL_SLOWDOWN),
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
                        va="center", color="#52514e", fontsize=14, annotation_clip=False)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # plain numbers (0.1, 1, 10) rather than powers of ten
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_xticks(pipe["Nh"].unique())
        ax.set_xticklabels(pipe["Nh"].unique())
        ax.minorticks_off()
        ax.set_xlim(None, pipe["Nh"].max() * 2.2)
        ax.set_title(label, loc="left", pad=12)
        ax.set_xlabel("Nh (number of hosts)")
    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle(f"Every krsim command against the number of hosts (T = {BASELINE['T']}, K = {BASELINE['K']})",
                 x=0.01, ha="left", fontsize=18)
    plt.tight_layout(h_pad=4)
    plt.savefig(os.path.join(IMAGES, "scaling_pipeline.png"), dpi=110)
    plt.close()

    print("wrote docs/images/scaling_pipeline.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, help_text in [("pipeline", "time every krsim command for one parameter set (appends to pipeline.csv)"),
                            ("run", "time the parts of krsim simulate for one parameter set (appends to benchmarks.csv)")]:
        s = sub.add_parser(name, help=help_text)
        for param, typ in [("Nh", int), ("T", int), ("K", int), ("chi", float), ("C", int)]:
            s.add_argument(f"--{param}", type=typ, default=BASELINE[param])
        s.add_argument("--random_seed", type=int, default=42)
        if name == "pipeline":
            s.add_argument("--csv", help="write to this CSV instead of pipeline.csv (e.g. for another machine)")
    sub.add_parser("plot", help="draw the figures from the CSVs")
    args = parser.parse_args()
    {"pipeline": pipeline, "run": run, "plot": plot}[args.command](args)
