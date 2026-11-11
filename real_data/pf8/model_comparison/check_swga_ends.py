"""
Evidence for the sWGA chromosome-end zone in krsim's depth model (end_width 30 kb, end_factor 0.0035).

Real: N random QC-pass Pf8 sWGA samples (window stats in $PF8/output/). For each sample, every 500 bp window is
divided by the sample's median; for each distance from a chromosome end (both ends of all 14 nuclear chromosomes)
the sample's median over those 28 ends gives its end profile, 0-200 kb.

  1. Where the drop is: the real median profile in 5 kb bins (as make_depth_model_defaults.py uses), the zone
     width each threshold would give (the run of bins from the end below it), and the steepest rise.
  2. The factor: the median depth of the zone, and the share of windows with no reads at all (real) against
     the share of sites with 0 reads that the factor gives in krsim (Poisson at the host's typical depth).
  3. Options compared on simulated hosts (krsim's sWGA model, same random tracks for every option): zone widths
     0-100 kb, each with its factor = the real median over that width, plus an end profile (one factor per 5 kb
     bin out to 150 kb). Score: mean |log10(simulated / real)| of the median profile per 500 bp window, in
     0-30, 30-100 and 0-200 kb from the end (0 = identical; 1 = 10x off on average).
  4. The whole genome: AR(8) alone against AR(8) with the 30 kb zone or the end profile. Per real sample and per
     simulated host, over every window: the share of windows below 0.01x, 0.1x and 0.5x the median, and the
     1st, 5th and 10th percentile window. Compared by the median over samples / hosts and the KS distance
     between real and simulated (0 = same distribution).

    python check_swga_ends.py [N]     # default 300; writes results/swga_ends/
"""

import os
import sys

import numpy as np
import pandas as pd

from krsim.output import depth_model as dm
from krsim.output.pf8_profiles import CHROMS, WINDOW_SIZE

PF8 = "/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "swga_ends")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
N_SIM = 200
SEED = 2026
W = 400            # windows from each end: 200 kb
BIN = 10           # windows per 5 kb bin
FLOOR = 1e-4       # depth floor for logs (a median of 0 reads)


def end_profile(chroms):
    """One genome's depth (relative to its median) by distance from the end: median over every chromosome end."""
    return np.median([e for c in chroms for e in (c[:W], c[::-1][:W])], axis=0)


def genome_stats(rel):
    """Whole-genome summary of one genome's relative depth (every window / the median)."""
    return {"share <0.01x": np.mean(rel < 0.01), "share <0.1x": np.mean(rel < 0.1), "share <0.5x": np.mean(rel < 0.5),
            "p1": np.percentile(rel, 1), "p5": np.percentile(rel, 5), "p10": np.percentile(rel, 10)}


def real_profiles(samples):
    profiles, zero_share, stats = [], [], []
    for s in samples:
        df = pd.read_csv(f"{PF8}/output/{s}.window_stats.tsv", sep="\t", usecols=["chrom", "mean"])
        df = df[~df["chrom"].str.contains("API|MIT")]
        m = df["mean"].median()
        if not m > 0:
            continue
        chroms = [g.to_numpy() / m for _, g in df.groupby("chrom", sort=False)["mean"]]
        profiles.append(end_profile(chroms))
        stats.append(genome_stats(np.concatenate(chroms)))
        ends = np.concatenate([np.r_[c[:60], c[::-1][:60]] for c in chroms])  # 0-30 kb
        zero_share.append(np.mean(ends == 0))
    return np.array(profiles), np.array(zero_share), pd.DataFrame(stats)


def bins(profile):
    return profile.reshape(-1, BIN).mean(axis=1)


def score(sim, real, a, b):
    s, r = np.maximum(sim[a:b], FLOOR), np.maximum(real[a:b], FLOOR)
    return np.mean(np.abs(np.log10(s / r)))


def main():
    os.makedirs(OUT, exist_ok=True)
    table = pd.read_csv(f"{PF8}/model_fits/depth_model_samples.tsv", sep="\t")
    rng = np.random.default_rng(SEED)
    picked = rng.choice(table.loc[table["sample type"] == "sWGA", "sample"].to_numpy(), N, replace=False)
    real, zero_share, real_stats = real_profiles(picked)
    real_med = np.median(real, axis=0)
    out = []

    # 1. where the drop is
    real_bins = np.array([np.median(np.median(real[:, i * BIN:(i + 1) * BIN], axis=1)) for i in range(W // BIN)])
    out.append(f"sWGA, {len(real)} real samples\n\n1. Real depth by distance from the chromosome end (median sample, 5 kb bins):")
    out += [f"   {i * 5:>3}-{i * 5 + 5:<3} kb  {v:.4f}x" for i, v in enumerate(real_bins[:24])]
    out.append("\n   Zone width each threshold gives (bins from the end below it):")
    for thr in (0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8):
        n = int(np.argmax(real_bins >= thr)) if (real_bins >= thr).any() else len(real_bins)
        out.append(f"   below {thr:<4}x -> {n * 5:>3} kb")
    rise = np.diff(np.log10(np.maximum(real_med, FLOOR)))
    k = int(np.argmax(rise))
    out.append(f"\n   Steepest rise (500 bp windows): {k * 0.5:.1f}-{(k + 1) * 0.5:.1f} kb from the end, x{10 ** rise[k]:.0f} "
               f"({real_med[k]:.4f}x -> {real_med[k + 1]:.3f}x)")
    for a, b in ((0, 25), (25, 30), (30, 50), (50, 100), (100, 150), (150, 200)):
        out.append(f"   {a:>3}-{b:<3} kb: median {np.median(real_med[a * 2:b * 2]):.4f}x")

    # 2. the factor
    s = dm.model_settings("sWGA")
    z = int(s["end_width"] // WINDOW_SIZE)
    out.append(f"\n2. Factor: krsim uses {s['end_factor']:.4f}x over {s['end_width'] / 1000:.0f} kb "
               f"(median of the 5 kb bin medians, all 17,330 samples)")
    out.append(f"   This sample set: median of 0-30 kb windows {np.median(real_med[:z]):.4f}x")
    out.append(f"   Real: share of 0-30 kb windows with no reads at all, median sample {np.median(zero_share):.0%} "
               f"(middle 50% of samples {np.percentile(zero_share, 25):.0%}-{np.percentile(zero_share, 75):.0%})")
    for depth in (14, 30, 63, 132, 300):
        out.append(f"   krsim, host typical depth {depth:>3}x: expected reads per site {depth * s['end_factor']:.2f}, "
                   f"sites with 0 reads {np.exp(-depth * s['end_factor']):.0%}")

    # 3. options on simulated hosts (the same tracks for every option)
    base = dict(s, end_width=0.0, end_factor=1.0, cap=np.inf)
    raw = [dm.relative_track(np.random.default_rng(SEED + h), base) for h in range(N_SIM)]
    # raw tracks before the median scaling can't be recovered, so each option rescales: track x end factors, then
    # divided by its own genome median (as relative_track does), then capped
    options = {}
    for kb in (0, 10, 20, 25, 30, 35, 40, 50, 75, 100):
        f = np.ones(W)
        if kb:
            f[:kb * 2] = np.median(real_med[:kb * 2])
        options[f"flat zone {kb} kb" + (" (krsim now)" if kb == 30 else "")] = f
    prof = np.ones(W)
    prof[:300] = np.repeat(np.maximum(real_bins[:30], FLOOR), BIN)  # one factor per 5 kb bin, 0-150 kb
    options["end profile (5 kb bins to 150 kb)"] = prof

    rows, sim_meds, sim_stats = [], {}, {}
    for name, f in options.items():
        host_profiles, host_stats = [], []
        for tr in raw:
            t = {c: v.copy() for c, v in tr.items()}
            for c in CHROMS:
                n = min(W, len(t[c]))
                t[c][:n] *= f[:n]
                t[c][len(t[c]) - n:] *= f[:n][::-1]
            med = np.median(np.concatenate(list(t.values())))
            capped = [np.minimum(v / med, s["cap"]) for v in t.values()]
            host_profiles.append(end_profile(capped))
            host_stats.append(genome_stats(np.concatenate(capped)))
        sim_stats[name] = pd.DataFrame(host_stats)
        sim = np.median(host_profiles, axis=0)
        sim_meds[name] = sim
        rows.append({"option": name, "0-30 kb": score(sim, real_med, 0, 60), "30-100 kb": score(sim, real_med, 60, 200),
                     "0-200 kb": score(sim, real_med, 0, W), "factor": f[0]})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "options.tsv"), sep="\t", index=False)
    out.append("\n3. Mean |log10(simulated / real)| of the median end profile (0 = identical, 1 = 10x off):")
    out.append(res.round(3).to_string(index=False))

    # 4. the whole genome: AR(8) alone vs with the ends
    from scipy.stats import ks_2samp
    compare = {"AR(8) alone": "flat zone 0 kb", "AR(8) + 30 kb zone (krsim now)": "flat zone 30 kb (krsim now)",
               "AR(8) + end profile": "end profile (5 kb bins to 150 kb)"}
    rows4 = [{"genome": "REAL", **{f"{c} median": real_stats[c].median() for c in real_stats}}]
    for label, key in compare.items():
        d = sim_stats[key]
        row = {"genome": label, **{f"{c} median": d[c].median() for c in d}}
        row.update({f"{c} KS": ks_2samp(d[c], real_stats[c]).statistic for c in d})
        row["mean KS"] = np.mean([row[f"{c} KS"] for c in d])
        rows4.append(row)
    res4 = pd.DataFrame(rows4).set_index("genome")
    res4.to_csv(os.path.join(OUT, "whole_genome.tsv"), sep="\t")
    out.append("\n4. Whole genome, every window (median over real samples / simulated hosts):")
    out.append(res4[[c for c in res4.columns if c.endswith("median")]].round(4).to_string())
    out.append("\n   KS distance to real (0 = same distribution over samples):")
    out.append(res4.drop(index="REAL")[[c for c in res4.columns if c.endswith("KS")]].round(3).to_string())

    np.savez(os.path.join(OUT, "profiles.npz"), real=real, **{k: v for k, v in sim_meds.items()})
    text = "\n".join(out)
    open(os.path.join(OUT, "summary.txt"), "w").write(text + "\n")
    print(text)
    plot(real, sim_meds)


def plot(real, sim_meds):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = (np.arange(W) + 0.5) * 0.5
    plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(12, 6.5))
    q25, q50, q75 = np.percentile(real, [25, 50, 75], axis=0)
    ax.fill_between(x, np.maximum(q25, FLOOR), np.maximum(q75, FLOOR), color="#0f2333", alpha=0.12, lw=0,
                    label="real sWGA: middle 50% of samples")
    ax.plot(x, np.maximum(q50, FLOOR), color="#0f2333", lw=3, label="real sWGA: median sample")
    show = {"flat zone 0 kb": ("#eb6834", "model, no end zone"),
            "flat zone 30 kb (krsim now)": ("#2a78d6", "model, 30 kb zone at 0.0035x (krsim now)"),
            "end profile (5 kb bins to 150 kb)": ("#1baf7a", "model, end profile (one factor per 5 kb)")}
    for k, (c, lab) in show.items():
        ax.plot(x, np.maximum(sim_meds[k], FLOOR), color=c, lw=2, label=lab)
    ax.axvline(30, color="#9a9890", lw=1, ls=":")
    ax.set_yscale("log")
    ax.set_xlim(0, 200)
    ax.set_yticks([0.001, 0.01, 0.1, 1])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}×"))
    ax.grid(color="#e6e5e0", lw=0.8)
    ax.set_axisbelow(True)
    ax.set_xlabel("kb from the chromosome end (both ends of all 14 chromosomes)")
    ax.set_ylabel("Depth ÷ sample's median (log scale)")
    ax.legend(frameon=False, loc="lower right", fontsize=11)
    ax.set_title("sWGA depth near chromosome ends: real Pf8 samples vs the depth model", loc="left", fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "swga_ends.png"), dpi=110)


if __name__ == "__main__":
    main()
