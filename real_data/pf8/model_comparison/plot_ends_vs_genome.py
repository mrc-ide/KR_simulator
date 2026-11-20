"""
Figure: read depth near chromosome ends against the rest of the genome, sWGA and gDNA (supplementary figure).

N random QC-pass Pf8 samples per type (window stats in $PF8/output/; nuclear chromosomes, 500 bp windows).
Each window is divided by its sample's median. Per sample, for each region, the median over both ends of all 14
chromosomes: 0-25, 25-50, 50-100 and 100-200 kb from the end, and the rest of the genome (more than 200 kb from
either end). One dot per sample; box = middle 50% of samples, whiskers = 1.5 IQR; grey band = the median
sample's 10-90% range of windows in the rest of the genome.

Shows that sWGA ends are near empty (the reason for krsim's sWGA end zone), while gDNA ends stay within the
usual range (no zone). The hatched area marks the 0-50 kb left out of the model fits (TRIM_KB in
compare_processes.py): in sWGA, depth there is below the usual range; beyond 50 kb it is within it.

Second figure (ends_profile): the same samples as a continuous profile, every 500 bp window out to 200 kb from
the end. Per sample, the median over both ends of all 14 chromosomes at each distance; line = median sample,
band = middle 50% of samples. Dashed lines: krsim's sWGA end zone (30 kb) and the fitting trim (50 kb).

    python plot_ends_vs_genome.py [N]     # default 300; writes results/swga_ends/ends_vs_genome.{png,pdf,tsv}
                                          # and ends_profile.{png,pdf}
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PF8 = "/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "swga_ends")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
SEED = 2026
ZONES = [("0–25 kb", 0, 50), ("25–50 kb", 50, 100), ("50–100 kb", 100, 200), ("100–200 kb", 200, 400)]  # 500 bp windows
TRIMMED = 2                                       # the first 2 regions (0-50 kb) are left out of the model fits
REST_FROM = 400                                                                   # 200 kb from either end
END_ZONE_KB, TRIM_KB = 30, 50                     # krsim's sWGA end zone; trimmed from each end for fitting
COLOUR = {"sWGA": "#2a78d6", "gDNA": "#eb6834"}
INK, MUTED, BAND, GRID = "#0f2333", "#52514e", "#9a9890", "#e6e5e0"


def sample_regions(sample):
    df = pd.read_csv(f"{PF8}/output/{sample}.window_stats.tsv", sep="\t", usecols=["chrom", "mean"])
    df = df[~df["chrom"].str.contains("API|MIT")]
    m = df["mean"].median()
    if not m > 0:
        return None
    ch = [g.to_numpy() / m for _, g in df.groupby("chrom", sort=False)["mean"]]
    row = {"sample": sample}
    for name, a, b in ZONES:
        row[name] = np.median(np.concatenate([np.r_[c[a:b], c[::-1][a:b]] for c in ch]))
    rest = np.concatenate([c[REST_FROM:len(c) - REST_FROM] for c in ch])
    row["rest of genome"] = np.median(rest)
    row["rest p10"], row["rest p90"] = np.percentile(rest, [10, 90])
    # per-window profile out to 200 kb: median over both ends of every chromosome
    row["profile"] = np.median([e for c in ch for e in (c[:REST_FROM], c[::-1][:REST_FROM])], axis=0)
    return row


def main():
    os.makedirs(OUT, exist_ok=True)
    table = pd.read_csv(f"{PF8}/model_fits/depth_model_samples.tsv", sep="\t")
    rng = np.random.default_rng(SEED)
    res = {}
    for t in ("sWGA", "gDNA"):
        picked = rng.choice(table.loc[table["sample type"] == t, "sample"].to_numpy(), N, replace=False)
        res[t] = pd.DataFrame([r for r in map(sample_regions, picked) if r is not None]).assign(type=t)
        print(t, len(res[t]), "samples;", res[t].drop(columns=["sample", "type", "profile"]).median().round(4).to_dict())
    profiles = {t: np.vstack(res[t].pop("profile").to_numpy()) for t in res}
    pd.concat(res.values()).to_csv(os.path.join(OUT, "ends_vs_genome.tsv"), sep="\t", index=False)

    cats = [z[0] for z in ZONES] + ["rest of genome"]
    plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.4), sharey=True)
    for ax, t in zip(axes, ("sWGA", "gDNA")):
        d = res[t]
        lo, hi = d["rest p10"].median(), d["rest p90"].median()
        ax.axhspan(lo, hi, color=BAND, alpha=0.15, lw=0)
        ax.axvspan(0.5, TRIMMED + 0.5, facecolor="none", edgecolor=GRID, hatch="///", lw=0, zorder=0)
        ax.text(TRIMMED / 2 + 0.5, 4.5, "left out of\nthe model fits", ha="center", va="center", fontsize=10.5,
                color=MUTED)
        ax.text(len(cats) + 0.45, hi * 1.08, f"usual range in the rest of the genome\n(10–90% of windows: {lo:.2g}×–{hi:.2g}×)",
                ha="right", va="bottom", fontsize=10.5, color=MUTED)
        data = [np.maximum(d[c].to_numpy(), 1e-4) for c in cats]   # log axis: 0 shown at the floor
        bp = ax.boxplot(data, widths=0.5, patch_artist=True, showfliers=False,
                        medianprops=dict(color=INK, lw=2))
        for p in bp["boxes"]:
            p.set(facecolor=COLOUR[t], alpha=0.55, edgecolor=COLOUR[t])
        for k, v in enumerate(data, 1):
            ax.scatter(np.full(len(v), k) + rng.uniform(-0.18, 0.18, len(v)), v, s=6, color=COLOUR[t],
                       alpha=0.35, lw=0, zorder=3)
            ax.text(k, np.median(v) * (1.25 if k < len(cats) else 0.62), f"{np.median(v):.2g}×", ha="center",
                    fontsize=12, fontweight="bold", color=INK, zorder=5)
        ax.axhline(1, color=MUTED, lw=1, ls="--")
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 8)
        ax.set_yticks([0.001, 0.01, 0.1, 1])
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}×"))
        ax.set_xticks(range(1, len(cats) + 1))
        ax.set_xticklabels(["0–25 kb\nfrom the end", "25–50 kb", "50–100 kb", "100–200 kb", "rest of\nthe genome"])
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.set_title(t, loc="left", fontweight="bold")
    axes[0].set_ylabel("Median depth in the region ÷ sample's median\n(one dot per sample, log scale)")
    fig.suptitle(f"Chromosome ends vs the rest of the genome, {N} random QC-pass Pf8 samples per type "
                 "(both ends of all 14 chromosomes)", x=0.01, ha="left", fontsize=14)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(OUT, f"ends_vs_genome.{ext}"), dpi=200)
    print("wrote", os.path.join(OUT, "ends_vs_genome.{png,pdf,tsv}"))
    plt.close(fig)
    plot_profile(res, profiles)


def plot_profile(res, profiles):
    x = (np.arange(REST_FROM) + 0.5) * 0.5   # kb from the chromosome end (window centres)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.4), sharey=True)
    for ax, t in zip(axes, ("sWGA", "gDNA")):
        d = res[t]
        lo, hi = d["rest p10"].median(), d["rest p90"].median()
        ax.axhspan(lo, hi, color=BAND, alpha=0.15, lw=0, label=f"usual range in the rest of the genome ({lo:.2g}×–{hi:.2g}×)")
        q25, q50, q75 = np.maximum(np.percentile(profiles[t], [25, 50, 75], axis=0), 1e-4)
        ax.fill_between(x, q25, q75, color=COLOUR[t], alpha=0.25, lw=0, label="middle 50% of samples")
        ax.plot(x, q50, color=COLOUR[t], lw=2, label="median sample")
        ax.axhline(1, color=MUTED, lw=1, ls="--")
        marks = [(TRIM_KB, "50 kb: trimmed\nfor fitting")] + ([(END_ZONE_KB, "30 kb: krsim\nend zone")] if t == "sWGA" else [])
        for kb, lab in marks:
            ax.axvline(kb, color=INK, lw=1.2, ls=(0, (4, 3)))
            ax.text(kb + 2, 2e-4 if kb == END_ZONE_KB else 2e-3, lab, fontsize=10.5, color=INK, va="bottom")
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 8)
        ax.set_xlim(0, 200)
        ax.set_yticks([0.001, 0.01, 0.1, 1])
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}×"))
        ax.set_xlabel("kb from the chromosome end")
        ax.grid(color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.set_title(t, loc="left", fontweight="bold")
        ax.legend(frameon=False, loc="lower right", fontsize=10.5)
    axes[0].set_ylabel("Depth ÷ sample's median\n(median over all chromosome ends, log scale)")
    fig.suptitle(f"Read depth by distance from the chromosome end, {len(res['sWGA'])} random QC-pass Pf8 samples "
                 "per type (500 bp windows)", x=0.01, ha="left", fontsize=14)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(OUT, f"ends_profile.{ext}"), dpi=200)
    print("wrote", os.path.join(OUT, "ends_profile.{png,pdf}"))


if __name__ == "__main__":
    main()
