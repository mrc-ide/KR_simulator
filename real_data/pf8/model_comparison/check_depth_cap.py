"""
Which cap on the depth model's relative depth makes simulated hosts closest to real Pf8 samples?

Real: N random QC-pass Pf8 samples per type (window stats in $PF8/output/), each window / its sample's median,
nuclear chromosomes. Simulated: N hosts per type from krsim's depth model (krsim.output.depth_model), the track
made without a cap and then capped at each candidate, so every cap is compared on the same tracks.

Per sample / host: the highest window, the 99.9th and 99th percentile window, and the share of windows above
10x and 100x the median. Each is compared between real and simulated by the median and 90th percentile over
samples, and by the KS distance (0 = the same distribution, 1 = no overlap).

    python check_depth_cap.py [N]        # default 300; writes results/depth_cap/
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

from krsim.output import depth_model as dm

PF8 = "/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "depth_cap")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
SEED = 2026


def stats(rel):
    """Summary of one genome's relative depth (window / median)."""
    return {"highest": rel.max(), "p99.9": np.percentile(rel, 99.9), "p99": np.percentile(rel, 99),
            "share >10x": np.mean(rel > 10), "share >100x": np.mean(rel > 100)}


def real_stats(samples):
    rows = []
    for s in samples:
        df = pd.read_csv(f"{PF8}/output/{s}.window_stats.tsv", sep="\t", usecols=["chrom", "mean"])
        df = df[~df["chrom"].str.contains("API|MIT")]
        m = df["mean"].median()
        if m > 0:
            rows.append(stats(df["mean"].to_numpy() / m))
    return pd.DataFrame(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    table = pd.read_csv(f"{PF8}/model_fits/depth_model_samples.tsv", sep="\t")
    rng = np.random.default_rng(SEED)
    lines, rows = [], []
    for t in ("sWGA", "gDNA"):
        spikes = table.loc[table["sample type"] == t, "highest spike"]
        caps = {"none": np.inf, "highest (now)": spikes.max(), "99th pct (before)": spikes.quantile(0.99),
                "90th pct": spikes.quantile(0.90), "median sample": spikes.median()}

        picked = rng.choice(table.loc[table["sample type"] == t, "sample"].to_numpy(), N, replace=False)
        real = real_stats(picked)

        s = dm.model_settings(t)
        s["cap"] = np.inf
        tracks = [np.concatenate(list(dm.relative_track(rng, s).values())) for _ in range(N)]

        lines.append(f"\n=== {t}: {len(real)} real samples, {N} simulated hosts ===")
        lines.append("caps: " + ", ".join(f"{k} {v:,.0f}x" for k, v in caps.items()))
        for name, cap in [("REAL", None)] + list(caps.items()):
            d = real if cap is None else pd.DataFrame([stats(np.minimum(tr, cap)) for tr in tracks])
            row = {"type": t, "cap": name, "cap value": None if cap is None else cap}
            for col in d.columns:
                row[f"{col} median"] = d[col].median()
                row[f"{col} 90th"] = d[col].quantile(0.9)
                if cap is not None:
                    row[f"{col} KS"] = ks_2samp(d[col], real[col]).statistic
            rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "cap_check.tsv"), sep="\t", index=False)

    pd.set_option("display.width", 250)
    for t in ("sWGA", "gDNA"):
        g = res[res["type"] == t].set_index("cap")
        print(f"\n=== {t} ===")
        show = g[[c for c in g.columns if c.endswith(("median", "90th"))]]
        print(show.round(4).T.to_string())
        ks = g[[c for c in g.columns if c.endswith("KS")]].drop(index="REAL")
        ks["mean KS"] = ks.mean(axis=1)
        print("\nKS distance to real (lower = closer):")
        print(ks.round(3).to_string())


if __name__ == "__main__":
    main()
