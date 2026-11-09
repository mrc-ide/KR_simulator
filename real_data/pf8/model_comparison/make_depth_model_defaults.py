"""
The defaults for `krsim vcf-host-depth --depth_model sWGA|gDNA`, from every QC-pass Pf8 sample.

Model values: the median over samples of the chosen model's fitted values (shared_params.tsv, made by
fit_chosen_models.py merge): sWGA AR(8) + Student t kicks, gDNA ARMA(2,2) + Student t kicks.

From each sample's window stats (nuclear chromosomes, 500 bp windows), as in the fits:
  typical depth   median window depth                         -> log-normal of hosts' typical depths
  offset          mean of log(depth / median + 0.05), ends trimmed (the fits removed it)
  highest spike   highest window / median                     -> cap: 99th percentile over samples
  end depth       window / median by distance from each chromosome end (5 kb bins, up to 150 kb)
                  -> sWGA end zone: the bins whose median over samples is below END_LOW;
                     factor = median over samples of the zone's windows / median

    python make_depth_model_defaults.py [workers]  # writes depth_model_defaults.json into src/krsim/data/
                                                  # and a per-sample table into $PF8/model_fits/ (~15 min, 8 workers)
    python make_depth_model_defaults.py reuse      # the same from the saved per-sample table (seconds)
"""

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

PF8 = "/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
OUT_JSON = os.path.join(REPO, "src", "krsim", "data", "depth_model_defaults.json")
OUT_TSV = f"{PF8}/model_fits/depth_model_samples.tsv"

WINDOW = 500
PSEUDO = 0.05      # as in the fits (compare_processes.py)
TRIM = 50_000 // WINDOW
END_BIN = 10       # windows per end bin (5 kb)
END_BINS = 30      # bins from each end (150 kb)
END_LOW = 0.1      # an end bin is in the zone if its median depth / sample median is below this (almost no reads)
CAP_PCT = 99
MODELS = {"sWGA": ("AR(8) t", 8, 0), "gDNA": ("ARMA(2,2) t", 2, 2)}


def sample_stats(job):
    sample, sample_type = job
    df = pd.read_csv(f"{PF8}/output/{sample}.window_stats.tsv", sep="\t", usecols=["chrom", "mean"])
    df = df[~df["chrom"].str.contains("API|MIT")]
    m = df["mean"].median()
    if not m > 0:
        return None
    chroms = [g.to_numpy() / m for _, g in df.groupby("chrom", sort=False)["mean"]]
    core = np.concatenate([c[TRIM:len(c) - TRIM] for c in chroms])
    ends = np.zeros(END_BINS)
    for b in range(END_BINS):
        sl = slice(b * END_BIN, (b + 1) * END_BIN)
        ends[b] = np.median(np.concatenate([np.r_[c[sl], c[::-1][sl]] for c in chroms]))
    row = {"sample": sample, "sample type": sample_type, "median depth": m,
           "offset": np.log(core + PSEUDO).mean(), "highest spike": np.concatenate(chroms).max()}
    row.update({f"end bin {b}": v for b, v in enumerate(ends)})
    return row


def main():
    if sys.argv[1:2] == ["reuse"]:
        st = pd.read_csv(OUT_TSV, sep="\t")
        print(f"{len(st):,} samples from {OUT_TSV}")
    else:
        info = pd.read_csv(f"{PF8}/model_fits/samples.tsv", sep="\t")
        info = info[info["Sample type"].isin(MODELS)]
        jobs = list(zip(info["Sample"], info["Sample type"]))
        workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
        print(f"{len(jobs):,} samples, {workers} workers", flush=True)
        with ProcessPoolExecutor(workers) as pool:
            rows = [r for r in pool.map(sample_stats, jobs, chunksize=50) if r is not None]
        st = pd.DataFrame(rows)
        st.round(5).to_csv(OUT_TSV, sep="\t", index=False)

    shared = pd.read_csv(f"{PF8}/model_fits/shared_params.tsv", sep="\t").set_index("sample type")
    defaults = {"source": "MalariaGEN Pf8 QC-pass samples; made by "
                          "real_data/pf8/model_comparison/make_depth_model_defaults.py",
                "window_size": WINDOW, "pseudo": PSEUDO}
    for t, (name, p, q) in MODELS.items():
        s, g = shared.loc[t], st[st["sample type"] == t]
        assert s["model"] == name, (s["model"], name)
        log_m = np.log(g["median depth"])
        # zone: the run of end bins, from the end inwards, whose median depth / sample median is below END_LOW
        bins = np.array([g[f"end bin {b}"].median() for b in range(END_BINS)])
        n_low = int(np.argmax(bins >= END_LOW)) if (bins >= END_LOW).any() else END_BINS
        defaults[t] = {
            "model": name,
            "ar": [float(s[f"ar{i + 1}"]) for i in range(p)],
            "ma": [float(s[f"ma{i + 1}"]) for i in range(q)],
            "mu": float(s["mu"]), "scale": float(s["scale"]), "nu": float(s["nu"]),
            "offset": float(np.median(g["offset"])),
            "depth_mean": float(np.exp(np.median(log_m))),
            "depth_spread": float((log_m.quantile(0.75) - log_m.quantile(0.25)) / 1.349),  # robust SD of log depth
            "cap": float(np.percentile(g["highest spike"], CAP_PCT)),
            "end_width": n_low * END_BIN * WINDOW,
            "end_factor": float(np.median(bins[:n_low])) if n_low else 1.0,
            "n_samples": int(len(g)),
        }
        d = defaults[t]
        print(f"\n== {t} ({len(g):,} samples)")
        print(f"  typical depth: centre {d['depth_mean']:.1f}x, spread (SD of log) {d['depth_spread']:.3f}; "
              f"real 10-90%: {g['median depth'].quantile(0.1):.0f}-{g['median depth'].quantile(0.9):.0f}x")
        print(f"  offset {d['offset']:.3f} | cap {d['cap']:.1f}x the median "
              f"(median sample's highest spike {g['highest spike'].median():.1f}x)")
        print("  end bins (5 kb, from the end):", np.round(bins[:12], 3).tolist())
        print(f"  end zone: {d['end_width']:,} bp at {d['end_factor']:.3f}x")
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(defaults, f, indent=2)
        f.write("\n")
    print(f"\nwrote {OUT_JSON} and {OUT_TSV}")


if __name__ == "__main__":
    main()
