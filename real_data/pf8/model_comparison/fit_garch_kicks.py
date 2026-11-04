"""
GARCH t kicks on the chosen structure, for every sample fitted by fit_chosen_models.py.

Gaussian and t kicks have the same spread everywhere along the genome. GARCH t kicks are t kicks whose
scale follows the recent kicks: scale_t^2 = w + a * kick_{t-1}^2 + b * scale_{t-1}^2, so a big kick makes
the next ones bigger, then it settles back (calm and jumpy stretches). The fit starts from the sample's
t fit and is scored on the same windows, so its held-out likelihood and BIC compare directly with the
Gaussian and t rows in fits.tsv.

Clustering of big kicks: the lag-1 correlation of |kick| (kicks divided by the model's scale at that window).
0 = big kicks are scattered at random, as Gaussian and t kicks assume.

    python fit_garch_kicks.py fit <task> <n_tasks>   # one SLURM array task (fit_garch_kicks.sbatch)
    python fit_garch_kicks.py merge                  # -> OUT/garch/kick_compare.tsv (Gaussian, t, GARCH t per sample)
"""

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from best_simple_model import CLOSE, label
from compare_arma_kicks import arma_kicks, garch_s2, score, unpack
from fit_chosen_models import CHOSEN, OUT, depth_path, load, samples, tails

GARCH_OUT = f"{OUT}/garch"
KICKS = ["Gaussian", "t", "GARCH t"]


def standardised(ys, x, p, q, kick):
    """Kicks per chromosome, divided by the model's scale at each window."""
    d = unpack(kick, x[p + q:])
    es = [e - d["mu"] for e in arma_kicks(ys, x, p, q)]
    if kick == "GARCH t":
        return [e / np.sqrt(garch_s2(e, d["w"], d["a"], d["b"])) for e in es]
    return [e / d.get("sd", d.get("scale")) for e in es]


def clustering(zs, lag=1):
    """Correlation of |kick| with |kick| lag windows later, within chromosomes."""
    a = np.concatenate([np.abs(z[:-lag]) for z in zs])
    b = np.concatenate([np.abs(z[lag:]) for z in zs])
    return np.corrcoef(a, b)[0, 1]


def chosen_x():
    """x of the Gaussian and t fits on the chosen structure, per sample (from fit_chosen_models.py's task files)."""
    files = sorted(f for f in os.listdir(f"{OUT}/tasks") if f.startswith("fits_"))
    xs = {}
    for f in files:
        df = pd.read_csv(f"{OUT}/tasks/{f}", sep="\t", usecols=["sample", "sample type", "correlation", "kicks", "x"])
        for r in df.itertuples(index=False):
            (p, q), _ = CHOSEN[r[1]]
            if r.correlation == label(p, q):
                xs.setdefault(r.sample, {})[r.kicks] = np.array(r.x.split(","), dtype=float)
    return xs


def fit_sample(job):
    sample, sample_type, xs = job
    t0 = time.time()
    (p, q), _ = CHOSEN[sample_type]
    ys = load(depth_path(sample))
    x, r = score(ys, p, q, "GARCH t", xs["t"][:p + q])
    d = unpack("GARCH t", x[p + q:])
    row = {"sample": sample, "sample type": sample_type, "model": f"{label(p, q)} GARCH t", **r, **d,
           "persistence": d["a"] + d["b"], **tails(ys, x, p, q, "GARCH t"), "x": ",".join(f"{v:.8g}" for v in x)}
    for kick, xk in (("Gaussian", xs["Gaussian"]), ("t", xs["t"]), ("GARCH t", x)):
        zs = standardised(ys, xk, p, q, kick)
        row[f"clustering {kick}"] = clustering(zs)
        row[f"clustering lag 10 {kick}"] = clustering(zs, 10)
    return row, time.time() - t0


def run_task(task, n_tasks):
    os.makedirs(f"{GARCH_OUT}/tasks", exist_ok=True)
    out = f"{GARCH_OUT}/tasks/garch_{task:04d}.tsv"
    info = samples().iloc[task::n_tasks]
    done = set(pd.read_csv(out, sep="\t", usecols=["sample"])["sample"]) if os.path.exists(out) else set()
    xs = chosen_x()
    jobs = [(s, t, xs[s]) for s, t in zip(info["Sample"], info["Sample type"]) if s not in done and s in xs]
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", max(1, (os.cpu_count() or 2) - 1)))
    print(f"task {task} of {n_tasks}: {len(info)} samples, {len(done)} already done, {workers} workers", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as pool:
        futures = {pool.submit(fit_sample, j): j[:2] for j in jobs}
        for f in as_completed(futures):
            sample, sample_type = futures[f]
            try:
                row, secs = f.result()
            except Exception as e:                # one bad sample should not stop the task
                print(f"  {sample} ({sample_type}): FAILED {type(e).__name__}: {e}", flush=True)
                with open(f"{GARCH_OUT}/tasks/failed_{task:04d}.tsv", "a") as fh:
                    fh.write(f"{sample}\t{sample_type}\t{type(e).__name__}: {e}\n")
                continue
            pd.DataFrame([row]).to_csv(out, sep="\t", index=False, mode="a", header=not os.path.exists(out))
            print(f"  {sample} ({sample_type}): {secs:.0f} s", flush=True)
    print(f"task {task} done in {time.time() - t0:.0f} s", flush=True)


def merge():
    files = sorted(f for f in os.listdir(f"{GARCH_OUT}/tasks") if f.startswith("garch_"))
    g = pd.concat([pd.read_csv(f"{GARCH_OUT}/tasks/{f}", sep="\t") for f in files], ignore_index=True)
    g = g.drop_duplicates("sample").drop(columns="x")
    g["kicks"] = "GARCH t"
    res = pd.read_csv(f"{OUT}/fits.tsv", sep="\t")
    chosen_corr = res["sample type"].map({t: label(p, q) for t, ((p, q), _) in CHOSEN.items()})
    res = res[(res["correlation"] == chosen_corr) & res["sample"].isin(g["sample"])]
    clus = g.set_index("sample")
    for kick in ("Gaussian", "t"):
        m = res["kicks"] == kick
        res.loc[m, "clustering"] = res.loc[m, "sample"].map(clus[f"clustering {kick}"]).to_numpy()
        res.loc[m, "clustering lag 10"] = res.loc[m, "sample"].map(clus[f"clustering lag 10 {kick}"]).to_numpy()
    g["clustering"] = g["clustering GARCH t"]
    g["clustering lag 10"] = g["clustering lag 10 GARCH t"]
    keep = ["sample", "sample type", "model", "kicks", "LL", "k", "n", "held-out LL", "nu", "w", "a", "b", "persistence",
            *[c for c in res if "expect" in c or c == "KS distance"], "clustering", "clustering lag 10"]
    kc = pd.concat([res, g], ignore_index=True).reindex(columns=keep)
    kc["BIC"] = kc["k"] * np.log(kc["n"]) - 2 * kc["LL"]
    kc["Heldout"] = -2 * kc["held-out LL"]
    for s in ("BIC", "Heldout"):
        kc[f"d{s}"] = kc[s] - kc.groupby("sample")[s].transform("min")
    kc["close to best"] = kc["dHeldout"] <= CLOSE
    kc.round(5).to_csv(f"{GARCH_OUT}/kick_compare.tsv", sep="\t", index=False)

    summ = kc.groupby(["sample type", "kicks"]).agg(
        samples=("sample", "nunique"), best_heldout=("dHeldout", lambda v: int((v == 0).sum())),
        best_BIC=("dBIC", lambda v: int((v == 0).sum())), median_dHeldout=("dHeldout", "median"),
        median_nu=("nu", "median"), median_clustering=("clustering", "median"),
        median_drops=("lowest 0.1% (expect 0.1)", "median"), median_spikes=("highest 0.1% (expect 0.1)", "median"),
    ).round(3)
    summ.to_csv(f"{GARCH_OUT}/kick_summary.tsv", sep="\t")
    wide = kc.pivot(index="sample", columns="kicks", values="Heldout")
    gain = pd.DataFrame({"sample type": kc.drop_duplicates("sample").set_index("sample")["sample type"],
                         "t over Gaussian": wide["Gaussian"] - wide["t"],
                         "GARCH t over t": wide["t"] - wide["GARCH t"]})
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(summ.to_string())
        print("\n== held-out gain per sample (-2 LL scale; > 0 = the second model is better)")
        print(gain.groupby("sample type").describe(percentiles=[0.1, 0.5, 0.9]).round(1).T.to_string())
    print(f"\nwrote kick_compare.tsv and kick_summary.tsv in {GARCH_OUT}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["fit"] and len(sys.argv) == 4:
        run_task(int(sys.argv[2]), int(sys.argv[3]))
    elif sys.argv[1:] == ["merge"]:
        merge()
    else:
        sys.exit(__doc__)
