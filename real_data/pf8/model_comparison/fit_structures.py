"""
Which correlation structure is best, on every QC-pass Pf8 gDNA and sWGA sample.

The structures of best_simple_model.py plus two bigger ones (AR(15), ARMA(3,3)), each with Gaussian and t kicks
(the t fit starts from the Gaussian one). AR(8) and ARMA(2,2) are already fitted by fit_chosen_models.py and are
taken from its fits.tsv; everything is scored on the same windows, so held-out likelihood, BIC and AIC compare
directly. Structures are ranked separately for Gaussian and for t kicks, so we can see whether the choice depends
on the kicks.

    python fit_structures.py fit <task> <n_tasks>   # one SLURM array task (fit_structures.sbatch)
    python fit_structures.py merge                  # -> OUT/structures/structures.tsv and structure_summary.tsv
"""

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from best_simple_model import CLOSE, label
from compare_arma_kicks import score, unpack
from fit_chosen_models import OUT, depth_path, load, samples

STRUCT_OUT = f"{OUT}/structures"
NEW = [(1, 0), (2, 0), (4, 0), (15, 0), (1, 1), (2, 1), (3, 3)]
DONE = [(8, 0), (2, 2)]                                  # from fit_chosen_models.py
STRUCTURES = [(1, 0), (2, 0), (4, 0), (8, 0), (15, 0), (1, 1), (2, 1), (2, 2), (3, 3)]
KICKS = ["Gaussian", "t"]


def fit_sample(job):
    sample, sample_type = job
    t0 = time.time()
    ys = load(depth_path(sample))
    rows = []
    for p, q in NEW:
        start = None
        for kick in KICKS:
            x, r = score(ys, p, q, kick, start)
            start = x[:p + q]                    # the t fit starts from the Gaussian fit's correlation part
            rows.append({"sample": sample, "sample type": sample_type, "model": f"{label(p, q)} {kick}",
                         "correlation": label(p, q), "kicks": kick, **r,
                         "nu": unpack(kick, x[p + q:]).get("nu", np.nan)})
    return rows, time.time() - t0


def run_task(task, n_tasks):
    os.makedirs(f"{STRUCT_OUT}/tasks", exist_ok=True)
    out = f"{STRUCT_OUT}/tasks/structures_{task:04d}.tsv"
    info = samples().iloc[task::n_tasks]
    done = set()
    if os.path.exists(out):
        counts = pd.read_csv(out, sep="\t", usecols=["sample"])["sample"].value_counts()
        done = set(counts.index[counts == len(NEW) * len(KICKS)])
    jobs = [(s, t) for s, t in zip(info["Sample"], info["Sample type"]) if s not in done]
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", max(1, (os.cpu_count() or 2) - 1)))
    print(f"task {task} of {n_tasks}: {len(info)} samples, {len(done)} already done, {workers} workers", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as pool:
        futures = {pool.submit(fit_sample, j): j for j in jobs}
        for f in as_completed(futures):
            sample, sample_type = futures[f]
            try:
                rows, secs = f.result()
            except Exception as e:                # one bad sample should not stop the task
                print(f"  {sample} ({sample_type}): FAILED {type(e).__name__}: {e}", flush=True)
                with open(f"{STRUCT_OUT}/tasks/failed_{task:04d}.tsv", "a") as fh:
                    fh.write(f"{sample}\t{sample_type}\t{type(e).__name__}: {e}\n")
                continue
            pd.DataFrame(rows).to_csv(out, sep="\t", index=False, mode="a", header=not os.path.exists(out))
            print(f"  {sample} ({sample_type}): {secs:.0f} s", flush=True)
    print(f"task {task} done in {time.time() - t0:.0f} s", flush=True)


def merge():
    files = sorted(f for f in os.listdir(f"{STRUCT_OUT}/tasks") if f.startswith("structures_"))
    new = pd.concat([pd.read_csv(f"{STRUCT_OUT}/tasks/{f}", sep="\t") for f in files], ignore_index=True)
    old = pd.read_csv(f"{OUT}/fits.tsv", sep="\t")
    old = old[old["correlation"].isin([label(p, q) for p, q in DONE]) & old["sample"].isin(new["sample"])]
    cols = ["sample", "sample type", "model", "correlation", "kicks", "LL", "k", "n", "held-out LL", "nu"]
    res = pd.concat([new[cols], old[cols]], ignore_index=True).drop_duplicates(["sample", "model"])
    res["AIC"] = 2 * res["k"] - 2 * res["LL"]
    res["BIC"] = res["k"] * np.log(res["n"]) - 2 * res["LL"]
    res["Heldout"] = -2 * res["held-out LL"]
    for s in ("AIC", "BIC", "Heldout"):          # behind the best structure with the same kicks
        res[f"d{s}"] = res[s] - res.groupby(["sample", "kicks"])[s].transform("min")
    res["close to best"] = res["dHeldout"] <= CLOSE
    res.round(4).to_csv(f"{STRUCT_OUT}/structures.tsv", sep="\t", index=False)

    wins = lambda v: int((v == 0).sum())
    summ = res.groupby(["sample type", "kicks", "correlation"]).agg(
        samples=("sample", "nunique"), k=("k", "first"), best_heldout=("dHeldout", wins),
        close_to_best=("close to best", "sum"), median_dHeldout=("dHeldout", "median"),
        best_BIC=("dBIC", wins), median_dBIC=("dBIC", "median"), best_AIC=("dAIC", wins),
    ).round(1).reset_index()
    summ.to_csv(f"{STRUCT_OUT}/structure_summary.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        for (t, k), g in summ.groupby(["sample type", "kicks"]):
            print(f"\n== {t}, {k} kicks")
            print(g.drop(columns=["sample type", "kicks"]).sort_values("median_dHeldout").to_string(index=False))
    print(f"\nwrote structures.tsv and structure_summary.tsv in {STRUCT_OUT}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["fit"] and len(sys.argv) == 4:
        run_task(int(sys.argv[2]), int(sys.argv[3]))
    elif sys.argv[1:] == ["merge"]:
        merge()
    else:
        sys.exit(__doc__)
