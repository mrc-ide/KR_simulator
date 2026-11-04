"""
The chosen depth models fitted to every QC-pass Pf8 gDNA and sWGA sample.

Chosen by best_simple_model.py on 41 samples: gDNA ARMA(2,2) + t kicks, sWGA AR(8) + t kicks.
Every sample gets both structures, each with Gaussian and t kicks (the t fit starts from the Gaussian one,
as in best_simple_model.py), so we see whether the choice holds, what the t kicks gain, and how well the
chosen model is calibrated. Model choice: held-out likelihood (fit on odd chromosomes, score even, then swap),
checked by BIC; AIC is written but not used.

Then the test for what krsim will ship: one parameter set per sample type (the median over samples of the
chosen model's parameters) scored on every sample without refitting, next to the sample's own fit.

    python fit_chosen_models.py fit <task> <n_tasks>   # samples task, task + n_tasks, ... (one SLURM array task)
    python fit_chosen_models.py merge                  # every task's fits -> tables, shared-parameter scores

Paths are set below (PF8, OUT). A task appends each sample as it finishes and skips samples already done,
so a task that ran out of time can be run again.
"""

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy import stats

from best_simple_model import CLOSE, label
from compare_arma_kicks import arma_kicks, kick_cdf, kick_logpdf, load, score, unpack
from compare_processes import pacf_to_coefs

PF8 = "/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth"
OUT = f"{PF8}/model_fits"

STRUCTURES = [(2, 2), (8, 0)]
KICKS = ["Gaussian", "t"]
CHOSEN = {"gDNA": ((2, 2), "t"), "sWGA": ((8, 0), "t")}
META = ["Country", "Population", "Year", "Study", "% callable", "cram_bytes"]


def depth_path(sample):
    return f"{PF8}/output/{sample}.window_stats.tsv"


def samples():
    """QC-pass gDNA and sWGA samples that have window stats."""
    info = pd.read_csv(f"{PF8}/output/Pf8_sample_info.tsv", sep="\t")
    info = info[info["QC pass"] & info["Sample type"].isin(CHOSEN)]
    info = info[[os.path.exists(depth_path(s)) for s in info["Sample"]]]
    return info.sort_values("Sample").reset_index(drop=True)


def coefs(x, p, q):
    """AR coefficients (y_t = sum ar_i y_t-i + ...) and MA coefficients (... + e_t + sum ma_i e_t-i)."""
    return pacf_to_coefs(np.tanh(x[:p])), -pacf_to_coefs(np.tanh(x[p:p + q]))


def tails(ys, x, p, q, kick):
    """Where the kicks land in the model's own distribution: % in each tail (expect 0.1 and 1) and KS distance."""
    u = np.concatenate([kick_cdf(kick, e, x[p + q:]) for e in arma_kicks(ys, x, p, q)])
    return {"lowest 0.1% (expect 0.1)": 100 * np.mean(u < 0.001), "highest 0.1% (expect 0.1)": 100 * np.mean(u > 0.999),
            "lowest 1% (expect 1)": 100 * np.mean(u < 0.01), "highest 1% (expect 1)": 100 * np.mean(u > 0.99),
            "KS distance": stats.kstest(u, "uniform").statistic}


def fit_sample(job):
    sample, sample_type = job
    t0 = time.time()
    ys = load(depth_path(sample))
    rows = []
    for p, q in STRUCTURES:
        start = None
        for kick in KICKS:
            x, r = score(ys, p, q, kick, start)
            start = x[:p + q]                    # the t fit starts from the Gaussian fit's correlation part
            ar, ma = coefs(x, p, q)
            rows.append({"sample": sample, "sample type": sample_type, "model": f"{label(p, q)} {kick}",
                         "correlation": label(p, q), "kicks": kick, "chosen": ((p, q), kick) == CHOSEN[sample_type],
                         **r, **unpack(kick, x[p + q:]),
                         **{f"ar{i + 1}": c for i, c in enumerate(ar)}, **{f"ma{i + 1}": c for i, c in enumerate(ma)},
                         **tails(ys, x, p, q, kick), "x": ",".join(f"{v:.8g}" for v in x)})
    return rows, time.time() - t0


def run_task(task, n_tasks):
    os.makedirs(f"{OUT}/tasks", exist_ok=True)
    out = f"{OUT}/tasks/fits_{task:04d}.tsv"
    info = samples().iloc[task::n_tasks]
    done = set()
    if os.path.exists(out):
        n_models = len(STRUCTURES) * len(KICKS)
        counts = pd.read_csv(out, sep="\t", usecols=["sample"])["sample"].value_counts()
        done = set(counts.index[counts == n_models])
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
                with open(f"{OUT}/tasks/failed_{task:04d}.tsv", "a") as fh:
                    fh.write(f"{sample}\t{sample_type}\t{type(e).__name__}: {e}\n")
                continue
            df = pd.DataFrame(rows)
            df.to_csv(out, sep="\t", index=False, mode="a", header=not os.path.exists(out))
            print(f"  {sample} ({sample_type}): {secs:.0f} s", flush=True)
    print(f"task {task} done in {time.time() - t0:.0f} s", flush=True)


def shared_score(job):
    """The sample's LL and tails under fixed (shared) parameters, no refitting."""
    sample, x, p, q, kick = job
    ys = load(depth_path(sample))
    ll = sum(kick_logpdf(kick, e, x[p + q:]).sum() for e in arma_kicks(ys, x, p, q))
    return {"sample": sample, "shared LL": ll, **tails(ys, x, p, q, kick)}


def merge():
    files = sorted(f for f in os.listdir(f"{OUT}/tasks") if f.startswith("fits_"))
    res = pd.concat([pd.read_csv(f"{OUT}/tasks/{f}", sep="\t") for f in files], ignore_index=True)
    res = res.drop_duplicates(["sample", "model"])
    res["AIC"] = 2 * res["k"] - 2 * res["LL"]
    res["BIC"] = res["k"] * np.log(res["n"]) - 2 * res["LL"]
    res["Heldout"] = -2 * res["held-out LL"]
    for s in ("AIC", "BIC", "Heldout"):
        res[f"d{s}"] = res[s] - res.groupby("sample")[s].transform("min")
    res["close to best"] = res["dHeldout"] <= CLOSE
    res.drop(columns="x").round(5).to_csv(f"{OUT}/fits.tsv", sep="\t", index=False)
    info = samples().set_index("Sample")
    info.loc[info.index.isin(res["sample"]), ["Sample type", *META]].to_csv(f"{OUT}/samples.tsv", sep="\t")

    q90 = lambda v: v.quantile(0.9)
    summ = res.groupby(["sample type", "model"]).agg(
        samples=("sample", "nunique"), chosen=("chosen", "first"), k=("k", "first"),
        best_heldout=("dHeldout", lambda v: int((v == 0).sum())), close_to_best=("close to best", "sum"),
        median_dHeldout=("dHeldout", "median"), q90_dHeldout=("dHeldout", q90),
        best_BIC=("dBIC", lambda v: int((v == 0).sum())), median_dBIC=("dBIC", "median"),
        median_nu=("nu", "median"),
    ).reset_index()
    cal = res.groupby(["sample type", "model"])[[c for c in res if "expect" in c or c == "KS distance"]].median()
    summ = summ.merge(cal.reset_index(), on=["sample type", "model"]).round(3)
    summ.to_csv(f"{OUT}/summary.tsv", sep="\t", index=False)

    # what t kicks gain over Gaussian on the chosen structure (held-out, per sample)
    gain = []
    for t, ((p, q), _) in CHOSEN.items():
        g = res[(res["sample type"] == t) & (res["correlation"] == label(p, q))].pivot(index="sample", columns="kicks",
                                                                                      values="Heldout")
        gain.append(pd.DataFrame({"sample type": t, "t gain": g["Gaussian"] - g["t"]}))
    gain = pd.concat(gain)

    # one parameter set per type: median over samples of the chosen model's parameters (on the fitting scale,
    # so the median is still a stationary, invertible model), scored on every sample without refitting
    params, shared = [], []
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", max(1, (os.cpu_count() or 2) - 1)))
    for t, ((p, q), kick) in CHOSEN.items():
        own = res[(res["sample type"] == t) & res["chosen"]].set_index("sample")
        x = np.median(np.vstack([np.array(s.split(","), dtype=float) for s in own["x"]]), axis=0)
        ar, ma = coefs(x, p, q)
        params.append({"sample type": t, "model": f"{label(p, q)} {kick}", **unpack(kick, x[p + q:]),
                       **{f"ar{i + 1}": c for i, c in enumerate(ar)}, **{f"ma{i + 1}": c for i, c in enumerate(ma)},
                       "x": ",".join(f"{v:.8g}" for v in x)})
        with ProcessPoolExecutor(workers) as pool:
            sc = pd.DataFrame(pool.map(shared_score, [(s, x, p, q, kick) for s in own.index], chunksize=20))
        sc = sc.set_index("sample")
        sc.insert(0, "sample type", t)
        sc.insert(1, "own LL", own["LL"])
        sc.insert(3, "n", own["n"])
        sc["loss per window"] = (sc["own LL"] - sc["shared LL"]) / sc["n"]
        shared.append(sc)
    pd.DataFrame(params).to_csv(f"{OUT}/shared_params.tsv", sep="\t", index=False)
    shared = pd.concat(shared)
    shared.round(5).to_csv(f"{OUT}/shared_scores.tsv", sep="\t")

    with pd.option_context("display.width", 260, "display.max_columns", 30):
        for t, g in summ.groupby("sample type"):
            print(f"\n== {t}")
            print(g.drop(columns="sample type").sort_values("median_dHeldout").to_string(index=False))
        print("\n== t kicks vs Gaussian on the chosen structure (held-out gain; > 0 = t better)")
        print(gain.groupby("sample type")["t gain"].describe(percentiles=[0.1, 0.5, 0.9]).round(1).to_string())
        print("\n== shared parameters")
        print(pd.DataFrame(params).drop(columns="x").round(4).to_string(index=False))
        print("\n== shared parameters vs the sample's own fit (LL lost per window; tails under shared parameters)")
        print(shared.groupby("sample type")[["loss per window", "lowest 0.1% (expect 0.1)", "highest 0.1% (expect 0.1)",
                                             "KS distance"]].describe(percentiles=[0.1, 0.5, 0.9]).round(3).T.to_string())
    print(f"\nwrote fits.tsv, samples.tsv, summary.tsv, shared_params.tsv and shared_scores.tsv in {OUT}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["fit"] and len(sys.argv) == 4:
        run_task(int(sys.argv[2]), int(sys.argv[3]))
    elif sys.argv[1:] == ["merge"]:
        merge()
    else:
        sys.exit(__doc__)
