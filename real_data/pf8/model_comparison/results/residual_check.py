"""Does the best model actually fit? Check the kicks it leaves behind (should be Gaussian and independent)."""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import lfilter
sys.path.insert(0, "/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/TSKIT_MALARIA_SIM/git/KR_simulator/real_data/pf8/model_comparison")
from compare_processes import K, fit, lag_matrix, load, models

PF8 = "/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/pf8_samples"
info = pd.read_csv(f"{PF8}/samples.tsv", sep="\t")
rows = []
for s, t in zip(info["Sample"], info["Sample type"]):
    name = "AR(8)" if t == "sWGA" else "ARMA(2,2)"
    build, start, _ = models()[name]
    ys = load(f"{PF8}/{s}/depth_plots/{s}.window_stats.tsv")
    Ls = [lag_matrix(y) for y in ys]
    x, _ = fit(ys, Ls, build, start(ys))
    phi, v = build(x)
    z = []
    for y, L in zip(ys, Ls):                       # standardised kicks: (real - prediction) / expected spread
        e = lfilter(np.r_[1, -phi[K, 1:]], [1], y)
        e[:K] = y[:K] - (phi[:K, 1:] * L).sum(1)
        z.append(e / np.sqrt(np.r_[v[:K], np.full(len(y) - K, v[K])]))
    zc = np.concatenate(z)
    lag1 = np.mean([np.corrcoef(a[:-1], a[1:])[0, 1] for a in z])
    lag5 = np.mean([np.corrcoef(a[:-5], a[5:])[0, 1] for a in z])
    # does the size of the kicks cluster? (correlation of |kick| between neighbours; 0 if the spread is constant)
    vol = np.mean([np.corrcoef(np.abs(a[:-1]), np.abs(a[1:]))[0, 1] for a in z])
    rows.append({"sample": s, "type": t, "model": name,
                 "skew": stats.skew(zc), "excess kurtosis": stats.kurtosis(zc),
                 "beyond 3 SD (%)": 100 * np.mean(np.abs(zc) > 3), "beyond 3 SD if Gaussian (%)": 0.27,
                 "below -3 SD (%)": 100 * np.mean(zc < -3), "above +3 SD (%)": 100 * np.mean(zc > 3),
                 "kick corr lag 1": lag1, "kick corr lag 5": lag5, "|kick| corr lag 1": vol})
    print(s, flush=True)
r = pd.DataFrame(rows)
r.round(3).to_csv("residual_check.tsv", sep="\t", index=False)
with pd.option_context("display.width", 250, "display.max_columns", 20):
    print(r.drop(columns=["beyond 3 SD if Gaussian (%)"]).round(2).to_string(index=False))
    print(r.groupby("type").median(numeric_only=True).round(2).to_string())
