"""
Figure: each host's highest window (/ its median), real Pf8 samples vs krsim's depth model with and without the cap.

N random QC-pass Pf8 samples per type (as in check_depth_cap.py) and N simulated hosts from krsim's depth model,
made without a cap and then capped at the highest real spike (CAP_PCT = 100 in make_depth_model_defaults.py).
Cumulative curves: the share of hosts whose highest window is at or below each value.

Shows that the cap only trims the rare runaway host; typical simulated hosts are unchanged (and their spikes are
still larger than real, the model's known limitation).

    python plot_cap_compare.py [N]     # default 300; writes results/depth_cap/cap_compare.{png,pdf}
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import check_depth_cap as c
from krsim.output import depth_model as dm

FIGURES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "figures", "realism"))

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
COLOUR = {"Real Pf8 samples": "#0f2333", "Model, no cap": "#eb6834", "Model, max cap": "#2a78d6"}


def main():
    os.makedirs(FIGURES, exist_ok=True)
    table = pd.read_csv(f"{c.PF8}/model_fits/depth_model_samples.tsv", sep="\t")
    rng = np.random.default_rng(c.SEED)
    plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, t, xlim in zip(axes, ("sWGA", "gDNA"), ((10, 3e5), (3, 3e3))):
        cap = table.loc[table["sample type"] == t, "highest spike"].max()
        picked = rng.choice(table.loc[table["sample type"] == t, "sample"].to_numpy(), N, replace=False)
        real = c.real_stats(picked)["highest"].to_numpy()
        s = dm.model_settings(t)
        s["cap"] = np.inf
        raw = np.array([max(tr.max() for tr in dm.relative_track(rng, s).values()) for _ in range(N)])
        series = {"Real Pf8 samples": real, "Model, no cap": raw, "Model, max cap": np.minimum(raw, cap)}
        for name, v in series.items():
            x = np.sort(v)
            ax.step(x, np.arange(1, len(x) + 1) / len(x), where="post", color=COLOUR[name],
                    lw=3.5 if name.startswith("Real") else 2.2, zorder=5 if name.endswith("max cap") else 3, label=name)
            print(f"{t} {name:18s} typical {np.median(v):>10,.0f}x  top 10% from {np.percentile(v, 90):>10,.0f}x  "
                  f"most extreme {v.max():>12,.4g}x")
        ax.set_xscale("log")
        ax.set_xlim(*xlim)
        ax.set_ylim(0, 1.02)
        ax.grid(color="#e6e5e0", lw=0.8)
        ax.set_axisbelow(True)
        ax.set_title(f"{t} (max cap {cap:,.0f}×)", loc="left", fontweight="bold")
        ax.set_xlabel("Each host's highest window ÷ its median (log scale)")
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}×"))
    axes[0].set_ylabel("Share of hosts whose highest window\nis at or below this")
    axes[1].legend(frameon=False, loc="lower right", fontsize=12)
    fig.suptitle(f"Highest window per host: real Pf8 samples vs the model with and without the max cap ({N} each)",
                 x=0.01, ha="left", fontsize=14)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(FIGURES, f"cap_compare.{ext}"), dpi=200)
    print("wrote", os.path.join(FIGURES, "cap_compare.{png,pdf}"))


if __name__ == "__main__":
    main()
