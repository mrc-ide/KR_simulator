"""Mean and SD of read depth in each window, from mosdepth's per-base output.

mosdepth --by gives only the mean per window; the SD needs the depth at every
base, which is in <sample>.per-base.bed.gz (mosdepth run without -n).

Usage:
    python window_stats.py <sample>.per-base.bed.gz <output>.tsv [window_size, default 500]
"""
import sys

import numpy as np
import pandas as pd

per_base_bed, out_tsv = sys.argv[1], sys.argv[2]
window_size = int(sys.argv[3]) if len(sys.argv) > 3 else 500

# each line is a run of bases with the same depth: chromosome, start (0-based), end, depth
runs = pd.read_csv(per_base_bed, sep="\t", header=None, names=["chrom", "start", "end", "depth"])

tables = []
for chrom, r in runs.groupby("chrom", sort=False):
    depth = np.repeat(r["depth"].to_numpy(), (r["end"] - r["start"]).to_numpy())  # depth at every base
    starts = np.arange(0, len(depth), window_size)
    windows = [depth[s:s + window_size] for s in starts]  # the last window may be shorter
    tables.append(pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end": np.minimum(starts + window_size, len(depth)),
        "mean": [w.mean() for w in windows],
        "sd": [w.std() for w in windows],
    }))

stats = pd.concat(tables, ignore_index=True)
stats.round({"mean": 2, "sd": 2}).to_csv(out_tsv, sep="\t", index=False)
print(f"{len(stats)} windows -> {out_tsv}")
