"""Read-depth profiles for `krsim vcf-host-depth`, and zips of them made by `krsim depth-zip`.

A profile is the mean and SD of one sample's read depth in every 500 bp window of the genome:
<sample>.window_stats.tsv, with columns chrom, start, end, mean and sd. window_stats() makes one
from mosdepth's per-base output (<sample>.per-base.bed.gz).

krsim simulates the 3D7 reference genome used by MalariaGEN Pf8 (PlasmoDB-54), so a profile must
come from reads mapped to that reference: chromosomes Pf3D7_01_v3 ... Pf3D7_14_v3, at their 3D7 lengths.
"""

import os
import zipfile
import logging

import numpy as np
import pandas as pd

from krsim.sim.params import COMP_BREAKPOINTS, GENOME_LENGTH

log = logging.getLogger(__name__)

WINDOW_SIZE = 500
PROFILE_SUFFIXES = (".window_stats.tsv", ".window_stats.tsv.gz")
BED_SUFFIX = ".per-base.bed.gz"
COLUMNS = ["chrom", "start", "end", "mean", "sd"]

# krsim names chromosomes 1-14; profiles use the 3D7 names
CHROMS = [f"Pf3D7_{i:02d}_v3" for i in range(1, 15)]
CHROM_LENGTHS = dict(zip(CHROMS, np.diff(np.r_[0, COMP_BREAKPOINTS, GENOME_LENGTH]).tolist()))


def window_stats(per_base_bed, out_tsv, window_size=WINDOW_SIZE):
    """Mean and SD of depth in each window, from mosdepth's per-base output."""
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
    pd.concat(tables, ignore_index=True).round({"mean": 2, "sd": 2}).to_csv(out_tsv, sep="\t", index=False)


def check_profile(df, name):
    """
    {3D7 chromosome name: (mean array, SD array)}, indexed by window number, after checking that
    the table is a 500 bp window profile of the 3D7 reference. Stops with a plain message if not.
    """
    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is not a depth profile: missing column{'s' if len(missing) > 1 else ''} "
                         f"{', '.join(missing)} (needs {', '.join(COLUMNS)})")
    chroms = {chrom: g.sort_values("start") for chrom, g in df.groupby("chrom", sort=False)}
    absent = [c for c in CHROMS if c not in chroms]
    if absent:
        found = ", ".join(list(chroms)[:3]) + (", ..." if len(chroms) > 3 else "")
        raise ValueError(f"{name} has no {', '.join(absent[:3])}{', ...' if len(absent) > 3 else ''}; krsim needs "
                         f"the 3D7 reference used by Pf8 (chromosomes Pf3D7_01_v3 ... Pf3D7_14_v3). Found: {found}")
    profile = {}
    for chrom in CHROMS:
        g = chroms[chrom]
        n_needed = -(-CHROM_LENGTHS[chrom] // WINDOW_SIZE)
        if not np.array_equal(g["start"].to_numpy()[:2], np.arange(min(2, len(g))) * WINDOW_SIZE):
            raise ValueError(f"{name}: {chrom}'s windows must start at 0, 500, 1000, ... (500 bp windows)")
        if len(g) < n_needed:
            raise ValueError(f"{name}: {chrom} has {len(g):,} windows but the 3D7 reference used by Pf8 needs "
                             f"{n_needed:,} ({CHROM_LENGTHS[chrom]:,} bp); were the reads mapped to another reference?")
        if g[["mean", "sd"]].isna().any().any() or (g[["mean", "sd"]] < 0).any().any():
            raise ValueError(f"{name}: {chrom} has missing or negative mean or SD values")
        profile[chrom] = (g["mean"].to_numpy(dtype=float), g["sd"].to_numpy(dtype=float))
    return profile


def _profile_name(path):
    """The sample name of a profile file (<sample>.window_stats.tsv or .tsv.gz), or None."""
    base = os.path.basename(path)
    for suffix in PROFILE_SUFFIXES:
        if base.endswith(suffix) and len(base) > len(suffix):
            return base[:-len(suffix)]
    return None


def zip_samples(zip_path):
    """{sample: file name inside the zip} for every <sample>.window_stats.tsv in the zip."""
    try:
        with zipfile.ZipFile(zip_path) as z:
            names = z.namelist()
    except zipfile.BadZipFile as exc:
        raise ValueError(f"{zip_path} is not a zip file") from exc
    samples = {}
    for name in names:
        # skip folders, and the hidden extra files macOS adds when it makes a zip
        if name.endswith("/") or name.startswith("__MACOSX/") or os.path.basename(name).startswith("."):
            continue
        sample = _profile_name(name)
        if sample is None:
            continue
        if sample in samples:
            raise ValueError(f"{zip_path} has two depth files for {sample}: {samples[sample]} and {name}")
        samples[sample] = name
    if not samples:
        raise ValueError(f"{zip_path} has no <sample>.window_stats.tsv files (make the zip with `krsim depth-zip`)")
    return samples


def load_zip_profile(zip_path, name):
    """One profile from inside a zip: {3D7 chromosome name: (mean array, SD array)}."""
    with zipfile.ZipFile(zip_path) as z, z.open(name) as f:
        try:
            df = pd.read_csv(f, sep="\t", compression="gzip" if name.endswith(".gz") else None)
        except Exception as exc:
            raise ValueError(f"{name} in {zip_path} is not a readable table ({exc.__class__.__name__})") from exc
    return check_profile(df, f"{name} in {zip_path}")


def make_depth_zip(beds, out_zip, work_dir):
    """
    One profile per mosdepth per-base file (<sample>.per-base.bed.gz), each checked against the 3D7
    reference, all written into out_zip as <sample>.window_stats.tsv. Returns the sample names.
    """
    samples = {}
    for bed in beds:
        base = os.path.basename(bed)
        if not base.endswith(BED_SUFFIX) or len(base) == len(BED_SUFFIX):
            raise ValueError(f"{bed}: expected mosdepth's per-base output, named <sample>{BED_SUFFIX}")
        sample = base[:-len(BED_SUFFIX)]
        if sample in samples:
            raise ValueError(f"two per-base files for {sample}: {samples[sample]} and {bed}")
        if not os.path.exists(bed):
            raise FileNotFoundError(2, "No such file", bed)
        samples[sample] = bed

    os.makedirs(work_dir, exist_ok=True)
    part = out_zip + ".part"  # a zip that stops partway is never mistaken for a finished one
    with zipfile.ZipFile(part, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for i, (sample, bed) in enumerate(samples.items(), 1):
            log.info(f"[{i}/{len(samples)}] {sample}")
            tsv = os.path.join(work_dir, f"{sample}.window_stats.tsv")
            window_stats(bed, tsv)
            check_profile(pd.read_csv(tsv, sep="\t"), bed)
            z.write(tsv, arcname=os.path.basename(tsv))
            os.remove(tsv)
    os.replace(part, out_zip)
    return list(samples)


if __name__ == "__main__":
    # make one profile from mosdepth's per-base output (used by real_data/pf8/extraction/pf8_mosdepth.sh and all_samples/pf8_depth_task.sh):
    #   python -m krsim.output.pf8_profiles <sample>.per-base.bed.gz <sample>.window_stats.tsv [window size]
    import sys

    window_stats(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else WINDOW_SIZE)
    print(f"wrote {sys.argv[2]}")
