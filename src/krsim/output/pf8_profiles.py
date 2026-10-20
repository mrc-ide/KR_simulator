"""Read-depth profiles from MalariaGEN Pf8 samples, for `krsim vcf-host-depth`.

A profile is the mean and SD of read depth in every 500 bp window of one real
Pf8 sample (<sample>.window_stats.tsv, or .tsv.gz). Profiles are made once and
kept in <pf8_dir>/output/; making one needs mosdepth (`conda install -c bioconda
mosdepth`) and downloads the sample's CRAM (0.5-11 GB), which is deleted after.

Data: MalariaGEN Pf8 (Wellcome Open Research 10:325, 2025). Reads are the Pf8
CRAMs in ENA project PRJEB96783; the reference and sample metadata are on Zenodo.
"""

import os
import shutil
import logging
import subprocess
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

ZENODO = "https://zenodo.org/records/18681980/files"
REF_NAME = "PlasmoDB-54-Pfalciparum3D7-Genome.fasta"
ENA_PROJECT = "PRJEB96783"
WINDOW_SIZE = 500


def _download(url, path):
    """Download to path (via a .part file, so a broken download is never mistaken for a finished one)."""
    with urllib.request.urlopen(url) as response, open(path + ".part", "wb") as f:
        shutil.copyfileobj(response, f, length=1 << 20)
    os.replace(path + ".part", path)


def _find_cram_url(sample):
    """The URL of this sample's Pf8 CRAM in ENA (matched on sample name), or None."""
    query = urllib.parse.urlencode({
        "result": "read_run",
        "query": f'study_accession={ENA_PROJECT} AND sample_alias="{sample}"',
        "fields": "submitted_ftp",
        "format": "tsv",
    })
    with urllib.request.urlopen(f"https://www.ebi.ac.uk/ena/portal/api/search?{query}") as response:
        lines = response.read().decode().strip().split("\n")[1:]
    paths = [p for line in lines for p in line.split("\t")[1].split(";")]
    crams = [p for p in paths if p.endswith(".cram")]
    return f"https://{crams[0]}" if crams else None


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


def _profile_path(outputs, sample):
    """An existing profile for this sample (.tsv or .tsv.gz), or None."""
    for ext in (".window_stats.tsv", ".window_stats.tsv.gz"):
        path = os.path.join(outputs, sample + ext)
        if os.path.exists(path):
            return path
    return None


def _make_profile(sample, inputs, outputs):
    """CRAM from ENA -> mosdepth (depth at every base) -> mean and SD per window; the CRAM is deleted."""
    if shutil.which("mosdepth") is None:
        raise ValueError(f"no depth profile for {sample} yet, and making one needs mosdepth: "
                         "conda install -c bioconda mosdepth")

    ref = os.path.join(inputs, REF_NAME)
    if not os.path.exists(ref):
        log.info("downloading the Pf8 reference genome")
        _download(f"{ZENODO}/{REF_NAME}", ref)

    url = _find_cram_url(sample)
    if url is None:
        raise ValueError(f"{sample} has no Pf8 CRAM in ENA (e.g. Indonesian samples have no reads in the release)")

    cram = os.path.join(inputs, f"{sample}.cram")
    log.info(f"{sample}: downloading {url}")
    _download(url, cram)
    _download(url + ".crai", cram + ".crai")

    log.info(f"{sample}: mosdepth")
    prefix = os.path.join(outputs, sample)
    subprocess.run(["mosdepth", "--fasta", ref, prefix, cram], check=True)

    log.info(f"{sample}: window mean and SD")
    stats = os.path.join(outputs, f"{sample}.window_stats.tsv")
    window_stats(f"{prefix}.per-base.bed.gz", stats)

    os.remove(cram)
    os.remove(cram + ".crai")
    return stats


def get_profiles(samples, pf8_dir):
    """
    Make sure every sample has a depth profile in <pf8_dir>/output/, making any
    that are missing. Returns ({sample: profile path}, the Pf8 sample metadata).
    """
    inputs, outputs = os.path.join(pf8_dir, "input"), os.path.join(pf8_dir, "output")
    os.makedirs(inputs, exist_ok=True)
    os.makedirs(outputs, exist_ok=True)

    meta_path = os.path.join(inputs, "Pf8_samples.txt")
    if not os.path.exists(meta_path):
        log.info("downloading the Pf8 sample metadata")
        _download(f"{ZENODO}/Pf8_samples.txt", meta_path)
    meta = pd.read_csv(meta_path, sep="\t").set_index("Sample")

    unknown = [s for s in samples if s not in meta.index]
    if unknown:
        raise ValueError(f"not Pf8 sample IDs (see the Sample column of {meta_path}): {', '.join(unknown)}")

    profiles = {}
    for sample in samples:
        path = _profile_path(outputs, sample)
        if path is not None:
            log.debug(f"{sample}: profile {path}")
        else:
            if not meta.loc[sample, "QC pass"]:
                log.info(f"{sample}: warning, this sample failed Pf8 QC")
            path = _make_profile(sample, inputs, outputs)
        profiles[sample] = path
    return profiles, meta


def load_profile(path):
    """{Pf8 chromosome name: (mean array, SD array)}, indexed by window number."""
    df = pd.read_csv(path, sep="\t")
    return {chrom: (g["mean"].to_numpy(), g["sd"].to_numpy()) for chrom, g in df.groupby("chrom", sort=False)}


if __name__ == "__main__":
    # make one profile from mosdepth's per-base output (used by real_data/pf8/pf8_mosdepth.sh):
    #   python -m krsim.output.pf8_profiles <sample>.per-base.bed.gz <sample>.window_stats.tsv [window size]
    import sys

    window_stats(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else WINDOW_SIZE)
    print(f"wrote {sys.argv[2]}")
