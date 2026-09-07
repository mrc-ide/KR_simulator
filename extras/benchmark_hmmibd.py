"""
Benchmark simulated IBD ground truth against hmmIBD's inferred IBD.

Steps
-----
1. Convert output_{seed}_strain.vcf -> hmmIBD hap file + freq file
2. Run hmmIBD
3. Parse hmmIBD fract output -> inferred pairwise frac_IBD
4. Merge with relate_strain_{seed}.csv (simulation ground truth)
5. Save comparison CSV and scatter plot PNG

Installing hmmIBD
-----------------
    git clone https://github.com/glipsnort/hmmIBD.git
    cd hmmIBD
    gcc -o hmmibd hmmIBD.c -lm
    # then either put on PATH or pass --hmmibd /path/to/hmmibd

Usage (from repo root)
----------------------
    python extras/benchmark_hmmibd.py \\
        --strain_vcf  results/output_42_strain.vcf \\
        --relate_strain results/relate_strain_42.csv \\
        --random_seed 42 \\
        --output_dir  results/
"""

import argparse
import itertools
import logging
import os
import subprocess
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark simulated IBD against hmmIBD inferred IBD"
    )
    parser.add_argument("--strain_vcf",     required=True,
                        help="Per-strain VCF (output of run_vcf_strain.py)")
    parser.add_argument("--relate_strain",  required=True,
                        help="Simulated pairwise IBD CSV (output of extras/relate_t_strain.py)")
    parser.add_argument("--random_seed",    type=int, required=True,
                        help="Random seed (used in output filenames)")
    parser.add_argument("--output_dir",     required=True,
                        help="Directory to write all outputs")
    parser.add_argument("--hmmibd",         default="hmmibd",
                        help="Path to hmmIBD binary (default: 'hmmibd' on PATH)")
    parser.add_argument("--chrom",          default="1",
                        help="Chromosome to subset to before running hmmIBD (default: '1')")
    return parser.parse_args()


def _format_elapsed(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}m"
    else:
        return f"{seconds / 3600:.2f}h"


# ---------------------------------------------------------------------------
# Step 1 — build hmmIBD input files from strain VCF
# ---------------------------------------------------------------------------

def _parse_strain_vcf(vcf_path):
    """Return (samples, rows) where rows is a list of dicts with chrom/pos/gts."""
    samples = None
    rows = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = cols[9:]
                continue
            chrom   = cols[0]
            pos     = cols[1]
            gts     = [c.split(":")[0] for c in cols[9:]]
            # convert to 0 / 1 / -1 (missing)
            def to_int(g):
                if g in (".", "./.", ".|."):
                    return -1
                return int(g)
            rows.append({"chrom": chrom, "pos": pos, "gts": [to_int(g) for g in gts]})
    return samples, rows


def build_hap_file(vcf_path, out_path, chrom=None):
    """
    hmmIBD hap file format (tab-separated):
      chrom  pos  sample1  sample2  ...
      1      100  0        1        ...
    Genotypes: 0 = ref, 1 = alt, -1 = missing
    If chrom is given, only rows matching that chromosome are written.
    """
    samples, rows = _parse_strain_vcf(vcf_path)
    if chrom is not None:
        rows = [r for r in rows if r["chrom"] == str(chrom)]
    with open(out_path, "w") as f:
        f.write("\t".join(["chrom", "pos"] + samples) + "\n")
        for r in rows:
            f.write("\t".join([r["chrom"], r["pos"]] + [str(g) for g in r["gts"]]) + "\n")
    logging.info(f"hap file: {len(rows)} sites on chrom {chrom}, {len(samples)} samples -> {out_path}")
    return samples, rows


def build_freq_file(samples, rows, out_path):
    """
    hmmIBD freq file format (tab-separated):
      CHROM  POS  N  FREQ
    N    = total number of haplotypes being analysed
    FREQ = alt allele frequency estimated from the sample
    """
    n = len(samples)
    with open(out_path, "w") as f:
        for r in rows:
            called = [g for g in r["gts"] if g != -1]
            freq   = sum(called) / len(called) if called else 0.5
            f.write(f"{r['chrom']}\t{r['pos']}\t{n}\t{freq:.6f}\n")
    logging.info(f"freq file written -> {out_path}")


# ---------------------------------------------------------------------------
# Step 2 — run hmmIBD
# ---------------------------------------------------------------------------

def run_hmmibd(hmmibd_bin, hap_path, freq_path, prefix):
    cmd = [hmmibd_bin, "-i", hap_path, "-f", freq_path, "-o", prefix]
    logging.info(f"running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        logging.info(result.stdout.strip())
    if result.returncode != 0:
        logging.error(result.stderr)
        raise RuntimeError(
            f"hmmIBD exited with code {result.returncode}.\n"
            f"stderr: {result.stderr}"
        )
    logging.info("hmmIBD completed")


# ---------------------------------------------------------------------------
# Step 3 — parse hmmIBD fract output
# ---------------------------------------------------------------------------

def parse_hmmibd_fract(prefix):
    """
    hmmIBD writes a fract file at {prefix}.hmm_fract.txt.
    Expected columns: samp1, samp2, ..., fract_vit_IBD
    Returns a DataFrame with columns: sample_1, sample_2, inferred_frac_IBD
    """
    fract_path = prefix + ".hmm_fract.txt"
    if not os.path.exists(fract_path):
        raise FileNotFoundError(
            f"hmmIBD fract output not found at {fract_path}. "
            f"Check that hmmIBD ran successfully."
        )

    df = pd.read_csv(fract_path, sep="\t")
    logging.info(f"hmmIBD fract: {len(df)} pairs in {fract_path}")

    # find the IBD fraction column — robust to version differences
    ibd_col = None
    for candidate in ["fract_vit_IBD", "fract_IBD_hmm", "fract_sites_IBD"]:
        if candidate in df.columns:
            ibd_col = candidate
            break
    if ibd_col is None:
        # fall back to any column with 'ibd' in its name (case-insensitive)
        hits = [c for c in df.columns if "ibd" in c.lower() or "fract" in c.lower()]
        if not hits:
            raise ValueError(f"Cannot find IBD fraction column. Columns: {list(df.columns)}")
        ibd_col = hits[-1]

    logging.info(f"using hmmIBD column '{ibd_col}' as inferred IBD fraction")

    # find sample columns
    s1_col = next((c for c in df.columns if c.lower() in ("samp1", "sample1", "sample_1")), df.columns[0])
    s2_col = next((c for c in df.columns if c.lower() in ("samp2", "sample2", "sample_2")), df.columns[1])

    out = df[[s1_col, s2_col, ibd_col]].copy()
    out.columns = ["sample_1", "sample_2", "inferred_frac_IBD"]
    out["sample_1"] = out["sample_1"].astype(str)
    out["sample_2"] = out["sample_2"].astype(str)
    return out


# ---------------------------------------------------------------------------
# Step 4 — merge with ground truth and compute stats
# ---------------------------------------------------------------------------

def build_comparison(truth_csv, hmmibd_df, samples):
    """
    Merge ground truth with hmmIBD inferred values.
    Pairs absent from ground truth -> true IBD = 0.
    Pairs absent from hmmIBD      -> inferred IBD = 0.
    """
    truth = pd.read_csv(truth_csv)
    truth["sample_1"] = truth["sample_1"].astype(str)
    truth["sample_2"] = truth["sample_2"].astype(str)

    # build a canonical key (smaller id first) for both tables
    def canonical(df, c1, c2):
        df = df.copy()
        df["_key"] = df.apply(
            lambda r: (min(r[c1], r[c2]), max(r[c1], r[c2])), axis=1
        )
        return df

    truth_k   = canonical(truth,   "sample_1", "sample_2")
    hmmibd_k  = canonical(hmmibd_df, "sample_1", "sample_2")

    # all possible pairs (from sample list)
    all_pairs = [
        (min(str(a), str(b)), max(str(a), str(b)))
        for a, b in itertools.combinations(samples, 2)
    ]
    all_df = pd.DataFrame(all_pairs, columns=["s1", "s2"])
    all_df["_key"] = list(zip(all_df["s1"], all_df["s2"]))

    truth_map  = dict(zip(truth_k["_key"],  truth_k["frac_IBD"]))
    hmmibd_map = dict(zip(hmmibd_k["_key"], hmmibd_k["inferred_frac_IBD"]))

    all_df["true_frac_IBD"]     = all_df["_key"].map(truth_map).fillna(0.0)
    all_df["inferred_frac_IBD"] = all_df["_key"].map(hmmibd_map).fillna(0.0)

    comp = all_df[["s1", "s2", "true_frac_IBD", "inferred_frac_IBD"]].copy()
    comp.columns = ["sample_1", "sample_2", "true_frac_IBD", "inferred_frac_IBD"]

    r     = np.corrcoef(comp["true_frac_IBD"], comp["inferred_frac_IBD"])[0, 1]
    rmse  = np.sqrt(np.mean((comp["true_frac_IBD"] - comp["inferred_frac_IBD"]) ** 2))
    logging.info(f"benchmark: n_pairs={len(comp)}, Pearson r={r:.3f}, RMSE={rmse:.4f}")

    return comp, r, rmse


# ---------------------------------------------------------------------------
# Step 5 — scatter plot
# ---------------------------------------------------------------------------

def scatter_plot(comp, r, rmse, seed, out_path):
    has_truth = comp["true_frac_IBD"] > 0

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(
        comp.loc[~has_truth, "true_frac_IBD"],
        comp.loc[~has_truth, "inferred_frac_IBD"],
        color="steelblue", alpha=0.6, s=40, label="True IBD = 0",
    )
    ax.scatter(
        comp.loc[has_truth, "true_frac_IBD"],
        comp.loc[has_truth, "inferred_frac_IBD"],
        color="firebrick", alpha=0.8, s=50, label="True IBD > 0",
    )

    lim = max(comp[["true_frac_IBD", "inferred_frac_IBD"]].max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", linewidth=0.8, label="y = x")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)

    ax.set_xlabel("Simulated ground truth frac_IBD")
    ax.set_ylabel("hmmIBD inferred frac_IBD")
    ax.set_title(f"IBD benchmark (seed {seed})\nPearson r = {r:.3f}   RMSE = {rmse:.4f}")
    ax.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    logging.info(f"scatter plot saved to {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    start = time.time()

    os.makedirs(args.output_dir, exist_ok=True)

    log_path = os.path.join(args.output_dir, f"KR_sim_{args.random_seed}.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.FileHandler(log_path, mode="a"), logging.StreamHandler()],
    )

    logging.info("")
    logging.info("========================================")
    logging.info("==== RUNNING extras/benchmark_hmmibd.py ===")
    logging.info("========================================")

    # check hmmIBD binary
    check = subprocess.run(["which", args.hmmibd], capture_output=True)
    if check.returncode != 0:
        logging.error(
            f"hmmIBD binary not found: '{args.hmmibd}'\n"
            "Install with:\n"
            "  git clone https://github.com/glipsnort/hmmIBD.git\n"
            "  cd hmmIBD && gcc -o hmmibd hmmIBD.c -lm\n"
            "Then pass --hmmibd /path/to/hmmibd or add it to your PATH."
        )
        sys.exit(1)

    seed = args.random_seed
    hap_path  = os.path.join(args.output_dir, f"hmmibd_hap_{seed}.txt")
    freq_path = os.path.join(args.output_dir, f"hmmibd_freq_{seed}.txt")
    prefix    = os.path.join(args.output_dir, f"hmmibd_{seed}")

    # step 1 — input files
    samples, rows = build_hap_file(args.strain_vcf, hap_path, chrom=args.chrom)
    build_freq_file(samples, rows, freq_path)

    # step 2 — run hmmIBD
    run_hmmibd(args.hmmibd, hap_path, freq_path, prefix)

    # step 3 — parse output
    hmmibd_df = parse_hmmibd_fract(prefix)

    # step 4 — merge with ground truth
    comp, r, rmse = build_comparison(args.relate_strain, hmmibd_df, samples)

    comp_path = os.path.join(args.output_dir, f"ibd_benchmark_{seed}.csv")
    comp.to_csv(comp_path, index=False)
    logging.info(f"### comparison table - {comp_path}")

    # step 5 — scatter plot
    plot_path = os.path.join(args.output_dir, f"ibd_benchmark_{seed}.png")
    scatter_plot(comp, r, rmse, seed, plot_path)
    logging.info(f"### scatter plot - {plot_path}")

    logging.info(f"extras/benchmark_hmmibd.py COMPLETED in {_format_elapsed(time.time() - start)}")


if __name__ == "__main__":
    main()
