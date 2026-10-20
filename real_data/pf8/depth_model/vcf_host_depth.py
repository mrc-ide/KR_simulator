"""Host VCF with realistic read depth from real MalariaGEN Pf8 samples.

A more realistic alternative to krsim's host VCF (src/krsim/output/vcf_host_single.py,
which stays as it is): instead of one constant depth, each host's depth at each SNP
follows a real Pf8 sample's sequencing depth.

Part 1 - depth profiles (skipped for samples already done):
  for each Pf8 sample ID given: find its CRAM in ENA, download it, run mosdepth
  (depth at every base), run window_stats.py (mean and SD per 500 bp window),
  then delete the CRAM.

Part 2 - host VCF:
  once per host:  give it one sample's profile (in order: 1st sample -> 1st host)
                  and draw its strain mix from Dirichlet(alpha = 1)
  for every SNP in every host:
    1. find the 500 bp window the SNP is in -> mean m, SD s
    2. draw the host's total depth: negative binomial with that mean and SD
       (r = m^2 / (s^2 - m), p = r / (r + m)); Poisson(m) if s^2 <= m; 0 if m = 0
    3. split the depth between the host's strains: Multinomial(depth, mix)
    4. add up reads by allele -> GT:AD:DP (GT = the two alleles with the most reads)

Outputs (in --output_dir):
  output_<seed>_host_depth.vcf   host VCF
  host_profiles_<seed>.csv       which Pf8 sample each host's depth came from
  host_proportions_<seed>.csv    each host's strain mix

Usage (with `conda activate mosdepth`; needs numpy and pandas):
  python vcf_host_depth.py --strain_vcf output_123_strain.vcf \\
      --ind_node_map_csv sim_ind_node_map_123.csv \\
      --samples SPT52544 PF1066-Cx RCN08168 \\
      --pf8_dir /Users/jbansal/Desktop/IMPERIAL/pf8 \\
      --random_seed 123 --output_dir results/
"""
import argparse
import csv
import os
import subprocess
import sys
import urllib.parse
import urllib.request
from collections import defaultdict

import numpy as np
import pandas as pd

ZENODO = "https://zenodo.org/records/18681980/files"
REF_NAME = "PlasmoDB-54-Pfalciparum3D7-Genome.fasta"
ENA_PROJECT = "PRJEB96783"
WINDOW_SIZE = 500
WINDOW_STATS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "window_stats.py")


# ---------------------------------------------------------------------------
# Part 1: depth profiles
# ---------------------------------------------------------------------------

def download(url, path):
    subprocess.run(["curl", "-sS", "-L", "--fail", "-o", path, url], check=True)


def find_cram_url(sample):
    """Ask ENA for this sample's Pf8 CRAM (matched on sample name)."""
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


def get_profiles(samples, pf8_dir):
    """Make sure every sample has a window_stats.tsv; return {sample: path}."""
    inputs, outputs = os.path.join(pf8_dir, "input"), os.path.join(pf8_dir, "output")
    os.makedirs(inputs, exist_ok=True)
    os.makedirs(outputs, exist_ok=True)

    # sample metadata and reference genome (downloaded once)
    meta_path = os.path.join(inputs, "Pf8_samples.txt")
    if not os.path.exists(meta_path):
        print("downloading Pf8 sample metadata")
        download(f"{ZENODO}/Pf8_samples.txt", meta_path)
    meta = pd.read_csv(meta_path, sep="\t").set_index("Sample")

    unknown = [s for s in samples if s not in meta.index]
    if unknown:
        sys.exit(f"not Pf8 sample IDs (check the Sample column of {meta_path}): {', '.join(unknown)}")

    ref = os.path.join(inputs, REF_NAME)
    profiles = {}
    for sample in samples:
        stats = os.path.join(outputs, f"{sample}.window_stats.tsv")
        profiles[sample] = stats
        if os.path.exists(stats):
            print(f"{sample}: profile already made")
            continue
        if not meta.loc[sample, "QC pass"]:
            print(f"{sample}: warning, this sample failed Pf8 QC")

        if not os.path.exists(ref):
            print("downloading the reference genome")
            download(f"{ZENODO}/{REF_NAME}", ref)

        url = find_cram_url(sample)
        if url is None:
            sys.exit(f"{sample}: no Pf8 CRAM in ENA (e.g. Indonesian samples have no reads in the release)")
        cram = os.path.join(inputs, f"{sample}.cram")
        print(f"{sample}: downloading {url}")
        download(url, cram)
        download(url + ".crai", cram + ".crai")

        print(f"{sample}: mosdepth")
        prefix = os.path.join(outputs, sample)
        subprocess.run(["mosdepth", "--fasta", ref, prefix, cram], check=True)
        print(f"{sample}: window mean and SD")
        subprocess.run([sys.executable, WINDOW_STATS, f"{prefix}.per-base.bed.gz", stats, str(WINDOW_SIZE)],
                       check=True)

        os.remove(cram)
        os.remove(cram + ".crai")
    return profiles, meta


def load_profile(path):
    """{Pf8 chromosome name: (mean array, SD array)}, indexed by window number."""
    df = pd.read_csv(path, sep="\t")
    return {chrom: (g["mean"].to_numpy(), g["sd"].to_numpy()) for chrom, g in df.groupby("chrom", sort=False)}


# ---------------------------------------------------------------------------
# Part 2: host VCF
# ---------------------------------------------------------------------------

def draw_depth(rng, mean, sd):
    """One read count from a window with this mean and SD of depth."""
    if mean <= 0:
        return 0
    variance = sd * sd
    if variance <= mean:
        return int(rng.poisson(mean))
    r = mean * mean / (variance - mean)
    return int(rng.negative_binomial(r, r / (r + mean)))


def host_field(counts):
    """GT:AD:DP from one read count per allele (same GT rule as krsim's host VCF)."""
    present = [a for a, n in enumerate(counts) if n > 0]
    if not present:
        gt = "./."
    elif len(present) == 1:
        gt = f"{present[0]}/{present[0]}"
    else:
        # the two alleles with the most reads (ties go to the lower allele), written in allele order
        top_two = sorted(sorted(present, key=lambda a: -counts[a])[:2])
        gt = f"{top_two[0]}/{top_two[1]}"
    return f"{gt}:{','.join(str(n) for n in counts)}:{sum(counts)}"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--strain_vcf", required=True, help="strain VCF from `krsim vcf-strain`")
    parser.add_argument("--ind_node_map_csv", required=True, help="host-strain map from `krsim simulate`")
    parser.add_argument("--samples", required=True, nargs="+",
                        help="Pf8 sample IDs, one per host (1st sample -> 1st host)")
    parser.add_argument("--pf8_dir", required=True, help="folder for Pf8 data (input/) and profiles (output/)")
    parser.add_argument("--random_seed", required=True, type=int, help="random seed (also the simulation ID)")
    parser.add_argument("--output_dir", required=True, help="output folder")
    args = parser.parse_args()

    sim_id = args.random_seed
    rng = np.random.default_rng(args.random_seed)
    os.makedirs(args.output_dir, exist_ok=True)

    # --- hosts and their strains: node_id -> host, VCF column -> node ---
    node_to_ind = {}
    with open(args.ind_node_map_csv, newline="") as f:
        for row in csv.DictReader(f):
            if int(row["sim_id"]) == sim_id:
                node_to_ind[int(row["node_id"])] = int(row["ind_id"])
    if not node_to_ind:
        sys.exit(f"no hosts for seed {sim_id} in {args.ind_node_map_csv}")

    meta_lines, header_cols = [], None
    with open(args.strain_vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header_cols = line.strip().split("\t")
                break
            meta_lines.append(line)
    if header_cols is None:
        sys.exit(f"no #CHROM header line in {args.strain_vcf}")

    ind_to_cols = defaultdict(list)  # host -> its strains' VCF columns
    for col_idx, name in enumerate(header_cols[9:]):
        ind_id = node_to_ind.get(int(name))
        if ind_id is not None:
            ind_to_cols[ind_id].append(col_idx)
    ind_ids = sorted(ind_to_cols)

    # --- one Pf8 sample per host, in order ---
    if len(args.samples) != len(ind_ids):
        sys.exit(f"{len(ind_ids)} hosts in this simulation, so give exactly {len(ind_ids)} Pf8 samples "
                 f"(one per host); got {len(args.samples)}")

    profile_paths, meta = get_profiles(args.samples, args.pf8_dir)
    host_sample = dict(zip(ind_ids, args.samples))
    profiles = {s: load_profile(profile_paths[s]) for s in set(args.samples)}

    # --- strain mix, once per host ---
    mix = {ind_id: rng.dirichlet(np.ones(len(ind_to_cols[ind_id]))) for ind_id in ind_ids}

    # --- host VCF: one depth draw per SNP per host, split between strains ---
    out_vcf = os.path.join(args.output_dir, f"output_{sim_id}_host_depth.vcf")
    with open(out_vcf, "w") as fout:
        for line in meta_lines:
            if line.startswith("##source="):
                fout.write(line)
                fout.write("##source=vcf_host_depth.py (read depth from MalariaGEN Pf8 samples)\n")
                for ind_id in ind_ids:
                    fout.write(f"##depth_profile=<host={ind_id},sample={host_sample[ind_id]}>\n")
            elif line.startswith("##FORMAT=<ID=GT"):
                fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Host genotype from reads: '
                           '0/0=all ref reads, 1/1=all alt reads, 0/1=mixed reads; '
                           'at multi-allelic sites, the two alleles with the most reads">\n')
                fout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description='
                           '"Allelic depths for the ref and each alt allele">\n')
                fout.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description='
                           '"Read depth, drawn from the host\'s Pf8 depth profile">\n')
            else:
                fout.write(line)
        fout.write("\t".join(header_cols[:9] + [str(i) for i in ind_ids]) + "\n")

        with open(args.strain_vcf) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                # krsim names chromosomes 1-14; Pf8 names them Pf3D7_01_v3 ...; VCF positions start at 1
                chrom = f"Pf3D7_{int(cols[0]):02d}_v3"
                window = (int(cols[1]) - 1) // WINDOW_SIZE
                n_alleles = 1 + len(cols[4].split(","))

                fields = []
                for ind_id in ind_ids:
                    means, sds = profiles[host_sample[ind_id]][chrom]
                    depth = draw_depth(rng, means[window], sds[window])
                    strain_reads = rng.multinomial(depth, mix[ind_id])
                    counts = [0] * n_alleles
                    for reads, col_idx in zip(strain_reads, ind_to_cols[ind_id]):
                        allele = cols[9 + col_idx]
                        if allele.isdigit() and int(allele) < n_alleles:
                            counts[int(allele)] += int(reads)
                    fields.append(host_field(counts))

                cols[8] = "GT:AD:DP"
                fout.write("\t".join(cols[:9] + fields) + "\n")

    # --- which Pf8 sample each host used ---
    out_profiles = os.path.join(args.output_dir, f"host_profiles_{sim_id}.csv")
    with open(out_profiles, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sim_id", "ind_id", "pf8_sample", "sample_type", "country", "year"])
        for ind_id in ind_ids:
            s = host_sample[ind_id]
            w.writerow([sim_id, ind_id, s, meta.loc[s, "Sample type"], meta.loc[s, "Country"], int(meta.loc[s, "Year"])])

    # --- each host's strain mix ---
    out_props = os.path.join(args.output_dir, f"host_proportions_{sim_id}.csv")
    with open(out_props, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sim_id", "ind_id", "node_ids", "proportions"])
        for ind_id in ind_ids:
            nodes = [int(header_cols[9 + c]) for c in ind_to_cols[ind_id]]
            w.writerow([sim_id, ind_id, str(nodes), str([round(float(p), 6) for p in mix[ind_id]])])

    print(f"wrote {out_vcf}\nwrote {out_profiles}\nwrote {out_props}")


if __name__ == "__main__":
    main()
