import os
import csv
import logging

log = logging.getLogger(__name__)
import numpy as np
from collections import defaultdict

from krsim.output.pf8_profiles import WINDOW_SIZE, get_profiles, load_profile


def draw_depth(rng, mean, sd):
    """One read count from a window with this mean and SD of depth."""
    if mean <= 0:
        return 0
    variance = sd * sd
    if variance <= mean:
        return int(rng.poisson(mean))
    # negative binomial with this mean and variance
    r = mean * mean / (variance - mean)
    return int(rng.negative_binomial(r, r / (r + mean)))


def host_field(counts):
    """GT:AD:DP from one read count per allele (the same GT rule as `krsim vcf-host`)."""
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


def process_host_depth_vcf(sim_id, input_vcf, ind_node_map_csv, samples, pf8_dir, output_dir, random_seed=None):
    """
    Aggregates a node-level VCF to host level, with each host's read depth
    following a real MalariaGEN Pf8 sample (a more realistic alternative to
    process_host_vcf, which uses one constant depth).

    Once per host:
      1. it gets one Pf8 sample's depth profile, in order (1st sample -> 1st host);
         a profile is the mean and SD of depth in every 500 bp window
      2. strain proportions are sampled from Dirichlet(alpha=[1, ..., 1])
    At every variant site, for every host:
      3. the host's depth is drawn from its profile window's mean m and SD s:
         negative binomial (r = m^2 / (s^2 - m), p = r / (r + m)); Poisson(m)
         if s^2 <= m; 0 if m = 0
      4. the depth is split between the host's strains: Multinomial(depth, proportions)
      5. reads are tallied by allele -> GT:AD:DP (AD has one count per allele;
         GT lists the two alleles with the most reads; ./. if there are no reads)

    Writes output_<seed>_host_depth.vcf, host_depth_profiles_<seed>.csv (which
    Pf8 sample each host's depth came from) and host_depth_proportions_<seed>.csv.
    Returns their paths.
    """
    rng = np.random.default_rng(random_seed)

    # --- build node_id -> ind_id from CSV for this sim_id ---
    node_to_ind = {}
    with open(ind_node_map_csv, newline="") as f:
        for row in csv.DictReader(f):
            if int(row["sim_id"]) == sim_id:
                node_to_ind[int(row["node_id"])] = int(row["ind_id"])
    if not node_to_ind:
        raise ValueError(f"no hosts for seed {sim_id} in {ind_node_map_csv}")

    # --- read the VCF header (meta lines and #CHROM line); data lines are streamed later ---
    meta_lines, header_cols = [], None
    with open(input_vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header_cols = line.strip().split("\t")
                break
            meta_lines.append(line)
    if header_cols is None:
        raise ValueError(f"no #CHROM header line in {input_vcf}")

    # per host: its strains' VCF column indices (0-based within samples)
    ind_to_cols = defaultdict(list)
    for col_idx, name in enumerate(header_cols[9:]):
        ind_id = node_to_ind.get(int(name))
        if ind_id is not None:
            ind_to_cols[ind_id].append(col_idx)
    ind_ids = sorted(ind_to_cols)

    # --- one Pf8 sample per host, in order ---
    if len(samples) != len(ind_ids):
        raise ValueError(f"{len(ind_ids)} hosts in this simulation, so give exactly {len(ind_ids)} Pf8 samples "
                         f"(one per host); got {len(samples)}")
    if len(set(samples)) != len(samples):
        raise ValueError("each host needs a different Pf8 sample; a sample is listed more than once")

    profile_paths, meta = get_profiles(samples, pf8_dir)
    host_sample = dict(zip(ind_ids, samples))
    profiles = {s: load_profile(profile_paths[s]) for s in samples}
    log.debug(f"{len(ind_ids)} hosts | profiles: " + ", ".join(f"{i}={s}" for i, s in host_sample.items()))

    # --- strain proportions, once per host ---
    proportions = {ind_id: rng.dirichlet(np.ones(len(ind_to_cols[ind_id]))) for ind_id in ind_ids}

    # --- write host-level VCF: one depth draw per site per host, split between its strains ---
    out_vcf = os.path.join(output_dir, f"output_{sim_id}_host_depth.vcf")
    with open(out_vcf, "w") as fout:
        for line in meta_lines:
            if line.startswith("##source="):
                fout.write(line)
                fout.write("##source=krsim vcf-host-depth (read depth from MalariaGEN Pf8 samples)\n")
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

        with open(input_vcf) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                # krsim names chromosomes 1-14, Pf8 Pf3D7_01_v3 ...; VCF positions start at 1
                chrom = f"Pf3D7_{int(cols[0]):02d}_v3"
                window = (int(cols[1]) - 1) // WINDOW_SIZE
                n_alleles = 1 + len(cols[4].split(","))

                sample_fields = []
                for ind_id in ind_ids:
                    means, sds = profiles[host_sample[ind_id]][chrom]
                    depth = draw_depth(rng, means[window], sds[window])
                    strain_reads = rng.multinomial(depth, proportions[ind_id])
                    counts = [0] * n_alleles
                    for reads, col_idx in zip(strain_reads, ind_to_cols[ind_id]):
                        allele = cols[9 + col_idx]
                        if allele.isdigit() and int(allele) < n_alleles:
                            counts[int(allele)] += int(reads)
                    sample_fields.append(host_field(counts))

                cols[8] = "GT:AD:DP"
                fout.write("\t".join(cols[:9] + sample_fields) + "\n")

    if not (os.path.exists(out_vcf) and os.path.getsize(out_vcf) > 0):
        raise RuntimeError(f"host VCF not created properly: {out_vcf}")

    # --- which Pf8 sample each host's depth came from ---
    out_profiles = os.path.join(output_dir, f"host_depth_profiles_{sim_id}.csv")
    with open(out_profiles, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sim_id", "ind_id", "pf8_sample", "sample_type", "country", "year"])
        for ind_id in ind_ids:
            s = host_sample[ind_id]
            writer.writerow([sim_id, ind_id, s, meta.loc[s, "Sample type"], meta.loc[s, "Country"],
                             int(meta.loc[s, "Year"])])

    # --- each host's strain proportions ---
    out_props = os.path.join(output_dir, f"host_depth_proportions_{sim_id}.csv")
    with open(out_props, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sim_id", "ind_id", "node_ids", "proportions"])
        for ind_id in ind_ids:
            node_ids = [int(header_cols[9 + c]) for c in ind_to_cols[ind_id]]
            writer.writerow([sim_id, ind_id, str(node_ids), str([round(float(p), 6) for p in proportions[ind_id]])])

    return out_vcf, out_profiles, out_props
