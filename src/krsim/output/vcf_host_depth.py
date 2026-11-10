import os
import csv
import logging

log = logging.getLogger(__name__)
import numpy as np
from collections import defaultdict

from krsim.output.pf8_profiles import WINDOW_SIZE, load_zip_profile, zip_samples
from krsim.output import depth_model as dm


def draw_depths(rng, means, sds):
    """One read count per host, from each host's window mean and SD of depth: negative binomial with that
    mean and variance; Poisson(mean) if the variance is not above the mean; 0 if the mean is 0."""
    depth = np.zeros(len(means), dtype=np.int64)
    variance = sds * sds
    pois = (means > 0) & (variance <= means)
    nb = (means > 0) & (variance > means)
    depth[pois] = rng.poisson(means[pois])
    m = means[nb]
    r = m * m / (variance[nb] - m)
    depth[nb] = rng.negative_binomial(r, r / (r + m))
    return depth


def site_alleles(genotypes, n_strains):
    """Each strain's allele at one site, from the tab-separated genotype text of the strain VCF line,
    as an int array (-1 where it isn't a plain allele number)."""
    raw = np.frombuffer(genotypes.encode(), dtype=np.uint8)
    if len(raw) == 2 * n_strains - 1:  # the usual case: every allele is one digit
        digits = raw[::2].astype(np.int64) - 48
        if digits.min() >= 0 and digits.max() <= 9:
            return digits
    return np.array([int(a) if a.isdigit() else -1 for a in genotypes.split("\t")], dtype=np.int64)


def host_fields(counts):
    """GT:AD:DP for every host, from its read count per allele (hosts x alleles)."""
    if counts.shape[1] == 2:  # biallelic: the usual case, done for all hosts at once
        ref, alt = counts[:, 0].tolist(), counts[:, 1].tolist()
        gt = ["./.", "1/1", "0/0", "0/1"]  # indexed by 2 * (ref > 0) + (alt > 0)
        return [f"{gt[2 * (r > 0) + (a > 0)]}:{r},{a}:{r + a}" for r, a in zip(ref, alt)]
    return [host_field(c) for c in counts.tolist()]


def host_field(counts):
    """GT:AD:DP from one read count per allele."""
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


def check_error_rates(error_wt_mut, error_mut_wt):
    """Both rates must be numbers from 0 to 1."""
    for name, rate in (("--error_wt_mut", error_wt_mut), ("--error_mut_wt", error_mut_wt)):
        if not 0.0 <= rate <= 1.0:
            raise ValueError(f"{name} must be from 0 to 1 (a chance), got {rate}")


def apply_errors(rng, alleles, use, n_alleles, error_wt_mut, error_mut_wt):
    """
    Sample-level error, in place: each strain's allele in its host (hosts x strain slots, `use` marks real
    alleles) is read wrong with a chance that depends on its true allele. WT (allele 0) is read as MUT with
    chance error_wt_mut, as one of the site's ALT alleles picked with equal chance; MUT (any ALT) is read as
    WT with chance error_mut_wt. Returns the number of (WT -> MUT, MUT -> WT) flips.
    """
    u = rng.random(alleles.shape)
    wt_mut = use & (alleles == 0) & (u < error_wt_mut)
    mut_wt = use & (alleles > 0) & (u < error_mut_wt)
    n_wt_mut, n_mut_wt = int(wt_mut.sum()), int(mut_wt.sum())
    alleles[mut_wt] = 0
    if n_wt_mut:
        alleles[wt_mut] = 1 if n_alleles == 2 else rng.integers(1, n_alleles, n_wt_mut)
    return n_wt_mut, n_mut_wt


def process_host_depth_vcf(sim_id, input_vcf, ind_node_map_csv, output_dir, random_seed=None,
                           depth_zip=None, depth_model=None, error_wt_mut=0.0, error_mut_wt=0.0):
    """
    Aggregates a node-level VCF to host level with realistic read depth. The
    depth comes from one of two sources; give exactly one:

    depth_zip: real depth profiles, from a zip made by `krsim depth-zip`. A
      profile is the mean and SD of read depth in every 500 bp window, in the 3D7
      reference used by MalariaGEN Pf8. Each host gets one profile; if the zip has
      more profiles than hosts, they are picked at random (set by random_seed),
      each used at most once. At each site the depth is drawn from the window's
      mean m and SD s: negative binomial (r = m^2 / (s^2 - m), p = r / (r + m));
      Poisson(m) if s^2 <= m; 0 if m = 0.

    depth_model: settings from depth_model.model_settings (sWGA AR(8) t or gDNA
      ARMA(2,2) t). Each host gets a typical depth drawn from a log-normal and its
      own simulated depth track along the genome (see krsim.output.depth_model);
      at each site the depth is Poisson(window depth).

    Then, as for both:
      - strain proportions, once per host, from Dirichlet(alpha=[1, ..., 1])
      - at every variant site, the depth is split between the host's strains:
        Multinomial(depth, proportions), and reads are tallied by allele -> GT:AD:DP
        (AD has one count per allele; GT lists the two alleles with the most reads;
        ./. if there are no reads)
      - sample-level error (if error_wt_mut or error_mut_wt is above 0): before the reads are
        split, each strain's allele in its host is read wrong with that chance (see apply_errors),
        and all that strain's reads go to the wrong allele. The strain VCF is not changed, so it
        stays the truth. With both rates 0 nothing extra is drawn, so the output is unchanged.

    Writes output_<seed>_host_depth.vcf, host_depth_profiles_<seed>.csv (each
    host's depth source: profile, or typical depth for the model) and
    host_depth_proportions_<seed>.csv. Returns their paths.
    """
    if (depth_zip is None) == (depth_model is None):
        raise ValueError("give exactly one depth source: a depth zip or a depth model")
    check_error_rates(error_wt_mut, error_mut_wt)
    with_errors = error_wt_mut > 0 or error_mut_wt > 0
    n_wt_mut = n_mut_wt = 0
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

    if depth_zip is not None:
        # --- one profile per host, picked at random from the zip (sorted first, so a seed always picks the same) ---
        in_zip = zip_samples(depth_zip)
        if len(in_zip) < len(ind_ids):
            raise ValueError(f"{len(ind_ids)} hosts in this simulation, but {depth_zip} has only {len(in_zip)} "
                             f"depth files; it needs at least one per host")
        samples = [str(s) for s in rng.choice(sorted(in_zip), size=len(ind_ids), replace=False)]
        host_sample = dict(zip(ind_ids, samples))
        profiles = {s: load_zip_profile(depth_zip, in_zip[s]) for s in samples}
        log.debug(f"{len(ind_ids)} hosts, {len(in_zip)} depth files in {depth_zip} | "
                  + ", ".join(f"{i}={s}" for i, s in host_sample.items()))
        # per chromosome: hosts x windows arrays of mean and SD (every profile has at least the 3D7 windows)
        stacked = {}
        for chrom in profiles[samples[0]]:
            n = min(len(profiles[s][chrom][0]) for s in samples)
            stacked[chrom] = (np.vstack([profiles[host_sample[i]][chrom][0][:n] for i in ind_ids]),
                              np.vstack([profiles[host_sample[i]][chrom][1][:n] for i in ind_ids]))

        def site_depths(chrom, window):
            means, sds = stacked[chrom]
            return draw_depths(rng, means[:, window], sds[:, window])

    # --- strain proportions, once per host ---
    proportions = {ind_id: rng.dirichlet(np.ones(len(ind_to_cols[ind_id]))) for ind_id in ind_ids}

    if depth_model is not None:
        # --- modelled depth: a typical depth per host, then its own track along the genome, kept only at the
        # windows that hold a variant site (the whole track is simulated, so its correlation along the genome holds)
        typical = dict(zip(ind_ids, dm.host_depths(rng, depth_model, len(ind_ids))))
        site_windows = defaultdict(set)
        with open(input_vcf) as fin:
            for line in fin:
                if not line.startswith("#"):
                    chrom, pos = line.split("\t", 2)[:2]
                    site_windows[f"Pf3D7_{int(chrom):02d}_v3"].add((int(pos) - 1) // WINDOW_SIZE)
        site_windows = {c: np.array(sorted(w)) for c, w in site_windows.items()}
        where = {c: {int(w): i for i, w in enumerate(ws)} for c, ws in site_windows.items()}
        # per chromosome: hosts x (windows with a site) array of window depth
        kept = {c: np.empty((len(ind_ids), len(ws))) for c, ws in site_windows.items()}
        for h, ind_id in enumerate(ind_ids):
            track = dm.relative_track(rng, depth_model)
            for c, ws in site_windows.items():
                kept[c][h] = typical[ind_id] * track[c][ws]
        log.debug(f"{len(ind_ids)} hosts, modelled depth ({depth_model['model']}) | typical depth: "
                  + ", ".join(f"{i}={typical[i]:.1f}" for i in ind_ids))

        def site_depths(chrom, window):
            return rng.poisson(kept[chrom][:, where[chrom][window]])

    # all hosts at once: each host's strain columns and proportions, padded to the most strains any host has
    n_hosts = len(ind_ids)
    n_slots = max(len(ind_to_cols[i]) for i in ind_ids)
    slot_cols = np.full((n_hosts, n_slots), -1, dtype=np.int64)
    slot_props = np.zeros((n_hosts, n_slots))
    for h, ind_id in enumerate(ind_ids):
        k = len(ind_to_cols[ind_id])
        slot_cols[h, :k] = ind_to_cols[ind_id]
        slot_props[h, :k] = proportions[ind_id]
    is_slot = slot_cols >= 0
    n_strains = len(header_cols) - 9
    slot_host = np.repeat(np.arange(n_hosts), n_slots).reshape(n_hosts, n_slots)

    # --- write host-level VCF: one depth draw per site per host, split between its strains ---
    out_vcf = os.path.join(output_dir, f"output_{sim_id}_host_depth.vcf")
    with open(out_vcf, "w") as fout:
        for line in meta_lines:
            if line.startswith("##source="):
                fout.write(line)
                if depth_zip is not None:
                    fout.write("##source=krsim vcf-host-depth (read depth from real depth profiles)\n")
                    fout.write(f"##depth_zip={os.path.basename(depth_zip)}\n")
                    for ind_id in ind_ids:
                        fout.write(f"##depth_profile=<host={ind_id},sample={host_sample[ind_id]}>\n")
                else:
                    fout.write("##source=krsim vcf-host-depth (modelled read depth)\n")
                    fout.write(f"##depth_model=<{dm.describe(depth_model)}>\n")
                    for ind_id in ind_ids:
                        fout.write(f"##depth_host=<host={ind_id},typical_depth={typical[ind_id]:.3f}>\n")
                fout.write(f"##error_rates=<wt_mut={error_wt_mut:g},mut_wt={error_mut_wt:g},Description=\"Chance "
                           "that a strain's allele in a host is read wrong: WT as MUT (wt_mut), MUT as WT (mut_wt); "
                           "the strain VCF holds the true alleles\">\n")
            elif line.startswith("##FORMAT=<ID=GT"):
                fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Host genotype from reads: '
                           '0/0=all ref reads, 1/1=all alt reads, 0/1=mixed reads; '
                           'at multi-allelic sites, the two alleles with the most reads">\n')
                fout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description='
                           '"Allelic depths for the ref and each alt allele">\n')
                fout.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description='
                           '"Read depth, drawn from the host\'s ' + ("depth profile" if depth_zip is not None else "modelled depth")
                           + '">\n')
            else:
                fout.write(line)
        fout.write("\t".join(header_cols[:9] + [str(i) for i in ind_ids]) + "\n")

        with open(input_vcf) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                # krsim names chromosomes 1-14, profiles use the 3D7 names Pf3D7_01_v3 ...; VCF positions start at 1
                chrom = f"Pf3D7_{int(cols[0]):02d}_v3"
                window = (int(cols[1]) - 1) // WINDOW_SIZE
                n_alleles = 1 + len(cols[4].split(","))

                # every host at once: its depth, split between its strains, then reads tallied by allele
                strain_reads = rng.multinomial(site_depths(chrom, window), slot_props)
                alleles = site_alleles(line.rstrip("\n").split("\t", 9)[9], n_strains)[np.maximum(slot_cols, 0)]
                use = is_slot & (alleles >= 0) & (alleles < n_alleles)
                if with_errors:
                    flips = apply_errors(rng, alleles, use, n_alleles, error_wt_mut, error_mut_wt)
                    n_wt_mut += flips[0]
                    n_mut_wt += flips[1]
                counts = np.zeros((n_hosts, n_alleles), dtype=np.int64)
                np.add.at(counts, (slot_host[use], alleles[use]), strain_reads[use])
                sample_fields = host_fields(counts)

                cols[8] = "GT:AD:DP"
                fout.write("\t".join(cols[:9] + sample_fields) + "\n")

    if not (os.path.exists(out_vcf) and os.path.getsize(out_vcf) > 0):
        raise RuntimeError(f"host VCF not created properly: {out_vcf}")
    if with_errors:
        log.info(f"sample-level error: {n_wt_mut:,} strain alleles read WT -> MUT, {n_mut_wt:,} MUT -> WT "
                 f"(rates {error_wt_mut:g}, {error_mut_wt:g})")

    # --- each host's depth source: which profile, or its typical depth under the model ---
    out_profiles = os.path.join(output_dir, f"host_depth_profiles_{sim_id}.csv")
    with open(out_profiles, "w", newline="") as f:
        writer = csv.writer(f)
        if depth_zip is not None:
            writer.writerow(["sim_id", "ind_id", "sample", "depth_file"])
            for ind_id in ind_ids:
                s = host_sample[ind_id]
                writer.writerow([sim_id, ind_id, s, in_zip[s]])
        else:
            writer.writerow(["sim_id", "ind_id", "depth_model", "typical_depth"])
            for ind_id in ind_ids:
                writer.writerow([sim_id, ind_id, depth_model["model"], round(float(typical[ind_id]), 3)])

    # --- each host's strain proportions ---
    out_props = os.path.join(output_dir, f"host_depth_proportions_{sim_id}.csv")
    with open(out_props, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sim_id", "ind_id", "node_ids", "proportions"])
        for ind_id in ind_ids:
            node_ids = [int(header_cols[9 + c]) for c in ind_to_cols[ind_id]]
            writer.writerow([sim_id, ind_id, str(node_ids), str([round(float(p), 6) for p in proportions[ind_id]])])

    return out_vcf, out_profiles, out_props
