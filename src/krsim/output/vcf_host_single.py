import os
import csv
import logging

log = logging.getLogger(__name__)
import numpy as np
from collections import defaultdict

def process_host_vcf(sim_id, input_vcf, ind_node_map_csv, output_vcf, depth, random_seed=None, proportions_out=None):
    """
    Aggregates a node-level VCF to individual (host) level using a
    Dirichlet-Multinomial model to simulate within-host sequencing.

    For each host:
      1. Strain proportions are sampled once from Dirichlet(alpha=[1, ..., 1])
         where the length equals the number of strains in that host (K or 2K).
      2. Read counts are sampled once from Multinomial(depth, proportions).
         Each read comes from one strain and carries that strain's allele.

    At each variant site the reads are tallied by allele:
      0/0  all `depth` reads carry the ref allele
      1/1  all `depth` reads carry the alt allele
      0/1  reads carry both alleles (at least one ref and one alt read)
      ./.  no reads assigned (should not occur with depth > 0)

    Each sample column carries GT:AD:DP, where AD is "ref_reads,alt_reads"
    and DP is their sum.
    """
    rng = np.random.default_rng(random_seed)

    # --- build node_id -> ind_id from CSV for this sim_id ---
    node_to_ind = {}
    with open(ind_node_map_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row["sim_id"]) == sim_id:
                node_to_ind[int(row["node_id"])] = int(row["ind_id"])

    if not node_to_ind:
        raise ValueError(f"no hosts for seed {sim_id} in {ind_node_map_csv}")

    # --- read node-level VCF ---
    with open(input_vcf) as f:
        lines = f.readlines()

    # find #CHROM header and map column index -> node_id
    header_cols = None
    for line in lines:
        if line.startswith("#CHROM"):
            header_cols = line.strip().split("\t")
            break

    if header_cols is None:
        raise ValueError(f"no #CHROM header line in {input_vcf}")

    # col_idx (0-based within samples) -> node_id
    col_to_node = {idx: int(name) for idx, name in enumerate(header_cols[9:])}

    # per-individual: ordered list of VCF column indices (one per strain)
    ind_to_cols = defaultdict(list)
    for col_idx, node_id in col_to_node.items():
        ind_id = node_to_ind.get(node_id)
        if ind_id is not None:
            ind_to_cols[ind_id].append(col_idx)

    ind_ids = sorted(ind_to_cols.keys())
    log.debug(f"{len(col_to_node)} strains -> {len(ind_ids)} hosts | depth={depth}")

    # --- sample Dirichlet proportions and Multinomial read counts once per host ---
    # reads[ind_id][i] = number of reads from strain i (the i-th column in ind_to_cols[ind_id])
    reads = {}
    proportions_store = {}  # ind_id -> (node_ids, proportions array, reads array)
    for ind_id in ind_ids:
        strain_cols = ind_to_cols[ind_id]
        n_strains = len(strain_cols)
        node_ids = [col_to_node[c] for c in strain_cols]
        proportions = rng.dirichlet(np.ones(n_strains))
        reads[ind_id] = rng.multinomial(depth, proportions)
        proportions_store[ind_id] = (node_ids, proportions, reads[ind_id])
        log.debug(
            f"  ind {ind_id}: {n_strains} strains | "
            f"proportions={np.round(proportions, 3)} | reads={reads[ind_id]}"
        )

    # --- write proportions CSV if requested ---
    if proportions_out is not None:
        with open(proportions_out, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["sim_id", "ind_id", "node_ids", "proportions", "reads"])
            for ind_id in ind_ids:
                node_ids, proportions, read_counts = proportions_store[ind_id]
                writer.writerow([
                    sim_id,
                    ind_id,
                    str(list(node_ids)),
                    str([round(float(p), 6) for p in proportions]),
                    str(list(int(r) for r in read_counts)),
                ])

    # --- write host-level VCF ---
    new_header = header_cols[:9] + [str(ind_id) for ind_id in ind_ids]

    with open(output_vcf, "w") as fout:
        for line in lines:
            if line.startswith("#CHROM"):
                break
            if line.startswith("##FORMAT=<ID=GT"):
                fout.write(
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description='
                    '"Host genotype from Dirichlet-Multinomial read sampling: '
                    '0/0=all ref reads, 1/1=all alt reads, 0/1=mixed reads">\n'
                )
                fout.write(
                    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description='
                    '"Allelic depths for the ref and alt alleles">\n'
                )
                fout.write(
                    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description='
                    '"Read depth">\n'
                )
            else:
                fout.write(line)

        fout.write("\t".join(new_header) + "\n")

        for line in lines:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            sample_fields = []
            for ind_id in ind_ids:
                strain_cols = ind_to_cols[ind_id]
                strain_reads = reads[ind_id]

                alt_reads = 0
                ref_reads = 0
                for i, col_idx in enumerate(strain_cols):
                    allele = cols[9 + col_idx]
                    if allele == "1":
                        alt_reads += strain_reads[i]
                    elif allele == "0":
                        ref_reads += strain_reads[i]

                if alt_reads == 0 and ref_reads == 0:
                    gt = "./."
                elif alt_reads == 0:
                    gt = "0/0"
                elif ref_reads == 0:
                    gt = "1/1"
                else:
                    gt = "0/1"

                dp = ref_reads + alt_reads
                sample_fields.append(f"{gt}:{ref_reads},{alt_reads}:{dp}")

            cols[8] = "GT:AD:DP"
            fout.write("\t".join(cols[:9] + sample_fields) + "\n")

    if not (os.path.exists(output_vcf) and os.path.getsize(output_vcf) > 0):
        raise RuntimeError(f"host VCF not created properly: {output_vcf}")
