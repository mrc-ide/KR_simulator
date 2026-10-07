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

    # --- read the VCF header (meta lines and #CHROM line); data lines are streamed later ---
    meta_lines = []
    header_cols = None
    with open(input_vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header_cols = line.strip().split("\t")
                break
            meta_lines.append(line)

    if header_cols is None:
        raise ValueError(f"no #CHROM header line in {input_vcf}")

    # col_idx (0-based within samples) -> node_id
    col_to_node = {idx: int(name) for idx, name in enumerate(header_cols[9:])}
    n_cols = len(col_to_node)

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

    # --- precompute, for a block of sites at once, each host's ref and alt reads ---
    # Columns are reordered so each host's strains are next to each other; each
    # strain is weighted by its read count, and each host's block is summed.
    order = np.array([c for ind_id in ind_ids for c in ind_to_cols[ind_id]], dtype=np.int64)
    weights = np.concatenate([reads[ind_id] for ind_id in ind_ids]).astype(np.int64)
    starts = np.cumsum([0] + [len(ind_to_cols[ind_id]) for ind_id in ind_ids[:-1]])
    max_reads = int(depth) + 1

    def host_fields(ref_reads, alt_reads):
        """Format one host's GT:AD:DP field (same rules as the per-strain loop)."""
        if alt_reads == 0 and ref_reads == 0:
            gt = "./."
        elif alt_reads == 0:
            gt = "0/0"
        elif ref_reads == 0:
            gt = "1/1"
        else:
            gt = "0/1"
        return f"{gt}:{ref_reads},{alt_reads}:{ref_reads + alt_reads}"

    field_cache = {}

    def slow_line(cols):
        """Original per-strain loop, used for any line the fast parser can't read."""
        sample_fields = []
        for ind_id in ind_ids:
            alt_reads = 0
            ref_reads = 0
            for i, col_idx in enumerate(ind_to_cols[ind_id]):
                allele = cols[9 + col_idx]
                if allele == "1":
                    alt_reads += reads[ind_id][i]
                elif allele == "0":
                    ref_reads += reads[ind_id][i]
            sample_fields.append(host_fields(ref_reads, alt_reads))
        return sample_fields

    def write_block(fout, block):
        """Fast path: genotypes as a byte array, host read sums with numpy."""
        fixed, genos = [], []
        for line in block:
            first9 = line.strip().split("\t", 9)
            fixed.append(first9[:9])
            genos.append(first9[9] if len(first9) == 10 else "")
        expected_len = 2 * n_cols - 1
        ok = [len(g) == expected_len and g.isascii() for g in genos]
        fast_idx = [i for i, good in enumerate(ok) if good]
        if fast_idx:
            # one byte per strain: single-character alleles separated by tabs
            raw = np.frombuffer("".join(genos[i] for i in fast_idx).encode("ascii"), dtype=np.uint8)
            raw = raw.reshape(len(fast_idx), expected_len)
            # a row whose separators are not all tabs is not one-character alleles: use the slow path
            for row_i in np.flatnonzero(~(raw[:, 1::2] == ord("\t")).all(axis=1)):
                ok[fast_idx[row_i]] = False
            alleles = raw[:, ::2][:, order]
            alt = np.add.reduceat((alleles == ord("1")) * weights, starts, axis=1)
            ref = np.add.reduceat((alleles == ord("0")) * weights, starts, axis=1)
            codes = ref * max_reads + alt
            code_row = dict(zip(fast_idx, range(len(fast_idx))))
        for line_i, (line, head, good) in enumerate(zip(block, fixed, ok)):
            if good:
                row = codes[code_row[line_i]]
                sample_fields = []
                for code in row.tolist():
                    field = field_cache.get(code)
                    if field is None:
                        field = field_cache[code] = host_fields(code // max_reads, code % max_reads)
                    sample_fields.append(field)
            else:
                sample_fields = slow_line(line.strip().split("\t"))
            head[8] = "GT:AD:DP"
            fout.write("\t".join(head + sample_fields) + "\n")

    # --- write host-level VCF ---
    new_header = header_cols[:9] + [str(ind_id) for ind_id in ind_ids]

    with open(output_vcf, "w") as fout:
        for line in meta_lines:
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

        # stream the data lines in blocks, so memory stays small for any run size
        block_size = 2000
        with open(input_vcf) as fin:
            block = []
            for line in fin:
                if line.startswith("#"):
                    continue
                block.append(line)
                if len(block) == block_size:
                    write_block(fout, block)
                    block = []
            if block:
                write_block(fout, block)

    if not (os.path.exists(output_vcf) and os.path.getsize(output_vcf) > 0):
        raise RuntimeError(f"host VCF not created properly: {output_vcf}")
