import os
import logging

log = logging.getLogger(__name__)

# -----------------------------
# CHROMOSOME BOUNDARIES
# -----------------------------
chr_ends = [
    657463, 1627369, 2716924, 3903212, 5308593, 6707872,
    8126958, 9592986, 11123805, 12825256, 14864994,
    17201217, 20178769, 23568951
]
chr_starts = [1] + [x + 1 for x in chr_ends[:-1]]
chr_lengths = [end - start + 1 for start, end in zip(chr_starts, chr_ends)]
N_CHR = len(chr_ends)

# -----------------------------
# CHR + LOCAL POSITION
# -----------------------------
def assign_chr_and_pos(pos):
    prev_end = 0
    for i, end in enumerate(chr_ends):
        if pos <= end:
            return str(i + 1), pos - prev_end
        prev_end = end
    return str(N_CHR), pos - chr_ends[-1]

# -----------------------------
# MAIN FUNCTION
# -----------------------------
def process_vcf(ts, input_vcf, output_vcf):
    """
    Takes a raw tskit VCF and:
      - renames sample columns to node IDs
      - assigns correct chromosome from global position
      - resets position to be relative to chromosome start
      - fixes contig header lines
      - deletes the raw input VCF on success
    """

    # verify node order matches VCF column order
    sample_nodes = list(ts.samples())
    log.debug(f"sample nodes (VCF column order): {sample_nodes}")

    # read raw VCF
    with open(input_vcf) as fin:
        lines = fin.readlines()

    # parse #CHROM header line
    for line in lines:
        if line.startswith("#CHROM"):
            header = line.strip().split("\t")
            break

    # log mapping of original column names to node IDs
    original_samples = header[9:]
    for orig, node_id in zip(original_samples, sample_nodes):
        log.debug(f"  column {orig} -> node_{node_id}")

    # build new sample column names from node IDs
    new_samples = [f"{node_id}" for node_id in sample_nodes]
    new_header = header[:9] + new_samples

    # write output VCF
    with open(output_vcf, "w") as fout:

        # write meta lines, replacing old contig lines
        for line in lines:
            if line.startswith("##contig="):
                continue        # drop old contig lines
            if line.startswith("#CHROM"):
                break
            fout.write(line)

        # write one correct contig line per chromosome
        for i, length in enumerate(chr_lengths):
            fout.write(f"##contig=<ID={i+1},length={length}>\n")

        # write column header
        fout.write("\t".join(new_header) + "\n")

        # write variants
        for line in lines:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            pos = int(cols[1])

            # assign chromosome and reset position
            chrom, pos_local = assign_chr_and_pos(pos)
            cols[0] = chrom
            cols[1] = str(pos_local)

            # expand haplotypes
            expanded = []
            for gt in cols[9:]:
                if gt in {".", "./.", ".|."}:
                    expanded.append(".")
                else:
                    expanded.extend(gt.split("|"))

            fout.write("\t".join(cols[:9] + expanded) + "\n")

    # confirm output and clean up
    if os.path.exists(output_vcf) and os.path.getsize(output_vcf) > 0:
        os.remove(input_vcf)
        log.debug(f"deleted raw VCF: {input_vcf}")
    else:
        raise RuntimeError(f"strain VCF not created properly; raw VCF kept: {input_vcf}")