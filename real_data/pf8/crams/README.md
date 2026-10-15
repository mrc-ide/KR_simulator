# CRAM files (not in git)

The aligned reads for each sample are CRAM files of 0.5–4 GB, too large for GitHub (100 MB limit per file). They are stored permanently in ENA, in MalariaGEN's Pf8 project [PRJEB96783](https://www.ebi.ac.uk/ena/browser/view/PRJEB96783), one per sample, named by the Pf8 sample ID. Download them into this folder.

Each sample needs two files: `<sample>.cram` (the reads) and `<sample>.cram.crai` (its index). Keep them together in this folder.

## Samples analysed so far

`cram_urls.tsv` lists the CRAMs behind the results in `../mosdepth/`:

| Sample | Sample type | ENA run | CRAM size |
|---|---|---|---:|
| FP0008-C | gDNA | ERR15606936 | 467 MB |

To download them all (run from `real_data/pf8`):

```bash
tail -n +2 crams/cram_urls.tsv | while IFS=$'\t' read sample type run size cram crai; do
    curl -L -o "crams/$sample.cram" "$cram"
    curl -L -o "crams/$sample.cram.crai" "$crai"
done
```

## Any other sample

**Easiest:** the script finds, downloads and runs mosdepth in one go:

```bash
bash scripts/run_mosdepth.sh FP0010-CW
```

**By hand:**

1. Find the sample's CRAM in ENA (replace `FP0010-CW` with any `Sample` ID from `../metadata/Pf8_samples.txt`):

   ```bash
   curl -s -G "https://www.ebi.ac.uk/ena/portal/api/search" \
     --data-urlencode 'result=read_run' \
     --data-urlencode 'query=study_accession=PRJEB96783 AND sample_alias="FP0010-CW"' \
     --data-urlencode 'fields=run_accession,sample_alias,submitted_bytes,submitted_ftp' \
     --data-urlencode 'format=tsv'
   ```

   This prints the run, the file sizes in bytes, and the two file paths (`.cram` and `.cram.crai`), separated by `;`.

2. Download both, adding `https://` to the front of each path (run from `real_data/pf8`):

   ```bash
   curl -L -o crams/FP0010-CW.cram      "https://ftp.sra.ebi.ac.uk/vol1/run/ERR156/ERR15606937/FP0010-CW.cram"
   curl -L -o crams/FP0010-CW.cram.crai "https://ftp.sra.ebi.ac.uk/vol1/run/ERR156/ERR15606937/FP0010-CW.cram.crai"
   ```

## Without downloading

samtools can read a region straight from ENA, which is handy for a quick look:

```bash
samtools view -T reference/PlasmoDB-54-Pfalciparum3D7-Genome.fasta \
  "https://ftp.sra.ebi.ac.uk/vol1/run/ERR156/ERR15606936/FP0008-C.cram" Pf3D7_07_v3:403000-404000 | head
```

## Notes

- Use the sample name (`Sample` column), not the `ENA` column of `Pf8_samples.txt`: that column points to the original sequencing deposit, not the Pf8-aligned CRAM.
- Every Pf8 sample has exactly one CRAM in PRJEB96783 (33,190 of 33,325 samples; 135 have none).
