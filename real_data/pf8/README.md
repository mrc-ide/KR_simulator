# Pf8 read depth (mosdepth)

Read depth for real *P. falciparum* samples from MalariaGEN's Pf8 release (Wellcome Open Research 10:325, 2025, https://doi.org/10.12688/wellcomeopenres.24031.1).

The repo holds only the script, [`pf8_mosdepth.sh`](pf8_mosdepth.sh). All data goes in a folder outside the repo, set by `PF8=` at the top of the script:

```
pf8/
├── input/    Pf8_samples.txt, <sample>.cram + .crai, reference .fasta + .fai
└── output/   mosdepth results
```

## Run

```bash
conda create -n mosdepth -c conda-forge -c bioconda mosdepth samtools pandas   # once per laptop
conda activate mosdepth
bash real_data/pf8/pf8_mosdepth.sh
```

Steps: 1) sample information, 2) count gDNA and sWGA samples, 3) find a sample's CRAM in ENA, 4) download the CRAM, 5) download the reference, 6) explore the alignment (optional), 7) mosdepth: depth at every base, 8) mean and SD of depth per 500 bp window (`depth_model/window_stats.py`).

## What is in a CRAM

- **One CRAM per sample.** It holds all the sequencing reads for that sample, aligned to the 3D7 reference (PlasmoDB-54). Where a sample was sequenced in several runs, the reads were merged into the one file.
- **A sample is one blood draw from one patient at one time point, not always one patient.** Most patients have one sample, but some were sampled more than once (e.g. before and after treatment). The `All samples same case` column of `Pf8_samples.txt` links samples from the same patient:

  | Samples per patient | Patients |
  |---|---:|
  | 1 | 30,312 |
  | 2 | 1,018 |
  | 3–8 | 278 |

- **Mapped and unmapped reads are both in the file.** For FP0008-C, 88.6% of reads map to 3D7. The unmapped 11% are probably mostly leftover human DNA (not checked).
- **Duplicates and secondary alignments are flagged, not removed.** mosdepth skips them by default.
- **Reads from all parasites in the blood sample are mixed together.** An infection can contain several genetically different clones, and the CRAM does not separate them, so the depth is for all of them combined. `Pf8_fws.tsv` in the Pf8 release gives how mixed each infection is.
- **The reads are stored relative to the reference.** A CRAM keeps only how each read differs from the reference, which is why it is smaller than a BAM. Tools need the reference (`-T` for samtools, `--fasta` for mosdepth) to read it; mosdepth reads the CRAM directly, so no BAM is needed.
