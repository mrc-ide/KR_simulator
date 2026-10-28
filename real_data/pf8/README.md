# Pf8 read depth (mosdepth)

Read depth for real *P. falciparum* samples from MalariaGEN's Pf8 release (Wellcome Open Research 10:325, 2025, https://doi.org/10.12688/wellcomeopenres.24031.1).

The repo holds only scripts. All data goes in a folder outside the repo, set by `PF8=` at the top of each script:

- [`extraction/pf8_mosdepth.sh`](extraction/pf8_mosdepth.sh): every step for one sample, to run by hand.
- [`all_samples/`](all_samples/): the same steps for many samples as a SLURM array job, timed, keeping only `<sample>.window_stats.tsv` (see below).
- [`depth_model/`](depth_model/): plots and notebooks that model the window stats.

```
pf8/
├── input/    Pf8_samples.txt, <sample>.cram + .crai, reference .fasta + .fai, sample lists
├── output/   mosdepth results (<sample>.window_stats.tsv), Pf8_sample_info.tsv
└── cram/, logs/, timings/, figures/   used by all_samples/
```

## Run

```bash
conda create -n mosdepth -c conda-forge -c bioconda mosdepth samtools pandas   # once per laptop
conda activate mosdepth
bash real_data/pf8/extraction/pf8_mosdepth.sh
```

Steps: 1) sample information, 2) count gDNA and sWGA samples, 3) find a sample's CRAM in ENA, 4) download the CRAM, 5) download the reference, 6) explore the alignment (optional), 7) mosdepth: depth at every base, 8) mean and SD of depth per 500 bp window (krsim's `krsim.output.pf8_profiles`).

## All samples (SLURM array)

On the cluster, with the `mosdepth` conda env active and from `all_samples/` (`PF8=/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth`; `input/` also needs `Pf8_ena_runs.tsv`, see `pf8_crams.py`):

```bash
python pf8_crams.py      # sample lists in $PF8/input, and $PF8/output/Pf8_sample_info.tsv
mkdir -p $PF8/logs

# test: time each stage on a few samples, then extrapolate to all of Pf8
sbatch --array=1-16%16 pf8_depth_array.sbatch $PF8/input/timing_subset.tsv
python plot_timings.py   # figures in $PF8/figures

# all of Pf8: 4 arrays of at most 10000 tasks, 50 at once, QC-pass samples first
# (at most 20000 tasks can be queued, so submit part 3 once part 1 is done, chained after part 2)
sbatch --array=1-10000%50 pf8_depth_array.sbatch $PF8/input/Pf8_crams_part1.tsv
sbatch --dependency=afterany:<part 1 job ID> --array=1-10000%50 pf8_depth_array.sbatch $PF8/input/Pf8_crams_part2.tsv

# when nothing is queued: every sample still without window stats, and the commands to run them again
python pf8_missing.py
```

Each task downloads one CRAM, runs mosdepth, writes `$PF8/output/<sample>.window_stats.tsv`, then deletes the CRAM and all other mosdepth output. A sample already done is skipped, so any list can be run again. Stage times and peak memory go in `$PF8/timings/<sample>.tsv`; a sample that fails gets a line in `$PF8/failed.tsv` with the stage it failed at.

- **QC.** `$PF8/output/Pf8_sample_info.tsv` has every Pf8 sample's metadata (`QC pass`, `Exclusion reason`, `% callable`, `All samples same case`, ...) and CRAM size, for filtering the window stats later. All samples are extracted, QC pass or not.
- **ENA.** Downloads drop partway or are refused (403, empty reply) when many tasks run at once, so each task starts at a random time in its first minute, resumes the CRAM until its size matches ENA's, and retries with growing waits (10 tries).
- **Time** (tests on 116 samples, 50 at once): about 30 s + 0.6 min per GB of CRAM, ~90% of it downloading; at most 1.3 GB of memory. All of Pf8 (33,190 CRAMs, 76 TB) is about 1,040 task-hours, ~22 h at 50 at once.

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

## Depth models

- [`depth_model/`](depth_model/): plotting scripts and notebooks exploring depth along the genome (autocorrelation, candidate processes).
- [`model_comparison/`](model_comparison/): `compare_processes.py` fits every candidate process (AR, ARMA, Gaussian processes) to every sample and compares them by AIC, BIC and held-out likelihood, separately for sWGA and gDNA.
