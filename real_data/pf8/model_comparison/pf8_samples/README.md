# Pf8 read depth: 10 sWGA and 10 gDNA samples

Read depth along the genome for 20 real MalariaGEN Pf8 samples, for building and comparing depth models
(a third way of making the krsim host VCF).

- Samples: 10 sWGA and 10 gDNA, picked at random from QC-pass Pf8 samples (no filtering on file size), plus FP0008-C (gDNA) from earlier
- Source: Pf8 CRAMs from ENA project PRJEB96783 (aligned to PlasmoDB-54 3D7 with bwa-mem, duplicates marked)
- Tools: mosdepth 0.3.14 (default flags: duplicates, secondary, QC-fail and unmapped reads not counted; overlapping mates counted once)
- The CRAMs were deleted after processing; `crams.tsv` has their URLs if one is needed again

```
pf8_samples/
├── samples.tsv            the 20 samples: type, country, year, % callable, study
├── crams.tsv              sample, CRAM size (GB), CRAM URL
├── run_samples.sh         download -> mosdepth -> 500 bp windows -> plots, for every sample in crams.tsv
├── compare.sh             rebuilds compare/ from every finished sample
├── compare_samples.py     mean-vs-SD panels + summary table
├── compare_depth.py       box plot + depth along chr7 / crt per sample type
├── reference/             PlasmoDB-54 3D7 reference (.fasta + .fai)
├── <sample>/              one folder per sample (below)
├── compare/
│   ├── sWGA/              the 10 sWGA samples together
│   ├── gDNA/              the 10 gDNA samples together
│   └── both/              all 20 side by side
└── model/                 exploring the depth structure, and the list of models to test
```

## In each sample folder

| File | What it is |
|---|---|
| `mosdepth/<sample>.per-base.bed.gz` | Depth at every base: `chrom start end depth`. Bases with the same depth are merged into one row |
| `mosdepth/<sample>.per-base.bed.gz.csi` | Index, to read one region with tabix |
| `depth_plots/<sample>.window_stats.tsv` | Mean and SD of depth in each 500 bp window: `chrom start end mean sd` |
| `depth_plots/<sample>_depth.png` | Window depth along all 14 chromosomes (line = mean, band = ± SD) |
| `depth_plots/<sample>_depth_crt.png` | The same, zoomed into chr7:350,000-450,000 around *crt* |
| `depth_plots/<sample>_mean_sd.png` | SD against mean for every window, with the Poisson line (SD = √mean) |

Read one region of the per-base depth:

```bash
~/miniforge3/envs/mosdepth/bin/tabix RCN25280/mosdepth/RCN25280.per-base.bed.gz Pf3D7_07_v3:403000-404000
```

## compare/

| File | In | What it is |
|---|---|---|
| `depth_<type>.png` | sWGA, gDNA | The 10 samples of one type together, one colour each: window mean (line) and ± SD (band), along chr7 and around *crt* |
| `depth_boxplot.png` | all three | Spread of window depth per sample (box = middle 50% of windows, whiskers = 5-95%) |
| `mean_sd_all.png` | all three | One mean-vs-SD panel per sample, same axes |
| `summary.tsv` | all three | Per sample: median depth, fraction of windows under 5x, fraction above Poisson, median variance / mean |

Note: two gDNA samples from study 1195-PF-TRAC2-DONDORP (RCN12993, QC0637-C) have depth as uneven as sWGA.

## model/

| File | What it is |
|---|---|
| `explore_depth.py` | Splits log window depth into sample level + shared window effect + sample-specific residual; GC content, autocorrelation, tails. `python explore_depth.py sWGA` |
| `sWGA_explore.txt`, `sWGA_explore.png` | Its results for sWGA |
| `MODEL_LIST.md` | The depth models to test (M0-M8) and how to compare them |
