# Pf8 coverage analysis (mosdepth)

Real *P. falciparum* sequencing data from MalariaGEN's Pf8 release, for looking at read depth, e.g. comparing sWGA and gDNA samples. Work in progress on the `pf8-real-data` branch of KR_simulator; intended to become a krsim feature later.

- Pf8 paper: Wellcome Open Research 10:325 (2025), https://doi.org/10.12688/wellcomeopenres.24031.1
- Official page: https://malariagen.net/resource/36/

## Folders

| Folder | Contents | Source |
|---|---|---|
| `metadata/` | `Pf8_samples.txt`: one row per sample (`Sample` ID, `Sample type`, `QC pass`, country, year...) | https://zenodo.org/records/18681980/files/Pf8_samples.txt |
| `reference/` | PlasmoDB-54 *P. falciparum* 3D7 genome + `.fai` index (the reference the Pf8 reads are aligned to) | https://zenodo.org/records/18681980/files/PlasmoDB-54-Pfalciparum3D7-Genome.fasta |
| `crams/` | One `<sample>.cram` + `.crai` per sample: reads aligned by MalariaGEN for Pf8. **Not in git** (too large): download them with [`crams/README.md`](crams/README.md) | ENA project PRJEB96783, matched on sample name |
| `mosdepth/<sample>/` | mosdepth results for each sample | made by `scripts/run_mosdepth.sh` |
| `scripts/` | `run_mosdepth.sh` | |

## Run

```bash
conda activate mosdepth
cd KR_simulator/real_data/pf8
bash scripts/run_mosdepth.sh FP0008-C
```

For each sample it finds the CRAM in ENA, downloads it into `crams/` (skipped if already there), and runs:

```bash
mosdepth --fasta reference/PlasmoDB-54-Pfalciparum3D7-Genome.fasta -n -x --by 1000 --thresholds 1,5,10 \
    mosdepth/<sample>/<sample> crams/<sample>.cram
```

- Several samples at once: `bash scripts/run_mosdepth.sh FP0008-C FP0010-CW`
- Delete each CRAM afterwards to save space (0.5–4 GB each): `DELETE_CRAM=1 bash scripts/run_mosdepth.sh FP0008-C`

## On another computer

```bash
git clone https://github.com/mrc-ide/KR_simulator.git   # or, in an existing clone: git fetch
cd KR_simulator
git switch pf8-real-data
cd real_data/pf8
```

The metadata, reference genome, scripts and mosdepth results all come with the repo. The CRAMs are not in git (0.5–4 GB each, over GitHub's 100 MB file limit). To get them, see **[`crams/README.md`](crams/README.md)**: it lists the download links for the samples analysed so far, and how to fetch any other sample. `scripts/run_mosdepth.sh` also downloads a sample's CRAM when you run that sample.

## mosdepth output (in `mosdepth/<sample>/`)

| File | Contents |
|---|---|
| `<sample>.mosdepth.summary.txt` | **Start here.** Mean depth per chromosome; the `total` row is the whole genome |
| `<sample>.regions.bed.gz` | Mean depth in each 1 kb window along the genome |
| `<sample>.thresholds.bed.gz` | For each 1 kb window, bases with at least 1, 5 and 10 reads |
| `<sample>.mosdepth.global.dist.txt` | Fraction of the genome with at least 0, 1, 2... reads |
| `<sample>.mosdepth.region.dist.txt` | The same, within the 1 kb windows |

Look at them with `cat`, or `gzcat file.bed.gz | head` for the `.gz` files.

Note: `mosdepth/FP0008-C/` holds an earlier run made by hand with prefix `FP0008` (without `--thresholds`). Running the script for FP0008-C adds files named `FP0008-C.*` alongside them.
