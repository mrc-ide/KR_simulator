# Notes for Claude

## Simple depth model for the host VCF (status, 2026-11-03)

Full record: [real_data/pf8/model_comparison/RECORD.md](real_data/pf8/model_comparison/RECORD.md)
(written before the 2026-11-03 work below; its "t wins every structure in every sample" is slightly too strong).

**Chosen models** (log window depth along the genome, fitted to Pf8 samples by `best_simple_model.py`):

| Type | Model | ν (median) first 21 / new 20 | Best in first 21 / new 20 |
|---|---|---|---|
| gDNA | ARMA(2,2) + Student t kicks | ≈ 7.0 / 6.8 | 8 of 11 / 9 of 10 |
| sWGA | AR(8) + Student t kicks | ≈ 3.6 / 4.3 | 8 of 10 / 7 of 10 (ARMA(2,2) t in the other 3) |

- The t kick is location-scale (μ, scale, ν); ν is fitted as `2.05 + exp(θ)`
  (`compare_arma_kicks.py`, `unpack` / `t_logpdf`). ν sets tail weight: small = frequent big jumps, ∞ = Gaussian.
- t beats Gaussian in 146 of 147 fits (first 21) and 140 of 140 (new 20). The one exception is a tie
  (SPT50506 AR(2), ν = 178, i.e. Gaussian). Median held-out gain ≈ 1,200 gDNA; 9,100 / 6,250 sWGA.
- Kick shape mattered more than correlation shape on the median, and in 14 of 21 samples; not in near-Gaussian
  samples (SPT50506; the sWGA-like gDNA samples QC0637-C, RCN12993, and in the new set SPT00902, ν ≈ 16).
- Checks done: held-out likelihood + BIC (model choice); kicks uncorrelated; tail calibration
  (sWGA t is well calibrated; gDNA t still has too many big drops, too few big spikes). Same on the new 20.
- Known misses, left out on purpose: sWGA calm/jumpy stretches (GARCH; visible in figure 7b), gDNA lopsided
  dips (skewed t), the shared genome-wide map of low/high regions.

**Agreed goal (from the user):** the model does not need to be precise or statistically tight. It is a
realistic starting point, better than one depth value across the whole genome. Users who need detail
supply their own depth data (`krsim vcf-host-depth`). Don't push for further model refinement.

**In krsim since 2026-11-09** (branch `depth-model`): `krsim vcf-host-depth --depth_model sWGA|gDNA`; see the section below.

### Where things are

- `real_data/pf8/model_comparison/results/simple_models/`: fits for the first 21 samples (sample list in
  `pf8_samples/samples.tsv`).
- `real_data/pf8/model_comparison/results/simple_models_new20/`: the same for 20 new samples (10 gDNA, 10 sWGA),
  list in its `samples.tsv`. Made by
  `python best_simple_model.py <pf8_samples dir> results/simple_models_new20` (~30 s on a laptop).
- `real_data/pf8/model_comparison/results/kick_gain/`: Gaussian vs t figures, all from `plot_kick_gain.py`:
  1 kick shapes; 2 score gain per sample; 3 tail check; 4 / 4b every model shape × kick type (first 21 / new 20);
  5 simulated Gaussian vs t kicks and depth; 6 how the held-out score is built; 7a / 7b real sample (PD0048-C gDNA,
  SPT43924 sWGA) next to genomes simulated from its own Gaussian and t fits.
- The depth data (`<sample>/depth_plots/<sample>.window_stats.tsv`) is **not in git**. On the user's laptop:
  `~/pf8_data/pf8_samples/` (new 20). The first 21 are on the first laptop only (see RECORD.md).
  On the HPC, Pf8 window stats are written to `/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth/output/`
  by `real_data/pf8/all_samples/` (flat `<sample>.window_stats.tsv`, no `depth_plots/` folder).

### Running on the HPC (Apocrita, repo at `/data/home/bty208/IMPERIAL/MRC_IDE/git/KR_simulator`)

- Python needs pandas, matplotlib, scipy, numba, statsmodels.
- `plot_kick_gain.py`: figures 1–6 only need the tables in git. Figures 7a/7b need the depth data; set
  `PF8_SAMPLES` at the top of the script (they are skipped if the files are missing).
- `best_simple_model.py <pf8_dir> <out>` expects `<pf8_dir>/samples.tsv` (columns `Sample`, `Sample type`) and
  `<pf8_dir>/<sample>/depth_plots/<sample>.window_stats.tsv`; point it at a folder laid out that way.

## User depth zips for the host VCF (status, 2026-11-09, branch `depth-zip`)

`krsim vcf-host-depth --depth_zip FILE` replaces the old Pf8 download route (`--samples`, `--pf8_dir`, ENA
download and mosdepth during the run were removed). The zip holds `<sample>.window_stats.tsv` files (mean and SD per
500 bp window, 3D7 names and lengths, checked by `check_profile` in `src/krsim/output/pf8_profiles.py`). Each host
gets one file, picked at random by the seed without replacement; the zip needs at least Nh files. Users make the zip
with the helper `krsim depth-zip --bed <sample>.per-base.bed.gz ... --output my.zip` (listed last in `krsim --help`
as "Helper (not a pipeline step)"). Example zips of 100 sWGA and 100 gDNA QC-pass Pf8 samples are made by
`real_data/pf8/all_samples/make_example_zips.py` into `/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth/example_zips/`;
each zip also holds its own `samples.tsv` (sample, type, country, year, study), which krsim ignores. They are meant
for a GitHub Release named `example-depth-data-v1`, not git (about 94 MB together); not published yet. Pf8 data are
CC BY 4.0; the docs and release notes cite MalariaGEN et al. (2025), Wellcome Open Research 10:325. The module keeps the name `pf8_profiles` because the `real_data/pf8/` scripts call
`python -m krsim.output.pf8_profiles`.

## Future work

- **Other references / chromosome lengths.** krsim is fixed to the 3D7 chromosome lengths used by Pf8
  (PlasmoDB-54, `COMP_BREAKPOINTS` / `GENOME_LENGTH` in `src/krsim/sim/params.py`), and depth files must match them.
  Later: an option to give other chromosome names and lengths.

## Modelled depth for the host VCF (status, 2026-11-09, branch `depth-model`)

`krsim vcf-host-depth --depth_model sWGA|gDNA` (instead of `--depth_zip`). Code: `src/krsim/output/depth_model.py`;
defaults: `src/krsim/data/depth_model_defaults.json`, made by
`real_data/pf8/model_comparison/make_depth_model_defaults.py` (`reuse` rebuilds from the saved per-sample table
`$PF8/model_fits/depth_model_samples.tsv` in seconds). Design agreed with the user, one part at a time:

1. One shared model per type (sWGA AR(8) t, gDNA ARMA(2,2) t); defaults = Pf8 medians (`shared_params.tsv`).
2. Each host's typical (median) depth drawn once from a log-normal of Pf8 sample medians: centre 63.0x / 76.8x,
   SD of log 0.919 / 0.910 (robust, IQR / 1.349). `--depth_mean`, `--depth_spread` (0 = all hosts the same).
3. Window depth = typical depth x relative track; relative track = exp(y + offset) - 0.05, ends lowered, scaled to
   median 1 over the genome, capped. Site depth = Poisson(window depth). The whole track is simulated
   (2,000-window burn-in per chromosome) but only SNP windows are kept.
4. Cap = 99th percentile over samples of (highest window / median): 14,293x sWGA, 116x gDNA. The model's upper
   tail is far heavier than real (typical simulated sWGA host max ~1,900x vs real 175x), so the cap binds often,
   especially for gDNA.
5. sWGA end zone: 5 kb end bins with median depth / sample median < 0.1 -> 30 kb at 0.0035x. Real sWGA stays
   ~0.2x from 30 to 65 kb (not modelled). gDNA ends are slightly high (1.3-1.5x), so no zone.
6. `--depth_params FILE.json` overrides any of ar, ma, mu, scale, nu, offset, cap, end_width, end_factor,
   depth_mean, depth_spread; the defaults are listed in `vcf-host-depth --help` (no helper command, as the user chose).

Outputs: VCF header `##depth_model=<...>` (all values) and `##depth_host=<host,typical_depth>`;
`host_depth_profiles_<seed>.csv` columns sim_id, ind_id, depth_model, typical_depth. scipy is now a dependency.
Check vs 100 real samples per type: typical depth, correlation (lag 1 and 20), spread and share of windows <5x all
close; upper tail too heavy (99th percentile of window/median 26x vs 15x sWGA, 5.2x vs 2.0x gDNA); gDNA 10th
percentile 0.42 vs 0.32 (too few low windows). Tests: `tests/test_sim_mut.sh` runs both types (outputs renamed
`*_model_<type>.*`) plus checks of the cap, end zone and `--depth_spread 0`.

Open: the user may replace constant-depth `vcf-host` with this; revisit.
