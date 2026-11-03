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

**Not yet in krsim**: nothing in `src/` uses these models.

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

**Next steps (waiting for the user's go-ahead before changing code):**
1. Choose the parameters to ship (e.g. median over samples per type; maybe leave out the sWGA-like gDNA
   samples QC0637-C, RCN12993, SPT00902).
2. Decide how the depth level is set (user value or drawn from the real samples' range).
3. Build it into krsim as a third host-VCF option: per host, log depth from the model, reads per site drawn
   around that depth, split between strains as in `vcf-host-depth`.
4. Maybe update RECORD.md with the new-20 results and the corrected wording above.
