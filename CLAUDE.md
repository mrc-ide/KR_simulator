# Notes for Claude

## Simple depth model for the host VCF (status, 2026-10-30)

Full record: [real_data/pf8/model_comparison/RECORD.md](real_data/pf8/model_comparison/RECORD.md).

**Chosen models** (log window depth along the genome, fitted to Pf8 samples by `best_simple_model.py`):

| Type | Model | ν (median) | Best in |
|---|---|---|---|
| gDNA | ARMA(2,2) + Student t kicks | ≈ 7.0 | 8 of 11 samples |
| sWGA | AR(8) + Student t kicks | ≈ 3.6 | 8 of 10 samples |

- The t kick is location-scale (μ, scale, ν); ν is fitted as `2.05 + exp(θ)`
  (`compare_arma_kicks.py`, `unpack` / `t_logpdf`).
- t kicks beat Gaussian kicks for every structure in every sample.
- Checks done: held-out likelihood + BIC (model choice); kicks uncorrelated; tail calibration
  (sWGA t is well calibrated; gDNA t still has too many big drops, too few big spikes).
- Known misses, left out on purpose: sWGA calm/jumpy stretches (GARCH), gDNA lopsided dips (skewed t),
  the shared genome-wide map of low/high regions.

**Agreed goal (from the user):** the model does not need to be precise or statistically tight. It is a
realistic starting point, better than one depth value across the whole genome. Users who need detail
supply their own depth data (`krsim vcf-host-depth`). Don't push for further model refinement.

**Not yet in krsim**: nothing in `src/` uses these models.

**Next steps (waiting for the user's go-ahead before changing code):**
1. Optional quick visual check: simulate a few genomes from each fitted model and plot them next to real samples.
2. Choose the parameters to ship (e.g. median over samples per type; maybe leave out the two sWGA-like gDNA
   samples from 1195-PF-TRAC2-DONDORP, QC0637-C and RCN12993).
3. Decide how the depth level is set (user value or drawn from the real samples' range).
4. Build it into krsim as a third host-VCF option: per host, log depth from the model, reads per site drawn
   around that depth, split between strains as in `vcf-host-depth`.
