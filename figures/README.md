# Figures

Every figure in the repo has exactly one home here, by topic. The scripts write straight into these folders, so
re-running a script updates the figure in place. There are no copies to keep in step.

- **Website:** the figures the website uses are copied into `docs/images/` by `docs/copy_figures.py` when Quarto
  builds the site. `docs/images/` is not in git. To show a figure on a page, add it to `WEBSITE` in that script.
- **Paper:** mark a figure `main` or `supp` in the Paper column below. Nothing is copied.

Paths in "Made by" are relative to `real_data/pf8/model_comparison/` unless they start with `figures/`,
`benchmarks/` or `analysis/`.

## simulator/

| File | Shows | Made by | Website | Paper |
|---|---|---|---|---|
| `KR_model.png` | The transmission model: hosts over generations, and one transmission through a mosquito (oocysts, recombination) | drawn by hand | background | |
| `KR_tskit.png` | How a simulation is recorded in tskit (individuals, nodes, edges) | drawn by hand | background | |
| `pipeline.png` | The krsim commands and their outputs | `figures/simulator/make_pipeline.py` | home | |
| `scaling_pipeline.png` | Time, peak RAM and output size of the four main commands, 50 to 1,000 hosts | `benchmarks/scaling.py plot` | scaling | |

## pf8_data/

| File | Shows | Made by | Website | Paper |
|---|---|---|---|---|
| `depth_along_genome_PF1066-Cx.png` | Depth along the 14 chromosomes of one sWGA sample | `real_data/pf8/depth_model/plot_depth.py` | realistic-depth | |
| `depth_spread_sWGA_samples.png` | Window mean against SD for 11 sWGA samples, with the Poisson line (log-log): windows are overdispersed | `real_data/pf8/depth_model/plot_compare_samples.py` | realistic-depth | |
| `sWGA_explore.png` | Exploration of sWGA depth: variability, shared regions, correlation | `results/explore_depth.py sWGA` | | |
| `gDNA_explore.png` | The same for gDNA | `results/explore_depth.py gDNA` | | |

## model_screen/ (screen of correlation structures, Gaussian kicks, first 21 samples)

| File | Shows | Made by | Website | Paper |
|---|---|---|---|---|
| `kicks_diagram.png` | What a "kick" is: next value = prediction from the past + kick | `results/notes/kicks_diagram.py` | | |
| `ranking.png` | All 20 structures ranked by held-out score, per type | `plot_processes.py` | | |
| `ranks.png` | Each model's rank in each sample | `plot_processes.py` | | |
| `scores.png` | Held-out and BIC scores per model, sWGA against gDNA | `compare_processes.py` | | |
| `ar_lookback.png` | How far back an AR model needs to look (AR order) | `plot_processes.py` | | |
| `arma_gp.png` | ARMA against Gaussian process kernels | `plot_processes.py` | | |
| `examples.png` | Every sample: real autocorrelation next to AR(1) and the best model | `plot_processes.py` | | |
| `gdna_arma22.png` | gDNA with ARMA(2,2) | `results/gdna/gdna_arma22.py` | | |

## kicks/ (kick distributions)

| File | Shows | Made by | Website | Paper |
|---|---|---|---|---|
| `kick_models.png` | sWGA AR(8) with five kick models (Gaussian, t, GARCH, ...) | `compare_kick_models.py` | | |
| `1_kick_shapes.png` | Gaussian against t kick shapes | `results/kick_gain/plot_kick_gain.py` | | |
| `2_score_gain.png` | Held-out score gained by t kicks, per sample | `results/kick_gain/plot_kick_gain.py` | | |
| `3_tail_check.png` | Are extreme kicks as common as each model predicts? | `results/kick_gain/plot_kick_gain.py` | | |
| `4_all_models.png` | Every structure x kick type, first 21 samples | `results/kick_gain/plot_kick_gain.py` | | |
| `4b_all_models_new20.png` | The same on the 20 new samples (replication) | `results/kick_gain/plot_kick_gain.py` | | |
| `5_kicks_compared.png` | Simulated Gaussian against t kicks, and the depth they give | `results/kick_gain/plot_kick_gain.py` | | |
| `6_how_scored.png` | How the held-out score is built (odd/even chromosomes) | `results/kick_gain/plot_kick_gain.py` | | |
| `7a_real_gDNA_PD0048-C.png` | A real gDNA sample next to genomes simulated from its Gaussian and t fits (needs depth data) | `results/kick_gain/plot_kick_gain.py` | | |
| `7b_real_sWGA_SPT43924.png` | The same for a real sWGA sample (chromosome 7) | `results/kick_gain/plot_kick_gain.py` | | |

## model_choice/ (final model choice)

| File | Shows | Made by | Website | Paper |
|---|---|---|---|---|
| `simple_models.png` | 7 structures x Gaussian/t kicks: held-out scores and tail calibration, first 21 samples | `best_simple_model.py <pf8_dir> results/simple_models` | | |
| `simple_models_new20.png` | The same on 20 new samples (replication) | `best_simple_model.py <pf8_dir> results/simple_models_new20` | | |

## realism/ (cap and sWGA ends)

| File | Shows | Made by | Website | Paper |
|---|---|---|---|---|
| `cap_compare.{png,pdf}` | Each host's highest window: real samples against the model with no cap and with the max cap | `plot_cap_compare.py` | | |
| `ends_profile.{png,pdf}` | Depth by distance from the chromosome end (0-200 kb), with the 30 kb end zone and 50 kb fitting trim | `plot_ends_vs_genome.py` | | supp |
| `ends_vs_genome.{png,pdf}` | Boxplots of depth at 0-25, 25-50, 50-100 and 100-200 kb from the end against the rest of the genome | `plot_ends_vs_genome.py` | | |
| `swga_ends.png` | The end-zone check: real end profile, flat zones of 0-100 kb against an end profile | `check_swga_ends.py` | | |

## coalescence/ (coalescence-time sweep)

`paper/` uses the write-up's parameters (baseline Nh 100, K 5, chi 0.5, C 1); `real_scale/` the real-scale page's
(baseline Nh 100, K 10, chi 0.2, C 3). Both from `analysis/coalescence_sweep/coal_sweep.py plot` (`COAL_SET=real_scale`
for the second). Rows: pairs of strains in the same host and in different hosts; columns: Nh, K, C, chi.

| File | Shows | Website | Paper |
|---|---|---|---|
| `coal_sweep_density_20gen.{png,pdf}` | Share of the genome coalescing per generation, generations 0-20 (the early rise for different hosts) | | |
| `coal_sweep_density_100gen.{png,pdf}` | The same, generations 0-100, one point per generation | | |
| `coal_sweep_density.{png,pdf}` | The same, generations 0-1,000, 20-generation periods | | |
| `coal_sweep.{png,pdf}` | Cumulative: fraction of the genome coalesced by each generation (different hosts) | | |
