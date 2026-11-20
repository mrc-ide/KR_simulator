# Figures for the krsim paper

- **`all/`**: every figure in the repo, one copy each, named by topic so they sort together. A catalogue to pick
  from; not all are paper quality.
- **`main/`**: figures chosen for the main text (none yet).
- **`supplementary/`**: figures chosen for the supplement.

These are **copies**. Each figure is made by the script listed below, into its own folder; after re-running a script,
copy the new figure here again.

Prefixes: `01` simulator, `02` Pf8 depth data, `03` screen of correlation structures, `04` kick distributions,
`05` final model choice, `06` realism adjustments (cap, sWGA ends).

## Chosen

| File | Shows | Made by |
|---|---|---|
| `supplementary/S_ends_profile.{png,pdf}` | Depth by distance from the chromosome end (0-200 kb), 300 sWGA and 300 gDNA samples; marks the 30 kb sWGA end zone and the 50 kb fitting trim | `real_data/pf8/model_comparison/plot_ends_vs_genome.py` |

## All figures

Paths in "Made by" are relative to `real_data/pf8/model_comparison/` unless they start with `docs/` or `benchmarks/`.

### 01 Simulator

| File | Shows | Made by |
|---|---|---|
| `01_simulator_transmission_model.png` | The transmission model: hosts over generations, and one transmission through a mosquito (oocysts, recombination) | drawn by hand (`docs/images/KR_model.png`) |
| `01_simulator_tskit_recording.png` | How a simulation is recorded in tskit (individuals, nodes, edges) | drawn by hand (`docs/images/KR_tskit.png`) |
| `01_simulator_pipeline.png` | The krsim commands and their outputs | `docs/images/make_pipeline.py` |
| `01_simulator_scaling.png` | Time, peak RAM and output size of the four main commands, 50 to 1,000 hosts | `benchmarks/scaling.py` |

### 02 Pf8 depth data

| File | Shows | Made by |
|---|---|---|
| `02_data_depth_along_genome_PF1066-Cx.png` | Depth along the 14 chromosomes of one sWGA sample | `real_data/pf8/depth_model/plot_depth.py` |
| `02_data_window_mean_vs_sd_sWGA.png` | Window mean against SD for 11 sWGA samples, with the Poisson line (log-log): windows are overdispersed (why depth zips use a negative binomial) | `real_data/pf8/depth_model/plot_compare_samples.py` |
| `02_data_explore_sWGA.png` | Exploration of sWGA depth: variability, shared regions, correlation | `results/explore_depth.py` |
| `02_data_explore_gDNA.png` | The same for gDNA | `results/explore_depth.py` |

### 03 Screen of correlation structures (Gaussian kicks, first 21 samples)

| File | Shows | Made by |
|---|---|---|
| `03_screen_kicks_diagram.png` | What a "kick" is: next value = prediction from the past + kick | `results/notes/kicks_diagram.py` |
| `03_screen_ranking.png` | All 20 structures ranked by held-out score, per type (ARMA(2,2) for gDNA, AR(~8) for sWGA) | `plot_processes.py` |
| `03_screen_ranks.png` | Each model's rank in each sample | `plot_processes.py` |
| `03_screen_scores.png` | Held-out and BIC scores per sample | `compare_processes.py` / `plot_processes.py` |
| `03_screen_ar_lookback.png` | How far back an AR model needs to look (AR order) | `plot_processes.py` |
| `03_screen_arma_gp.png` | ARMA against Gaussian process kernels | `plot_processes.py` |
| `03_screen_examples.png` | Example samples with fitted models | `plot_processes.py` |
| `03_screen_gdna_arma22.png` | gDNA with ARMA(2,2) | `results/gdna/gdna_arma22.py` |

### 04 Kick distributions

| File | Shows | Made by |
|---|---|---|
| `04_kicks_swga_kick_models.png` | sWGA AR(8) with five kick models (Gaussian, t, GARCH, ...) | `compare_kick_models.py` |
| `04_kicks_1_shapes.png` | Gaussian against t kick shapes | `results/kick_gain/plot_kick_gain.py` |
| `04_kicks_2_score_gain.png` | Held-out score gained by t kicks, per sample | `results/kick_gain/plot_kick_gain.py` |
| `04_kicks_3_tail_check.png` | Are extreme kicks as common as each model predicts? | `results/kick_gain/plot_kick_gain.py` |
| `04_kicks_4_all_models.png` | Every structure x kick type, first 21 samples | `results/kick_gain/plot_kick_gain.py` |
| `04_kicks_4b_all_models_new20.png` | The same on the 20 new samples (replication) | `results/kick_gain/plot_kick_gain.py` |
| `04_kicks_5_compared.png` | Simulated Gaussian against t kicks, and the depth they give | `results/kick_gain/plot_kick_gain.py` |
| `04_kicks_6_how_scored.png` | How the held-out score is built (odd/even chromosomes) | `results/kick_gain/plot_kick_gain.py` |
| `04_kicks_7a_real_vs_sim_gDNA.png` | A real gDNA sample next to genomes simulated from its Gaussian and t fits | `results/kick_gain/plot_kick_gain.py` |
| `04_kicks_7b_real_vs_sim_sWGA.png` | The same for a real sWGA sample (chromosome 7) | `results/kick_gain/plot_kick_gain.py` |

### 05 Final model choice

| File | Shows | Made by |
|---|---|---|
| `05_selection_simple_models.png` | 7 structures x Gaussian/t kicks: held-out scores and tail calibration, first 21 samples (AR(8) t sWGA, ARMA(2,2) t gDNA) | `best_simple_model.py` |
| `05_selection_simple_models_new20.png` | The same on 20 new samples (replication) | `best_simple_model.py` |

### 06 Realism adjustments

| File | Shows | Made by |
|---|---|---|
| `06_realism_cap_compare.{png,pdf}` | Each host's highest window: real samples against the model with no cap and with the max cap | `plot_cap_compare.py` |
| `06_realism_ends_profile.{png,pdf}` | Depth by distance from the chromosome end, with the 30 kb end zone and 50 kb trim (also in `supplementary/`) | `plot_ends_vs_genome.py` |
| `06_realism_ends_vs_genome.{png,pdf}` | Boxplots of depth at 0-25, 25-50, 50-100 and 100-200 kb from the end against the rest of the genome | `plot_ends_vs_genome.py` |
| `06_realism_swga_ends_check.png` | The end-zone check: real end profile, flat zones of 0-100 kb against an end profile | `check_swga_ends.py` |
