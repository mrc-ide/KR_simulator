# Which process describes Pf8 read depth?

**Start with [RECORD.md](RECORD.md)**: the whole analysis, results and next steps in one place.

[`compare_processes.py`](compare_processes.py) fits candidate random processes to the read depth of real Pf8 samples along the genome,
and compares them, separately for sWGA and gDNA. The models and likelihood are those of
[`../depth_model/depth_process_models.ipynb`](../depth_model/depth_process_models.ipynb), run for every sample in one go.

## Run

```bash
python compare_processes.py <pf8_dir> <output_dir>          # --types sWGA gDNA (default both), --workers N
```

`<pf8_dir>` holds `samples.tsv` (columns `Sample`, `Sample type`, …) and `<sample>/depth_plots/<sample>.window_stats.tsv`
(mean and SD of depth per 500 bp window, from `krsim.output.pf8_profiles`). Needs numpy, pandas, scipy, statsmodels and matplotlib.

## What it fits

The series for each sample is `y = log(window mean depth / sample median + 0.05)` on the 14 nuclear chromosomes,
with 50 kb trimmed from each chromosome end and the sample's mean removed. Each chromosome is a separate series.
Every model is fitted to every sample by maximum likelihood, with its own parameters.

| Family | Models | Used? | Why |
|---|---|---|---|
| No correlation | iid Gaussian | No | Depth in neighbouring windows is strongly correlated, so it cannot fit |
| Autoregressive | AR(1), AR(2) … AR(10), AR(15), AR(20) | Yes | AR(1) is the simplest correlated process; higher orders test whether depth depends on windows further back |
| ARMA | ARMA(1,1), (2,1), (2,2) | Yes | A wandering level plus noise (1,1); a fast plus a slow process (2,1) |
| Gaussian process | exponential, Matérn 3/2, Matérn 5/2, squared exponential, two exponentials | Yes | Test the shape of the correlation: sharp (exponential = AR(1)) or smooth at short distances, one or two length-scales |
| Hidden Markov model | 2 or 3 states | No | Its correlation has the same exponential shape as AR(1); it would make depth jump between a few fixed levels |

## How models are compared

All three scores are on the −2 × log-likelihood scale. For each sample, the score of every model is given relative to the best model for that sample (0 = best, larger = worse).

| Score | What it is | Notes |
|---|---|---|
| Likelihood (LL) | How probable the observed depth is under the fitted model | Never used alone: a bigger model always fits at least as well |
| AIC | 2k − 2 LL | A small penalty per parameter (k); with ~44,000 windows per sample it tends to keep adding parameters |
| BIC | k ln(n) − 2 LL | About 5× AIC's penalty per parameter here; prefers smaller models |
| Held-out | −2 × LL on chromosomes the model was not fitted to (fit on odd chromosomes, score even ones, then swap) | No penalty needed; extra parameters only help if they generalise to new data |

Held-out likelihood is the most direct test, and BIC is the safer of the two penalised scores with this much data.
Where they agree, the result is solid.

## Outputs

| File | Contents |
|---|---|
| `fits.tsv` | Every sample × model: LL, k, n, AIC, BIC, held-out LL, the three scores (`dAIC`, `dBIC`, `dHeldout`) and ranks |
| `summary.tsv` | Per sample type × model: how many samples it wins, its median rank and median score under each criterion |
| `scores.png` | Each model's score in every sample (dots) and the median (bar), sWGA and gDNA side by side |

[`plot_processes.py`](plot_processes.py) `<pf8_dir> <output_dir>` adds easier plots to the same folder:

| File | Contents |
|---|---|
| `ranking.png` | Models from best to worst for each sample type (median held-out score), with how many samples each one wins |
| `ranks.png` | The rank of every model in every sample (1 = best), sWGA above gDNA |
| `examples.png` | Every sample: real autocorrelation over the first 5 kb next to AR(1) and the best model of its type, ordered by how far AR(1) is behind |
| `ar_lookback.png` | AR(p): median held-out score against how far back the model looks (p × 0.5 kb), sWGA and gDNA |
| `arma_gp.png` | The ARMA and Gaussian-process models: median held-out score, sWGA and gDNA |

## Better kicks, and the simplest good model

`compare_processes.py` assumes Gaussian kicks (the random surprise added at each window). Checking what the
winners leave over showed the correlation is captured but the kicks are not Gaussian: far too many extreme
jumps (heavy tails), lopsided dips in gDNA, and calm and jumpy stretches in sWGA. These scripts change the kicks.
They use a conditional likelihood (the first windows of each chromosome are only a run-in), so every model is
scored on the same windows; held-out and BIC as above.

| Script | What it compares |
|---|---|
| [`best_simple_model.py`](best_simple_model.py) | AR(1), AR(2), AR(4), AR(8), ARMA(1,1), ARMA(2,1), ARMA(2,2), each with Gaussian or Student t kicks (t = one extra parameter, ν, for heavier tails), for every sample, plus a calibration check. Writes `simple_fits.tsv`, `simple_summary.tsv`, `simple_calibration.tsv`, `simple_models.png` |
| [`compare_kick_models.py`](compare_kick_models.py) | sWGA: AR(8) with Gaussian, t, GARCH (jumpiness that changes along the genome), GARCH t, or a 2-state σ (calm / jumpy). Needs numba |
| [`compare_arma_kicks.py`](compare_arma_kicks.py) | ARMA(p, q) for p, q ≤ 3, then skewed t, GARCH and dropout kicks; also the fitting code that `best_simple_model.py` uses |

```bash
python best_simple_model.py <pf8_dir> <output_dir>
```

Results for the 10 sWGA and 11 gDNA samples:

- **t kicks beat Gaussian kicks for every correlation structure in every sample**: about 1,170 better (held-out) for gDNA
  and about 9,100 for sWGA, more than any choice of correlation structure. The extreme kicks then happen about as
  often as the model says (sWGA: 1.4% beyond 3 SD, as a t with ν ≈ 3.6 predicts).
- **gDNA: ARMA(2,2) with t kicks** (ν ≈ 7), best in 8 of 11 samples. Its dips are still lopsided (more down than up).
- **sWGA: AR(8) with t kicks** (ν ≈ 3.6), best in 8 of 10 samples. GARCH t fits better still, but is left out to keep the model simple.
- Deliberately left out of the simple models: lopsided dips, changing jumpiness and the genome-wide map of low and
  high regions shared by all samples. Real depth profiles (`krsim vcf-host-depth`) cover those.
