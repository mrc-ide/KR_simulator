# Depth models to test (sWGA first, then gDNA)

Data: 10 sWGA samples, 500 bp windows, y = log(window mean depth + 1).
Findings that drive the list (`sWGA_explore.txt`, `sWGA_explore.png`):

- Depth level differs a lot between samples: median 23x to 240x.
- 70% of the variance in log depth is **shared** by all samples (same windows low/high); 30% is sample-specific.
- The shared effect is lower in GC-rich windows (r = -0.45). A cluster at about 32% GC sits around 30x below average.
- The sample-specific part is correlated along the genome, and an **AR(1) fits its autocorrelation closely**
  (0.86 at 500 bp, ~0 by 10 kb). A small extra tail (~0.04 at 10-50 kb) is the only miss.
- Its tails are heavier than Gaussian in some samples (excess kurtosis up to 3).
- Within a window, variance / mean is about 2 (Poisson is 1).

Every model has the same three layers; the models differ in which parts they include.

```
log mean depth (sample s, window w) = sample level a_s + shared effect b_w + sample-specific noise z_sw
depth at a site in window w         ~ Poisson or negative binomial around exp(log mean)
```

| Model | Sample level a_s | Shared effect b_w | Sample-specific z | Site depth | Role |
|---|---|---|---|---|---|
| M0 | fixed (`--depth`) | none | none | Poisson | current `vcf-host`, baseline |
| M1 | one real sample | real profile | real profile | NB from window mean/SD | current `vcf-host-depth`, realistic reference |
| M2 | drawn from the 10 samples | none | independent Gaussian | NB | spread only, no position |
| M3 | drawn | mean profile across samples | independent Gaussian | NB | adds the shared effect |
| **M4** | drawn | mean profile | **AR(1)** | NB | main candidate |
| M5 | drawn | mean profile | AR(1), Student-t steps | NB | heavier tails |
| M6 | drawn | mean profile | sum of two AR(1) (or GP, Matérn) | NB | two length scales |
| M7 | drawn | **predicted from GC** (+ chromosome-end term) | AR(1) | NB | needs no stored profile |
| M8 | drawn | mean profile | 2-state HMM (amplified / dropped) + AR(1) | NB | dropouts and spikes |

Options to test within a model:
- Site depth: Poisson vs negative binomial (dispersion from the within-window variance / mean of 2).
- Sites in the same window: independent draws (as `vcf-host-depth` now) vs one continuous process at site resolution.

## How to compare

Leave one sample out: fit on 9 samples, test on the 10th, repeat for all 10.

1. **Likelihood** of the held-out sample's window depths (for M2-M8), plus AIC/BIC for model size.
2. **Summary statistics**, simulated vs real held-out sample:
   - distribution of window depth (Wasserstein distance on log depth)
   - fraction of windows below 5x and at 0
   - autocorrelation curve out to 100 kb
   - mean-vs-SD relationship
3. **Downstream** (what matters for krsim): run the host VCF with each model on the same simulation and compare
   missing genotypes, heterozygous sites called as homozygous, and allele-depth spread against M1.
