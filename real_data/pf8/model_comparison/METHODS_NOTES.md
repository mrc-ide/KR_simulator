# Methods notes: decisions behind krsim's depth model (Nov 2026)

What was decided after the model fits, why, and the key numbers, for the paper's methods and limitations.
How the models were fitted and chosen (held-out likelihood, BIC, AR(8) t / ARMA(2,2) t) is in
[RECORD.md](RECORD.md). Guiding rule throughout: **keep the model simple; anyone who needs every detail uses real
depth files** (`--depth_zip`).

## 1. Transformation (before fitting)

`y = log(window mean / sample median + 0.05)`, then the sample's mean of y is removed
(`compare_processes.py`, `load`).

- **Divide by the median**: depth becomes relative to each sample's own level, so samples sequenced to very
  different depths (10x to 240x) can be fitted and compared on one scale.
- **Log**: a halving and a doubling become equal-sized steps, so ups and downs are symmetric and the kicks are
  closer to the shapes the models assume. On raw depth, deep samples and big spikes would dominate the fit.
- **+0.05**: windows with zero depth would give log(0) = -inf.
- **Mean removed**: the model describes variation around zero; krsim adds the mean back as `offset`.
- krsim undoes it: relative depth = exp(y + offset) - 0.05.

## 2. The 50 kb trim (fitting only)

- The 50 kb at each end of the 14 nuclear chromosomes (`TRIM_KB = 50`) are left out **when fitting**; apicoplast
  and mitochondrion are left out entirely. About 3% of windows; fitting has ~44,000 per sample.
- Why: sWGA subtelomeric windows are near empty (0.001-0.009x in the first 25 kb). A stationary AR/ARMA process
  cannot describe a long run at log(0.05); left in, it would inflate the variance and long-range correlation.
- Why 50 kb and not 30: the drop is not a clean step (0.09x at 25-30 kb, 0.15-0.3x at 30-60 kb), so 50 kb is a
  buffer. A round number, not tuned. The same trim for gDNA keeps the two fits consistent.
- **No genome is missed in simulation**: krsim simulates every window; the sWGA ends are then lowered by the end
  zone (section 4).
- Methods wording: "The 50 kb at each chromosome end were excluded when fitting, to keep the near-empty
  subtelomeric windows out of the fit. The whole genome was simulated, with the sWGA ends handled separately."

## 3. The cap

- Relative depth (window / median) is capped at the highest window of any QC-pass Pf8 sample of that type:
  **37,623x sWGA, 610x gDNA** (`CAP_PCT = 100` in `make_depth_model_defaults.py`; was the 99th percentile,
  14,293x / 116x). Applied last, after median scaling.
- Why: it prevents unrealistic results. Converting back from log can blow up rare values; without a cap one
  simulated sWGA host reached about 10^16x. With the real maximum, nothing beyond what real data has shown can be
  simulated, but real-sized spikes still can (the user's choice, 2026-11-11).
- It is a safety guard, not part of the fitted model, so no model comparison is needed for it. The check
  (`check_depth_cap.py`, `results/depth_cap/cap_check.tsv`) shows it does not change typical hosts: no cap and
  the max cap give identical KS distances to real data.
- Methods wording: "Relative depth was capped at the highest window observed in any QC-pass Pf8 sample of that
  type (37,623x sWGA, 610x gDNA), to prevent unrealistic values."

## 4. The sWGA end zone

Evidence: `check_swga_ends.py`, `results/swga_ends/` (`summary.txt`, `options.tsv`, `whole_genome.tsv`,
`swga_ends.png`); 300 random real sWGA samples, 200 simulated hosts.

- **What**: the first and last 30 kb of each chromosome are multiplied by 0.0035, before median scaling.
- **Why sWGA only**: real sWGA ends are near empty (median sample: 0.0019x at 0-5 kb, 0.0086x at 20-25 kb,
  0.088x at 25-30 kb, ~0.15-0.3x at 30-60 kb, ~1x by 115 kb). gDNA ends are within its normal range
  (~1.3x at 0-20 kb, ~0.77x at 25-100 kb, against a genome range of 0.32-1.7x), so no zone.
- **Why 30 kb**: the 5 kb end bins where the median over samples is below 0.1x the sample median. The width is
  stable to the threshold: 0.01-0.05x give 25 kb, 0.1-0.15x give 30 kb, 0.2x 35 kb, 0.3x 40 kb.
- **Why 0.0035**: the median real depth of those bins (all 17,330 samples; 0.0033 in the 300-sample check).
  At a typical host depth of 63x that is 0.22 expected reads per site, so ~80% of sites have 0 reads
  (95% at 14x, 35% at 300x); real samples have a median 33% of 0-30 kb windows with no reads at all.
- **Flat zones compared** (score = mean |log10(simulated / real)| of the median end profile; 0 = identical):

  | Zone | Factor | 0-30 kb | 30-100 kb | 0-200 kb |
  |---|---:|---:|---:|---:|
  | none | 1 | 2.33 | 0.42 | 0.53 |
  | 25 kb | 0.0029 | 0.46 | 0.44 | 0.26 |
  | **30 kb (krsim)** | 0.0033 | 0.50 | 0.45 | 0.27 |
  | 50 kb | 0.015 | 0.76 | 0.60 | 0.36 |
  | 100 kb | 0.27 | 1.83 | 0.19 | 0.38 |
  | end profile (one factor per 5 kb bin to 150 kb) | | 0.20 | 0.08 | 0.09 |

  - 25 and 30 kb are near-tied; wider flat zones are worse, because one factor must sit between ~0.003x
    (0-30 kb) and ~0.2x (30-50 kb), so neither part is right.
  - The end profile fits best but was **rejected to keep the model simple** (the user's decision); real depth
    files cover that detail.
- **Whole-genome effect** (median over samples / hosts; KS distance to real, 0 = same distribution):

  | | Real | AR(8) alone | AR(8) + 30 kb zone |
  |---|---:|---:|---:|
  | Share of windows < 0.01x | 3.5% | 1.4% | 4.2% |
  | KS for that share | | 0.80 | 0.57 |
  | KS for the 5th percentile window | | 0.89 | 0.55 |

  - Without the zone, depth within 30 kb of the ends is ~200x too high; with it, about the real level.
- Methods wording: "In sWGA, the 30 kb at each chromosome end, where real median depth is below 0.1x of the
  genome median, were multiplied by 0.0035, the median real depth there. A flat 30 kb zone was chosen over wider
  zones (up to 100 kb) because it gave the closest match to real subtelomeric depth; depth between 30 and 100 kb
  remains somewhat overestimated."

## 5. Comparing simulated with real: summary statistics and KS

1. Each genome (real sample or simulated host) is divided by its median and summarised by a few statistics, e.g.
   the share of windows below 0.01x, 0.1x and 0.5x, and the 1st, 5th and 10th percentile window.
2. For each statistic, the real values and the simulated values are compared with the two-sample
   Kolmogorov-Smirnov distance: the largest vertical gap between their cumulative curves (0 = same
   distribution, 1 = no overlap).

- **Why KS**: it compares whole distributions (typical value and between-sample spread); it is non-parametric
  (the statistics are bounded and skewed); it puts every statistic on a common 0-1 scale; and it gives an effect
  size, not a p-value, which with hundreds of samples would flag trivial differences.
- **Weaknesses**: less sensitive to the extremes than the middle; odd with many tied values (e.g. many zeros).
- **Do not average the KS values**: averaging mixes statistics a change should affect with ones it should not
  (the end zone's mean KS 0.61 -> 0.52 undersells it). Report the statistics that answer the question.
- The KS left over (~0.5) is mostly not the end zone: simulated hosts are more alike than real samples.
- Methods wording: "To compare simulated hosts with real samples, we summarised each genome by the share of
  windows below 0.01x, 0.1x and 0.5x its median depth, and by its 1st, 5th and 10th percentile window. For each
  statistic, the difference between the real and simulated distributions was measured with the two-sample
  Kolmogorov-Smirnov distance, chosen because it is non-parametric, compares whole distributions rather than
  means, and gives an effect size on a common 0-1 scale rather than a p-value."

## 6. Simulating a host's depth in krsim (order of steps)

1. Draw the host's typical (median) depth from a log-normal of Pf8 sample medians (centre 63.0x sWGA / 76.8x
   gDNA).
2. Simulate log depth along every chromosome from the fitted model (AR(8) t sWGA, ARMA(2,2) t gDNA).
3. Undo the transformation: exp(y + offset) - 0.05.
4. sWGA only: apply the end zone.
5. Scale to median 1 over the genome, then cap.
6. Window depth = typical depth x relative depth; reads per site ~ Poisson(window depth).
7. Draw strain proportions from a Dirichlet.
8. Apply the sample-level error rates (WT->MUT, MUT->WT; default 0).
9. Split each site's reads between strains with a multinomial.

## 7. Limitations

- **Spikes too extreme**: a typical simulated host's highest window is ~12x too high (2,545x vs 216x real sWGA;
  133x vs 11x gDNA), and the top 1% of windows are 1.5-2x too high. The cap cannot fix this (the excess is spread
  over the top windows). Still to revisit.
- **Hosts too alike**: all simulated hosts share one parameter set; real samples vary more between each other.
- **sWGA 30-100 kb overestimated**: real depth is ~0.15-0.75x there; the model gives ~1x.
- **Left out on purpose**: gDNA's lopsided dips, sWGA's changing jumpiness (GARCH), and the shared genome map of
  low/high regions (see [RECORD.md](RECORD.md)).

## Feedback on the methods outline (2026-11-11)

- Data: say nuclear chromosomes only, the 50 kb trim, and the number of QC-pass samples per type.
- Transformation: add the mean removal and the full formula.
- Fitting: state the chosen models at the end, not only the candidates.
- Biological realism: the cap and end zone are simulation steps; describe them in the simulation section and
  justify them here (sections 3-4).
- Validation: "reproduces the spread between hosts" contradicts the finding that hosts are too alike; move it to
  the limitations.
