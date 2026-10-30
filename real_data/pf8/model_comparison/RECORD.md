# Record: modelling Pf8 read depth along the genome (Oct 2026)

Goal: a third way of making the krsim host VCF, between constant depth (`krsim vcf-host`) and copying a real
Pf8 sample's depth profile (`krsim vcf-host-depth`): a **simple model** of read depth along the genome, fitted
to real Pf8 data, separately for sWGA and gDNA. Capture the main features, not every detail; anyone who needs
every detail uses real depth profiles (Pf8, or their own CRAMs).

## Where everything is

| What | Where |
|---|---|
| Sample list (21 samples) and CRAM URLs | [`pf8_samples/samples.tsv`](pf8_samples/samples.tsv), [`pf8_samples/crams.tsv`](pf8_samples/crams.tsv) |
| Scripts that made the depth data and comparison plots | [`pf8_samples/`](pf8_samples/) (`run_samples.sh`, `compare.sh`, `compare_depth.py`, `compare_samples.py`) |
| Model comparison scripts | this folder: `compare_processes.py`, `plot_processes.py`, `best_simple_model.py`, `compare_kick_models.py`, `compare_arma_kicks.py` (see [README](README.md)) |
| All results (tables, figures) | [`results/`](results/) |
| All figures with explanations, one PDF | [`results/figures_guide/figures_guide.pdf`](results/figures_guide/figures_guide.pdf) |
| Depth data (per-base and 500 bp windows, ~1.6 GB) | **not in git**; on the first laptop in `RA_POSTDOC_IMPERIAL/pf8_samples/`; rebuild with the steps below |

The scripts in `pf8_samples/` and `results/` have the first laptop's paths written in
(`/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/...`); change them at the top of each script on a new laptop.

## 1. Data

- **Samples**: 10 sWGA and 10 gDNA, picked at random from QC-pass Pf8 samples (seeds 2026 and 2027, no
  filtering on file size), plus FP0008-C (gDNA) from earlier: 10 sWGA + 11 gDNA. Listed in `pf8_samples/samples.tsv`.
- **Per sample**: the Pf8 CRAM (ENA project PRJEB96783) → mosdepth per-base depth → mean and SD of depth in
  500 bp windows (`krsim.output.pf8_profiles`) → plots. CRAMs deleted after use (0.2–4.3 GB each, 43 GB total).
- **What a CRAM's size reflects**: mostly how many reads the sample has (3–94 million), which sets its depth.

### Rebuilding the depth data on a new laptop

```bash
# once: conda env with mosdepth + samtools, and the repo's Python env (numpy, pandas, scipy, statsmodels, matplotlib, numba)
mkdir -p <pf8_samples>/reference
# reference: PlasmoDB-54-Pfalciparum3D7-Genome.fasta from Zenodo record 18681980, then samtools faidx
cp real_data/pf8/model_comparison/pf8_samples/* <pf8_samples>/
# edit P, REPO, BIN, PY at the top of run_samples.sh, then:
bash <pf8_samples>/run_samples.sh          # ~5 min per sample, mostly download; keeps per-base depth + window stats + plots
```

## 2. What real depth looks like

- **Depth level** differs a lot between samples (median 10× to 240×), set by how many reads were sequenced.
  gDNA falls into two groups (~10–20× in older studies, ~80–180×).
- **sWGA is much more variable than gDNA** (variance of log depth ~2.0 vs ~0.9): sWGA DNA is amplified with
  primers (phi29), which over-copies some stretches and under-copies others, differently in every sample.
- **~70% of the variation is shared by all samples of a type** (the same regions are low or high everywhere);
  linked to GC content (opposite signs: GC-rich is low in sWGA, high in gDNA).
- **Nearby windows are similar.** sWGA: stretches of a few kb. gDNA: mostly ~0.7 kb ups and downs plus a weak
  slow trend (15–35 kb); adjacent windows extra alike (probably reads crossing window boundaries).
- **Two "gDNA" samples from study 1195-PF-TRAC2-DONDORP (QC0637-C, RCN12993) behave like sWGA**: 4–7× more
  variable, best fitted by sWGA-like models. Probably not standard gDNA; leave them out of a gDNA model or treat separately.
- Exploration: `results/explore_depth.py`, `results/{sWGA,gDNA}_explore.{png,txt}`.

## 3. How models were compared

Each model is fitted to each sample's log depth (`log(window mean / sample median + 0.05)`, 14 nuclear
chromosomes, 50 kb trimmed from each end, mean removed) by maximum likelihood, with its own parameters.

- **Held-out likelihood (main score)**: fit on odd chromosomes, score even ones, then swap. Tests prediction on
  unseen data; no penalty formula needed.
- **BIC (backup)**: penalty ln(n) ≈ 10.7 per parameter. Mostly agrees with held-out.
- **AIC (reported, not relied on)**: penalty 2 per parameter, too small for ~44,000 windows per sample; it
  drifts to the biggest model.
- Scores are on the −2 × log-likelihood scale and shown relative to each sample's best model (0 = best).
- **Absolute checks**: are the leftover kicks uncorrelated, and do extreme kicks happen as often as the model says?

## 4. Results

### a) Correlation structure (Gaussian kicks) — `compare_processes.py`, `results/process_comparison/`

20 models (AR(1)–AR(20), ARMA(1,1)/(2,1)/(2,2), GP kernels) on all 21 samples. iid and HMM left out (no
correlation / same exponential shape as AR(1)).

- **gDNA: ARMA(2,2)**, best in 9 of 9 typical samples under held-out, BIC and AIC; ARMA(3,3) barely better.
- **sWGA: AR(~8)** (looks back ~4 kb); big gains stop there, AR(8)–AR(20) near-tied; best order varies by sample.
- Smooth GP kernels are worst for both (depth is rough, not smooth); AR(1) much worse than the winners.

### b) The kicks are not Gaussian — `results/residual_check.*`

The winners capture the correlation (leftover kicks uncorrelated), but sWGA has ~5× too many extreme kicks
(kurtosis ~11) and calm/jumpy stretches; gDNA's dips are lopsided (big drops ~14× more common than big spikes).

### c) Better kicks for sWGA — `compare_kick_models.py`, `results/kick_models/`

AR(8) with GARCH t kicks (t tails + jumpiness that changes along the genome) is best in 9 of 10 sWGA samples
and well calibrated; ~12,000 better than Gaussian. Left out of the simple model to keep it simple.

### d) The simplest good model — `best_simple_model.py`, `results/simple_models/` (**the current answer**)

7 correlation structures × Gaussian or Student t kicks (t = one extra parameter ν for heavier tails):

- **t kicks beat Gaussian kicks for every structure in every sample**: ~1,170 better for gDNA, ~9,100 for
  sWGA — more than any choice of structure. Extreme kicks then match (sWGA: 1.4% beyond 3 SD, as predicted).
- **gDNA: ARMA(2,2) + t kicks** (ν ≈ 7), best in 8 of 11; next AR(8) t (39 behind), AR(4) t (48).
- **sWGA: AR(8) + t kicks** (ν ≈ 3.6), best in 8 of 10; next ARMA(2,2) t (97 behind), AR(4) t (183).
- One structure for both: ARMA(2,2) + t (best for gDNA, 97 behind for sWGA).

**Deliberately left out**: gDNA's lopsided dips (skewed t), sWGA's changing jumpiness (GARCH), and the shared
genome map of low/high regions. Real depth profiles cover those.

## 5. Concepts (short)

- **Kick**: the random surprise added at each window: new value = prediction from the past + kick.
  Diagram: `results/notes/kicks_diagram.png`.
- **AR(p)**: predict from the previous p depths. **ARMA(p, q)**: p previous depths + q previous kicks.
- **t kicks**: like Gaussian but with heavier tails (more big jumps); ν sets how heavy (small ν = heavier).
- **Held-out log-likelihood**: how unsurprised a model is by chromosomes it was not fitted to.

## 6. Next steps

1. Build the simple depth model into krsim as a third host-VCF option: per host, a depth level, then log depth
   along the genome from ARMA(2,2)+t (gDNA) or AR(8)+t (sWGA) with fitted parameters; reads per site drawn
   around that level (Poisson / negative binomial) and split between strains as in `vcf-host-depth`.
2. Decide the parameter values to ship (e.g. medians over samples) and how the depth level is chosen.
3. Check simulated genomes against real ones (fraction of windows below 5×, spike sizes).
4. Optionally: let `vcf-host-depth` take the user's own CRAM / window stats.
5. The all-Pf8 SLURM run (`real_data/pf8/all_samples/`) could refit the chosen models on many more samples.

Note: `results/MODEL_LIST.md` is an early plan (M0–M8) from before the comparisons; section 4 supersedes it.
