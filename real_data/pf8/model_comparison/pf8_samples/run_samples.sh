#!/usr/bin/env bash
# For each sample in crams.tsv: download the CRAM, get the depth at every base (mosdepth),
# the mean and SD in 500 bp windows, and the depth plots; then delete the CRAM.
# Then rebuild the comparison plots. Re-running skips samples that are already done.
#   bash run_samples.sh
set -euo pipefail
P=/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/pf8_samples
REPO=/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/TSKIT_MALARIA_SIM/git/KR_simulator
BIN=~/miniforge3/envs/mosdepth/bin
PY=$REPO/.venv/bin/python
REF=$P/reference/PlasmoDB-54-Pfalciparum3D7-Genome.fasta
PLOTS=$REPO/real_data/pf8/depth_model

while IFS=$'\t' read -r S GB URL; do
  D=$P/$S
  echo "== $S ($GB GB) $(date +%H:%M)"
  if [ -s $D/depth_plots/$S.window_stats.tsv ]; then echo "  already done"; continue; fi
  mkdir -p $D/mosdepth $D/depth_plots
  cd $D

  # full CRAM and index; retry/resume until the size matches ENA
  until [ "$(stat -f %z $S.cram 2>/dev/null || echo 0)" = "$(curl -sI https://$URL | awk 'tolower($1)=="content-length:"{print $2+0}')" ]; do
    curl -s -C - -o $S.cram "https://$URL" || sleep 5
  done
  curl -s -o $S.cram.crai "https://$URL.crai"

  # depth at every base (only the per-base file is kept)
  $BIN/mosdepth -t 4 --fasta $REF mosdepth/$S $S.cram
  find mosdepth -type f ! -name "$S.per-base.bed.gz" ! -name "$S.per-base.bed.gz.csi" -delete

  # mean and SD in 500 bp windows, and the plots
  PYTHONPATH=$REPO/src $PY -m krsim.output.pf8_profiles mosdepth/$S.per-base.bed.gz depth_plots/$S.window_stats.tsv 500
  $PY $PLOTS/plot_depth.py depth_plots/$S.window_stats.tsv depth_plots/${S}_depth.png
  $PY $PLOTS/plot_depth.py depth_plots/$S.window_stats.tsv depth_plots/${S}_depth_crt.png Pf3D7_07_v3:350000-450000
  $PY $PLOTS/plot_mean_sd.py depth_plots/$S.window_stats.tsv depth_plots/${S}_mean_sd.png

  # the CRAM is no longer needed (crams.tsv has its URL)
  rm -f "$D/${S:?}.cram" "$D/${S:?}.cram.crai"
done < $P/crams.tsv

# comparison plots: compare/sWGA, compare/gDNA, compare/both
bash $P/compare.sh
echo "ALL DONE $(date +%H:%M)"
