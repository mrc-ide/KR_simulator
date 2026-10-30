#!/usr/bin/env bash
# Comparison plots from every finished sample: compare/sWGA, compare/gDNA, and compare/both.
#   bash compare.sh
set -euo pipefail
P=/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/pf8_samples
PY=/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/TSKIT_MALARIA_SIM/git/KR_simulator/.venv/bin/python
cd $P
rm -rf compare
for T in sWGA gDNA both; do
  F=()
  for S in $(awk -F'\t' -v t=$T 'NR>1 && (t=="both" || $2==t) {print $1}' samples.tsv); do
    f=$S/depth_plots/$S.window_stats.tsv
    [ -s $f ] && F+=($f)
  done
  mkdir -p compare/$T
  $PY compare_samples.py compare/$T/mean_sd_all.png compare/$T/summary.tsv "${F[@]}" > /dev/null
  $PY compare_depth.py compare/$T "${F[@]}" > /dev/null
  # the per-type depth plots already live in compare/sWGA and compare/gDNA
  [ $T = both ] && rm -f compare/both/depth_sWGA.png compare/both/depth_gDNA.png
  echo "compare/$T: ${#F[@]} samples"
done
