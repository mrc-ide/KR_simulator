#!/usr/bin/env bash
# Steps 3, 4, 7 and 8 of pf8_mosdepth.sh for every sample in a list, one after another.
# Needs steps 1 and 5 done first (Pf8_samples.txt and the reference in $PF8/input).
# The CRAM is deleted after each sample; the per-base depths and window stats are kept.
#
# Usage (with `conda activate mosdepth`):
#   bash run_samples.sh random10_sWGA.tsv      # a list in $PF8/input, first column = sample ID

# Folder for the data, and the KR_simulator repo (change these on another laptop)
PF8=/Users/jbansal/Desktop/IMPERIAL/pf8
REPO=/Users/jbansal/Desktop/IMPERIAL/KR_tskit/git/KR_simulator

cd $PF8/input
REF=PlasmoDB-54-Pfalciparum3D7-Genome.fasta

for SAMPLE in $(tail -n +2 $1 | cut -f1); do
  echo "== $SAMPLE"
  if [ -s ../output/$SAMPLE.window_stats.tsv ]; then echo "  already done"; continue; fi

  # Step 3: find the CRAM for this sample in ENA
  URL=$(curl -s -G "https://www.ebi.ac.uk/ena/portal/api/search" \
    --data-urlencode 'result=read_run' \
    --data-urlencode "query=study_accession=PRJEB96783 AND sample_alias=\"$SAMPLE\"" \
    --data-urlencode 'fields=submitted_ftp' \
    --data-urlencode 'format=tsv' | tail -n +2 | cut -f2 | tr ';' '\n' | grep '\.cram$')

  # Step 4: download the CRAM and its index
  curl -s -O "https://$URL"
  curl -s -O "https://$URL.crai"

  # Step 7: depth at every base
  mosdepth --fasta $REF ../output/$SAMPLE $SAMPLE.cram

  # Step 8: mean and SD of depth in 500 bp windows
  PYTHONPATH=$REPO/src python -m krsim.output.pf8_profiles \
    ../output/$SAMPLE.per-base.bed.gz ../output/$SAMPLE.window_stats.tsv 500

  # the CRAM is no longer needed
  rm $SAMPLE.cram $SAMPLE.cram.crai
done
