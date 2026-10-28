#!/usr/bin/env bash
# One Pf8 sample, timed: download the CRAM -> mosdepth -> mean and SD per 500 bp window -> delete the rest.
# Run as a SLURM array by pf8_depth_array.sbatch: array task i does line i of the sample list
# (a table with header, made by pf8_crams.py: sample, type, qc_pass, cram_bytes, url).
#
# Keeps only $PF8/output/<sample>.window_stats.tsv, and $PF8/timings/<sample>.tsv with, for each
# stage, the wall-clock seconds, peak memory (GB) and exit code. A sample that fails gets a line in
# $PF8/failed.tsv (its list row, the stage that failed, job and date); pf8_missing.py lists every sample
# still without window stats (this also catches tasks SLURM killed), ready to run again.

# Folder for the data, and the KR_simulator repo (change these on another computer)
PF8=/data/home/bty208/scratch/IMPERIAL/pf8_mosdepth
REPO=/data/home/bty208/IMPERIAL/MRC_IDE/git/KR_simulator

LIST=$1
THREADS=${SLURM_CPUS_PER_TASK:-1}
REF=$PF8/input/PlasmoDB-54-Pfalciparum3D7-Genome.fasta

read -r SAMPLE TYPE QC BYTES URL < <(tail -n +2 $LIST | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "== $SAMPLE ($TYPE, $((BYTES / 1000000)) MB) on $(hostname), $THREADS threads"
if [ -s $PF8/output/$SAMPLE.window_stats.tsv ]; then echo "  already done"; exit 0; fi

mkdir -p $PF8/cram $PF8/output $PF8/timings
CRAM=$PF8/cram/$SAMPLE.cram
PREFIX=$PF8/output/$SAMPLE
TIMES=$PF8/timings/$SAMPLE.tsv
echo -e "sample\ttype\tcram_bytes\tthreads\tnode\tjob\tstage\tseconds\tmax_rss_gb\texit" > $TIMES

# run a stage under /usr/bin/time and add its line to $TIMES; later stages still run if one fails
stage() {
  local name=$1; shift
  /usr/bin/time -f "%e\t%M\t%x" -o $TIMES.tmp "$@"
  awk -v pre="$SAMPLE\t$TYPE\t$BYTES\t$THREADS\t$(hostname)\t$SLURM_ARRAY_JOB_ID\t$name" -F'\t' -v OFS='\t' \
    'END {print pre, $1, sprintf("%.3f", $2 / 1048576), $3}' $TIMES.tmp >> $TIMES
  tail -1 $TIMES
}

# add this sample to $PF8/failed.tsv, then stop
fail() {
  echo "  failed at $1"
  echo -e "$SAMPLE\t$TYPE\t$QC\t$BYTES\t$URL\t$1\t${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}\t$(date +%F)" >> $PF8/failed.tsv
  exit 1
}

# ENA sometimes closes the connection partway through a large CRAM, and refuses new connections when
# many tasks start at once, so: start at a random time in the first minute, then resume (-C -) until the
# size matches ENA's, waiting longer after each failed try (30 s, 60 s, ...), up to 10 tries;
# without the whole CRAM there is no point going on
sleep $((RANDOM % 60))
stage download bash -c '
  for try in $(seq 10); do
    curl -sS --fail -C - -o $1 $2
    [ "$(stat -c %s $1 2>/dev/null)" = "$3" ] && { echo "  downloaded in $try tries"; exit 0; }
    sleep $((30 * try))
  done
  exit 1' _ $CRAM $URL $BYTES
if [ "$(tail -1 $TIMES | cut -f10)" != 0 ]; then
  rm -f $CRAM $TIMES.tmp; fail download
fi
# the index too (ENA refuses some requests when busy, with a 403 or an empty reply)
stage download_crai bash -c '
  for try in $(seq 10); do
    curl -sS --fail -o $1 $2 && [ -s $1 ] && exit 0
    sleep $((30 * try))
  done
  exit 1' _ $CRAM.crai $URL.crai
if [ "$(tail -1 $TIMES | cut -f10)" != 0 ]; then
  rm -f $CRAM $CRAM.crai $TIMES.tmp; fail download_crai
fi
stage mosdepth mosdepth --threads $((THREADS - 1)) --fasta $REF $PREFIX $CRAM
stage window_stats env PYTHONPATH=$REPO/src python -m krsim.output.pf8_profiles \
  $PREFIX.per-base.bed.gz $PREFIX.window_stats.tsv 500

# keep only the window stats
stage cleanup rm -f $CRAM $CRAM.crai $PREFIX.per-base.bed.gz $PREFIX.per-base.bed.gz.csi \
  $PREFIX.mosdepth.global.dist.txt $PREFIX.mosdepth.summary.txt
rm -f $TIMES.tmp
[ -s $PREFIX.window_stats.tsv ] || fail $(awk -F'\t' 'NR > 1 && $10 != 0 {print $7; exit}' $TIMES)
