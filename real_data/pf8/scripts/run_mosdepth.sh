#!/usr/bin/env bash
# Coverage (mosdepth) for Pf8 samples.
#
# For each sample ID given (the `Sample` column of metadata/Pf8_samples.txt):
#   1. finds its Pf8 CRAM in ENA project PRJEB96783 (matched on sample name)
#   2. downloads the CRAM and its index into crams/ (skipped if already there)
#   3. runs mosdepth, writing results to mosdepth/<sample>/
#
# Usage (from anywhere, with `conda activate mosdepth`):
#   bash scripts/run_mosdepth.sh FP0008-C
#   bash scripts/run_mosdepth.sh FP0008-C FP0010-CW
#
# Set DELETE_CRAM=1 to delete each CRAM after mosdepth (they are 0.5-4 GB each):
#   DELETE_CRAM=1 bash scripts/run_mosdepth.sh FP0008-C

set -euo pipefail

# folders, relative to this script (so it works from any directory)
PF8="$(cd "$(dirname "$0")/.." && pwd)"
REF="$PF8/reference/PlasmoDB-54-Pfalciparum3D7-Genome.fasta"
CRAMS="$PF8/crams"
OUT="$PF8/mosdepth"

# mosdepth settings
WINDOW=1000          # also report mean depth in windows of this many bp
THRESHOLDS=1,5,10    # bases per window with at least this many reads (Pf8 "callable" uses 5)

if [ $# -eq 0 ]; then
    echo "Usage: bash scripts/run_mosdepth.sh SAMPLE [SAMPLE ...]   e.g. FP0008-C" >&2
    exit 1
fi

for tool in curl samtools mosdepth; do
    command -v "$tool" >/dev/null || { echo "Missing $tool: run 'conda activate mosdepth' first" >&2; exit 1; }
done

# reference genome (downloaded once)
if [ ! -f "$REF" ]; then
    echo "Downloading reference genome..."
    mkdir -p "$(dirname "$REF")"
    curl -L -o "$REF" "https://zenodo.org/records/18681980/files/PlasmoDB-54-Pfalciparum3D7-Genome.fasta"
fi
[ -f "$REF.fai" ] || samtools faidx "$REF"

mkdir -p "$CRAMS" "$OUT"

for SAMPLE in "$@"; do
    echo "== $SAMPLE"
    CRAM="$CRAMS/$SAMPLE.cram"

    if [ ! -f "$CRAM" ] || [ ! -f "$CRAM.crai" ]; then
        # ask ENA for this sample's file paths in the Pf8 project
        PATHS=$(curl -s -G "https://www.ebi.ac.uk/ena/portal/api/search" \
            --data-urlencode "result=read_run" \
            --data-urlencode "query=study_accession=PRJEB96783 AND sample_alias=\"$SAMPLE\"" \
            --data-urlencode "fields=submitted_ftp" \
            --data-urlencode "format=tsv" | tail -n +2 | cut -f2)
        CRAM_URL=$(echo "$PATHS" | tr ';' '\n' | grep '\.cram$' || true)
        if [ -z "$CRAM_URL" ]; then
            echo "  no Pf8 CRAM found in ENA for $SAMPLE (check the name, or it may not have one); skipping" >&2
            continue
        fi
        echo "  downloading https://$CRAM_URL"
        curl -L --fail -o "$CRAM" "https://$CRAM_URL"
        curl -L --fail -o "$CRAM.crai" "https://$CRAM_URL.crai"
    else
        echo "  CRAM already downloaded"
    fi

    mkdir -p "$OUT/$SAMPLE"
    echo "  running mosdepth"
    mosdepth --fasta "$REF" -n -x --by "$WINDOW" --thresholds "$THRESHOLDS" \
        "$OUT/$SAMPLE/$SAMPLE" "$CRAM"
    echo "  done: $OUT/$SAMPLE/"
    awk '$1 == "total"' "$OUT/$SAMPLE/$SAMPLE.mosdepth.summary.txt" \
        | awk '{printf "  mean depth across the genome: %.1fx\n", $4}'

    if [ "${DELETE_CRAM:-0}" = "1" ]; then
        rm -f "$CRAM" "$CRAM.crai"
        echo "  deleted the CRAM (DELETE_CRAM=1)"
    fi
done
