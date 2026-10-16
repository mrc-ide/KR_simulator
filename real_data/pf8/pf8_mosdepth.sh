#!/usr/bin/env bash
# Pf8 read depth with mosdepth. All data goes outside the git repo.

# Folder for the data (change this on another laptop)
PF8=/Users/jbansal/Desktop/IMPERIAL/pf8

# Step 1: sample information
mkdir -p $PF8/input
cd $PF8/input
curl -L -o Pf8_samples.txt "https://zenodo.org/records/18681980/files/Pf8_samples.txt"

# Step 2: how many gDNA and sWGA samples
head -1 Pf8_samples.txt | tr '\t' '\n' | cat -n                                    # numbered column names
tail -n +2 Pf8_samples.txt | cut -f16 | sort | uniq -c                             # samples per Sample type
tail -n +2 Pf8_samples.txt | awk -F'\t' '$14=="True" {print $16}' | sort | uniq -c  # QC-pass only

# Step 3: find the CRAM for one sample in ENA
curl -s -G "https://www.ebi.ac.uk/ena/portal/api/search" \
  --data-urlencode 'result=read_run' \
  --data-urlencode 'query=study_accession=PRJEB96783 AND sample_alias="PF1066-Cx"' \
  --data-urlencode 'fields=run_accession,sample_alias,submitted_bytes,submitted_ftp' \
  --data-urlencode 'format=tsv'

# Step 4: download the alignment (CRAM) and its index
curl -O "https://ftp.sra.ebi.ac.uk/vol1/run/ERR156/ERR15615973/PF1066-Cx.cram"
curl -O "https://ftp.sra.ebi.ac.uk/vol1/run/ERR156/ERR15615973/PF1066-Cx.cram.crai"

# Step 5: download the reference genome and index it
REF=PlasmoDB-54-Pfalciparum3D7-Genome.fasta
curl -L -O "https://zenodo.org/records/18681980/files/PlasmoDB-54-Pfalciparum3D7-Genome.fasta"
samtools faidx PlasmoDB-54-Pfalciparum3D7-Genome.fasta

# Step 6 (optional): explore the alignment
samtools view -H PF1066-Cx.cram | head        # header
samtools idxstats PF1066-Cx.cram              # reads per chromosome
samtools flagstat PF1066-Cx.cram              # mapping summary
samtools view -T $REF PF1066-Cx.cram Pf3D7_07_v3:403000-404000 | head
# ^ reads around crt

# Step 7: depth in windows
mosdepth --fasta $REF -n --by 500 ../output/PF1066-Cx PF1066-Cx.cram     # coverage
