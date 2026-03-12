#!/bin/bash

# Execute this shell after installing Docker and Nextflow

# ====================
# Fetch reference data
# ====================

# Download Human genome
wget -L ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
# Download annotation
wget -L ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip *.gz

# Limit the target genome to only chromosome 20
# (also: samtools faidx GRCh38.fasta && samtools faidx GRCh38.fasta outChr20 > outChr20.fa)
perl scripts/extract-selected-chr-from-fa.pl Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa 20
# Limit the target annotation to only chromosome 20
sh scripts/extract-selected-chr-from-gtf.sh Homo_sapiens.GRCh38.112.gtf 20
# Now outChr20.fa and outChr20.gtf exist in the current directory

# ====================
# Run Nextflow (fetchngs, rnaseq, methylseq)
# ====================

nextflow pull nf-core/fetchngs
nextflow pull nf-core/rnaseq
nextflow pull nf-core/methylseq

# Fetch SRA data in ids.csv
nextflow run nf-core/fetchngs -profile docker --input ids.csv --outdir ./fetchngs_output

# output csv did not have the field 'strandedness', so add it
awk 'BEGIN{FS=OFS=","} NR==1{print $0, "\"strandedness\""} NR>1{print $0, "\"auto\""}' ./fetchngs_output/samplesheet/samplesheet.csv > ./fetchngs_output/samplesheet/samplesheet-with-strandedness.csv

# Run RNA-seq pipeline
nextflow run nf-core/rnaseq -profile docker --input ./fetchngs_output/samplesheet/samplesheet-with-strandedness.csv --outdir ./rnaseq_output --fasta outChr20.fa --gtf outChr20.gtf --aligner star_salmon

# Run Bisulfite-seq(RRBS) pipeline
nextflow run nf-core/methylseq -profile docker --input ./fetchngs_output/samplesheet/samplesheet-with-strandedness.csv --outdir ./methylseq_output --fasta outChr20.fa --rrbs --aligner bismark
