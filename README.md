# epigenome-integration

## Overview

This is an analysis similar to my previous research, but on a smaller dataset.

This project implements an integrative epigenomics analysis pipeline combining:

- RRBS methylation analysis
- RNA-seq expression analysis
- Integrate two analysis above

The goal is to identify **tissue-specific Differentially Methylated Regions (DMRs)** on CpG island shores and investigate their association with **Differentially Expressed Genes (DEGs)**.

--- 

### Background

DNA methylation in CpG island shores is known to be associated with tissue-specific gene regulation.  
This project focuses on identifying DMRs located in promoter regions and evaluating their relationship with downstream gene expression.

---

### Dataset
- Source: PRJNA880812  
- Samples:
  - Primary colorectal tumors
  - Liver metastases  
- Replicates: 2 biological replicates per condition (limited sample size)

Reference:  
Rodger EJ et al., *An epigenetic signature of advanced colorectal cancer metastasis.*, iScience (2023)  
https://doi.org/10.1016/j.isci.2023.106986

--- 

### Result

Three genes were identified in which DMR on CpG island shores in the promoter region is associated with downstream gene expression.

| Gene | diff.Methylation | log2FC |
|------|----------------|--------|
| PIEZO2 | 0.147 | +3.63 |
| MME | 0.197 | -4.01 |
| LOC124904613 | 0.177 | +3.40 |

Among the candidates, **MME (membrane metalloendopeptidase)** showed:


- increased methylation in CpG island shore of the promoter region
- decreased gene expression

This suggests that DNA methylation **may be associated with gene repression** of MME.

(Fig. 1 is a scatter plot of methylation degree and expression level of three genes. Fig. 2 is the location of DMR, CpGisland, and shore around the MME gene.)

| Fig. 1 | Fig. 2 |
|--------|--------|
| <img width="600" height="700" alt="Fig. 1" src="https://github.com/user-attachments/assets/d7bd443f-2ee5-4116-8531-96715e5cb363" /> | <img width="880" height="720" alt="MME_integration_plot" src="https://github.com/user-attachments/assets/3d0c555a-a4b6-4d7b-8249-f0e425bb9ef5" /> |


---

### Note
* Small sample size limits statistical power
* Results should be interpreted as exploratory
* The following commands was run on GCE (n2-highmem-16 (16 vCPUs, 128 GB Memory, 500GB SSD, Ubuntu 25.10 Minimal)
* It seems the pipelines used in this workflow use up to 12 CPUs and 72 GB of memory without a nextflow.config
  
---

## Setup 
```bash
sudo apt-get update

# Install Docker
sudo apt-get install -y zip unzip docker.io
sudo chmod 666 /var/run/docker.sock

# Install SDKMAN!
curl -s "https://get.sdkman.io" | bash

# Install Java
sdk install java 25.0.2-tem

# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

## Run
```bash
git clone https://github.com/tonyu0/epigenome-integration.git
cd epigenome-integration

# To make sure nextflow will keep running after the session is disconnected, use tmux 
tmux new -s bio-analysis

sh run-nextflow.sh

# if session is disconnected, reconnect like:
# tmux attach -t bio-analysis
```

## Data transfer 
```bash
# GCP example
tar -czvf rnaseq_output.tar.gz rnaseq_output
tar -czvf methylseq_output.tar.gz methylseq_output
gcloud init
gsutil cp rnaseq_output.tar.gz gs://[bucket name]/
gsutil cp methylseq_output.tar.gz gs://[bucket name]/
```

## Data QC
Check MultiQC reports:

* rnaseq_output/multiqc/star_salmon/multiqc_report.html

* methylseq_output/multiqc/bwameth/multiqc_report.html

## ADownstream Analysis with R
```bash
# if R is not installed:
sudo apt-get install r-base r-base-dev
# these are also required for installing package in R
sudo apt-get install libxml2-dev libcurl4-openssl-dev

# find DEGs (saved to DEGs.csv)
Rscript scripts/rnaseq_to_DEGs.R

# find DMRs (saved to DMRs.bed)
Rscript scripts/methylseq_to_DMRs.R
# match format with UCSC's CpGisland.bed
awk  'BEGIN{FS=OFS="\t"}{$1="chr"$1;print$0}' DMRs.bed > DMRs_chr.bed

Rscript scripts/DEGs_DMRs_integration.R
# output a final result (DMR_regulated_genes.csv)

Rscript scripts/visualize_result.R
# output visualized images

```

