# epigenome-integration

## Overview

Tissue-specific Differentially Methylated Region (DMR) on CpG island shores and Tissue-specific Differentially Expressed Genes (DEGs) integration analysis

Epigenome regulation analysis pipeline

- RRBS methylation analysis
- RNA-seq expression analysis
- Integrate two analysis above

### Result
One gene was identified in which DMR on CpG island shores in the promoter region appears to influence downstream gene expression, i.e. MME (membrane metalloendopeptidase).

### Data Source
RNA-seq & RRBS Data: PRJNA880812 (Liver, Primary samples)
Reference Paper: Rodger EJ et al., "An epigenetic signature of advanced colorectal cancer metastasis.", iScience, 2023 Jun 16;26(6):106986. doi:
10.1016/j.isci.2023.106986
https://doi.org/10.1016/j.isci.2023.106986

### Note
The following commands was run on GCE (n2-highmem-16 (16 vCPUs, 128 GB Memory, 500GB SSD, Ubuntu 25.10 Minimal)
It seems the pipelines used in this workflow use up to 12 CPUs and 72 GB of memory without a nextflow.config

## Setup 
```
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
```
git clone https://github.com/tonyu0/epigenome-integration.git
cd epigenome-integration

# To make sure nextflow will keep running after the session is disconnected, use tmux 
tmux new -s bio-analysis

sh run-nextflow.sh

# if session is disconnected, reconnect like:
# tmux attach -t bio-analysis
```

## Data transfer 
```
# GCP example
tar -czvf rnaseq_output.tar.gz rnaseq_output
tar -czvf methylseq_output.tar.gz methylseq_output
gcloud init
gsutil cp rnaseq_output.tar.gz gs://[bucket name]/
gsutil cp methylseq_output.tar.gz gs://[bucket name]/
```

## Data check
* rnaseq_output/multiqc/star_salmon/multiqc_report.html

* methylseq_output/multiqc/bwameth/multiqc_report.html

## Analysis using R
```R
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

## Interpretation of results

In above analysis, three genes are extracted.
(See "DMRs_regulated_genes.csv")

ENSG00000154864(diff.Methy: 0.14668812964347, log2FoldChange: 3.63034872976175)
PIEZO2 (piezo type mechanosensitive ion channel component 2)

ENSG00000196549(diff.Methy: 0.196966182641586, log2FoldChange: -4.00606354053707)
MME (membrane metalloendopeptidase)

ENSG00000273768(diff.Methy: 0.176937585802774, log2FoldChange: 3.39530230478619)
LOC124904613 (U1 spliceosomal RNA)

Among these, MME shows a negative correlation between methylation level and expression level, suggesting that DMR may be involved in its expression regulation. 
