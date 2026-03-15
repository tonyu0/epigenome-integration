# epigenome-integration
Tissue-specific DNA methylation and gene expression integration analysis

Epigenome regulation analysis pipeline

- RRBS methylation analysis
- RNA-seq expression analysis
- integrative epigenomics

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