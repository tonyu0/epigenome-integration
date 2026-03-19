# ========== Load libraries ==========
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "SummarizedExperiment", "dplyr"), ask=F)

library(DESeq2)
library(SummarizedExperiment)
library(dplyr)

# ---------- Memo for how to deal with output data from star_salmon ---------- 
# rnaseq_output/star_salmon/salmon.merged.gene.SummarizedExperiment.rds is helpful for DEG analysis.
# R can directly read its output and lead to DEG analysis with DSS.

# ========== Step 1. Load SummarizedExperiment and convert it to DESeqDataSet ==========
convert_se_to_dds <- function () {
  se <- readRDS("./rnaseq_output/star_salmon/salmon.merged.gene.SummarizedExperiment.rds")
  # se_transcript <- readRDS("./rnaseq_output/star_salmon/salmon.merged.transcript.SummarizedExperiment.rds")

  # "SRX17589504:"Liver#314, "SRX17589506";Primary#87, "SRX17589513";Primary#314, "SRX17589516";Liver#87
  colData(se)$tissue <- factor(c("Liver", "Primary", "Primary", "Liver"))
  colData(se)$individual <- factor(c("#314", "#87", "#314", "#87"))

  { # DESeqDataSet requires assay(se)$counts to be numeric (and integer), so convert the mode of this assay to numeric
    raw_counts <- assay(se, "counts")
    # Force it into a numeric matrix
    numeric_counts <- as.matrix(raw_counts)
    # Force it into a integer matrix
    numeric_counts <- round(numeric_counts)
    # Verify the mode is now numeric
    mode(numeric_counts) 
    # Replace the counts assay
    assays(se)$counts <- numeric_counts
  }
  # Including individual differences in the design formula
  DESeqDataSet(se, design = ~ individual + tissue)
}

# ========== Step 2. Extract DEGs from DDS of step 1 ==========
extract_deg <- function(dds) {
  dds <- DESeq(dds)
  # Extracting results (difference between Liver and Primary)
  res <- results(dds, contrast = c("tissue", "Liver", "Primary"))
  
  # Convert results to a data frame and extract significant DEGs (p < 0.05)
  # (use dplyr functions to construct target data frame efficiently)
  deg_results <- as.data.frame(res) %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    mutate(gene_id = rownames(.)) # add gene_id column from row name
  # LFC is log2(treat / control) (in this analysis, Liver / Primary)

  # Need col.names to combine with DMRs info by gene_id
  write.csv(deg_results, "DEGs.csv")
}

{ # main
  dds <- convert_se_to_dds()
  extract_deg(dds)
  # Output to DEGs.csv
}
