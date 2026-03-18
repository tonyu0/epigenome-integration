# ========== Load libraries ==========
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DSS", "bsseq"), ask=F)
library(DSS)
library(bsseq)

# ----- Memo for how to deal with output data from bismark ----- 
# Data in methylseq_output/bismark/methylation_calls/methylation_calls is raw data, which is massive and containing every data
# (such like: CHG_OB_SRX17589483_trimmed_bismark_bt2.txt.gz)
# Data in methylseq_output/bismark/methylation_calls/methylation_coverage is summarized data, which is useful for subsequent R analysis.
# (such like: SRX17589483_trimmed_bismark_bt2.bismark.cov.gz)
# BAM file in methylseq_output/bismark/alignments can be also readable in R, but it is too massive to use.

# there is only two replicates (low-replicate) RRBS data in this analysis. 
# DSS is appropriate for this analysis (uses a Bayesian hierarchical model to shrink variance estimates)

# ========== Step 1. Load coverage files from bismark to BSseq data ==========
load_cov_to_bsseq <- function() {
  liver_files <- c(
    "./methylseq_output/bismark/methylation_calls/methylation_coverage/SRX17589493_trimmed_bismark_bt2.bismark.cov.gz",
    "./methylseq_output/bismark/methylation_calls/methylation_coverage/SRX17589500_trimmed_bismark_bt2.bismark.cov.gz"
    )
  primary_files <- c(
    "./methylseq_output/bismark/methylation_calls/methylation_coverage/SRX17589483_trimmed_bismark_bt2.bismark.cov.gz", 
    "./methylseq_output/bismark/methylation_calls/methylation_coverage/SRX17589490_trimmed_bismark_bt2.bismark.cov.gz"
    )
  # "SRX17589483":Primary#87, "SRX17589490":Primary#314, "SRX17589493":Liver#87, "SRX17589500":Liver#314
  sample_names <- c("Liver#87", "Liver#314", "Primary#87", "Primary#314")

  # Bismark .cov format: chr, start, end, meth%, count_meth, count_unmeth
  # DSS needs: chr, pos, N (total), X (methylated)
  read_dss <- function(f) {
    dat <- read.table(f, header=FALSE)
    data.frame(chr=dat$V1, pos=dat$V2, N=dat$V5+dat$V6, X=dat$V5)
  }
  
  # Load all samples
  liver_data <- lapply(liver_files, read_dss)
  primary_data <- lapply(primary_files, read_dss)
  makeBSseqData(c(liver_data, primary_data), sample_names)
}

# ========== Step 2. Call DMRs ==========

call_dmrs <- function(bsseq_data) {
  # Perform DML test (Differential Methylated Loci)
  dml_test <- DMLtest(bsseq_data, group1=c("Liver#87", "Liver#314"), group2=c("Primary#87", "Primary#314"), smoothing=TRUE)
  
  # Call DMRs (Differential Methylated Regions)
  # delta=0.1 means at least 10% difference in methylation
  callDMR(dml_test, delta=0.1, p.threshold=0.05, minlen=50, minCG=3, dis.merge=100)
}

# ========== Step 3. Output DMRs ==========
output_dmrs <- function(dmrs) {
  # DSS output columns: chr, start, end, length, nCG, meanMethy1, meanMethy2, diff.Methy, areaStat
  bed_output <- data.frame(
    chrom = dmrs$chr,
    chromStart = dmrs$start - 1, # BED file is 0-based
    chromEnd = dmrs$end,
    name = paste0("DMR_", seq_len(nrow(dmrs))),
    score = dmrs$diff.Methy # Using methylation difference as the score
  )
  
  # BED format: chr, start, end, name, score, strand
  write.table(bed_output, "DMRs.bed", sep="\t", quote=FALSE, row.names=T, col.names=T)
}

{ # main
  setwd("~/epigenome-integration")
  paste("Move to", getwd())
  
  bsseq_data <- load_cov_to_bsseq()
  dmrs <- call_dmrs(bsseq_data)
  output_dmrs(dmrs)
  # Output DMRs.bed
}

