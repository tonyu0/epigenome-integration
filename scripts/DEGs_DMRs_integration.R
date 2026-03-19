# ========== Load libraries ==========
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"), ask=F)
# For region operations
library(GenomicRanges)
library(rtracklayer)
# For DMRs 
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(org.Hs.eg.db)

# ========== Step 1. Define CGI shore ==========
define_cgi_shores <- function() {
  cgi <- import("Human_CpGisland.bed")
  
  # Use 'flank' to get regions outside the island boundaries (CpG island shores)
  # <-up [5' -> 3'] down->
  shores_up <- flank(cgi, width = 2000, start = T) # start = T means up side flank, F means down side flank
  shores_down <- flank(cgi, width = 2000, start = F)
  reduce(c(shores_up, shores_down)) # Combine and merge overlaps
  # flank: start~end
  #   start=T: start-width ~ start-1
  #       both=T: start-width ~ start+width
  #.  start=F: end+1 ~ end+width
  #       both=T: end-width ~ end+width
}

# ========== Step 2. Annotate DMRs in CGI shores with nearby genes using ChIPseeker =========

annotate_DMRs_in_shores <- function (cgi_shores) {
  dmrs <- import("DMRs_chr.bed")
  # Find DMRs that overlap with CGI shores
  # Ensure that the chr_column of dmrs and chr_column of cgi_shores have the same format (they should look like "chr1" or "chrX" respectively).
  dmrs_in_shores <- subsetByOverlaps(dmrs, cgi_shores)
  
  # Annotate DMRs with nearest TSS/Gene
  dmrs_anno <- annotatePeak(dmrs_in_shores, 
                           tssRegion = c(-3000, 3000), # what range of distance is called neighbors (explicitly setting default value)
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                           annoDb = "org.Hs.eg.db") # For Gene Symbols
  
  as.data.frame(dmrs_anno)
}

# ========== Step 3. Integrate with DEGs ==========
integrate_DEGs_with_DMRs <- function(dmrs_anno) {
  degs <- read.csv("DEGs.csv")
  # Merge by gene_id
  final_genes <- merge(dmrs_anno, degs, by.x = "ENSEMBL", by.y = "gene_id")
  
  # Filter for biologically interesting patterns 
  # (e.g., Hypermethylation in Liver vs Primary with Downregulated expression)
  interest_genes <- subset(final_genes, padj < 0.05 & abs(log2FoldChange) > 1)
  # ENSG00000196549: Liver vs Primary. methylation score is high (means liver methylation is higher than Primary, and log2FoldChange is lower.)
  write.csv(interest_genes, "DMRs_regulated_genes.csv")
}

{ # main
  cgi_shores <- define_cgi_shores()
  dmrs_anno <- annotate_DMRs_in_shores(cgi_shores)
  integrate_DEGs_with_DMRs(dmrs_anno)
}


