if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer", "RColorBrewer", "TxDb.Hsapiens.UCSC.hg38.knownGene"), ask = F)

library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# ========== 1. Scatter plot with ggplot ==========
draw_expression_methlylation_scatter <- function() { # Scatter plot of expression and methylation of extracted gene
  dmr_regulated_genes <-  read.csv("DMRs_regulated_genes.csv")
  
  png("expression_methylation_scatter.png", width = 1000, height = 1200, res = 150)
  p <- ggplot(data = dmr_regulated_genes, aes(x = score, y = log2FoldChange, color = SYMBOL)) +
    geom_point() +
    labs(title = "Expression level vs. Methylation level (Liver vs Primary)", 
         x = "Methylation level (diff %)", y = "Expression level (Log2 Fold Change)", color = "Gene symbol")
  print(p) # this is required when executed from function
  dev.off()
}

# ========== 2. Visualize tracks of gene model (annotation) with Gviz ==========
draw_gene_with_gviz <- function() {
  # DMR: 155079868,155080372
  # Gene: 155079932,155183678,
  target_chr <- "chr3" 
  target_start <- 155078000
  target_end   <- 155084000
  target_gene <- "ENSG00000196549"
  gen <- "hg38"
  
  axis_track <- GenomeAxisTrack(chromosome=target_chr)
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  gene_track <- GeneRegionTrack(txdb, name = "Gene",
                                genome = gen, chromosome = target_chr,
                  		from = target_start, to = target_end, 
				filter = list(gene_id = target_gene), # only target_gene to draw
				collapseTranscripts = "meta", # all genes(ENSG) and transcripts(ENST) to 1 line
				showId = TRUE,
    				geneSymbol = TRUE,
    				showExonId = FALSE,
                                transcriptAnnotation = "symbol", background.title = "brown")
  
  # ========== CGI track ==========
  cgi <- import("Human_CpGisland.bed")
  shores_up <- flank(cgi, width = 2000, start = T) # start = T means up side flank, F means down side flank
  shores_down <- flank(cgi, width = 2000, start = F)
  cgi_shores <- reduce(c(shores_up, shores_down)) # Combine and merge overlaps
  cgi_track   <- AnnotationTrack(cgi, name = "CGI", fill = "green", chromosome=target_chr)
  shore_track <- AnnotationTrack(cgi_shores, name = "Shore", fill = "lightgreen", chromosome=target_chr)
  # ========== DMR track ==========
  dmrs <- import("DMRs_chr.bed")
  dmrs_track   <- AnnotationTrack(dmrs, name = "DMR", fill = "red", alpha = 0.8, chromosome=target_chr)
  
  # Use data only from #314
#options(ucscChromosomeNames=FALSE) # bismark bedGraph somehow output error
  # ========== RNA-seq track (use bigWig from star_salmon) ==========
#liver_expression_track <- DataTrack(range = "rnaseq_output/star_salmon/bigwig/SRX17589504.forward.bigWig", name = "Liver RNA", 
#                                      type = "gradient", col = "darkred", ylim = c(0, 1), chromosome=target_chr)
#primary_expression_track <- DataTrack(range = "rnaseq_output/star_salmon/bigwig/SRX17589513.forward.bigWig" , name = "Primary RNA", 
#                                        type = "gradient", col = "darkblue", ylim = c(0, 1), chromosome=target_chr)
  # ========== Methylation track (use bedGraph from bismark) ==========
#liver_methylation_track <- DataTrack(range = "methylseq_output/bismark/methylation_calls/bedGraph/SRX17589500_trimmed_bismark_bt2.bedGraph.gz", name = "Liver Meth", 
#                                       type = "histogram", fill = "salmon", col = "salmon", chromosome=target_chr)
#primary_methylation_track <- DataTrack(range = "methylseq_output/bismark/methylation_calls/bedGraph/SRX17589490_trimmed_bismark_bt2.bedGraph.gz", name = "Primary Meth", 
#                                         type = "histogram", fill = "skyblue", col = "skyblue", chromosome=target_chr)
  # ========== Plot ==========
  png("MME_integration_plot.png", width = 1920, height = 1080, res = 150)
  plotTracks(list(axis_track, 
#                  liver_expression_track, liver_methylation_track, 
#                  primary_expression_track, primary_methylation_track, 
                  dmrs_track, shore_track, cgi_track, 
                  gene_track),
             from = target_start, to = target_end, chromosome=target_chr,
             main = "Epigenetic Regulation of MME (ENSG00000196549)", transcriptAnnotation = "symbol")
  dev.off()
}

{ # main
  print("========== Drawing expression_methylation_scatter.png ==========")
  draw_expression_methlylation_scatter()
  print("========== Drawing MME_integration_plot.png ==========")
  draw_gene_with_gviz()
}
