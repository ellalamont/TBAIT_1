# DESeq2 Cough 75 vs Cough 25
# E. Lamont 
# 2.10.26

# This might be better than Duffytools because I think there is a way for DESeq2 to incorporate replicates
# Also I think there is a way to do batch correction as well which we need to do for these samples

# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering


source("Import_data.R")


################################################
################ LOAD PACKAGES #################

# https://github.com/JulioLeonIncio/Tutorial-from-DEseq2-to-GSEA-in-R/blob/main/DEseq2_to_GSEA_JL.Rmd
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html
# https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-DESeq2/

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# pkgs <- c("DESeq2","sva","fgsea","clusterProfiler","GSEABase","tidyverse","pheatmap", "apeglm")
# for (p in pkgs) {
#   if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask=FALSE)
# }
library(DESeq2)
library(sva) # ComBat_seq
library(fgsea)
# library(clusterProfiler)
library(GSEABase)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(msigdbr)
library(apeglm)


################################################
############### BATCH CORRECTION ###############
# With combat-seq

source("BatchCorrection_CombatSeq.R")
# my_RawReads_bc


################################################
################## PRE-FILTER ##################

# Remove the H37Rv
my_RawReads_bc2 <- my_RawReads_bc %>% 
  column_to_rownames("X") %>% 
  dplyr::select(!contains("Rv"))

# Remove genes with <5 reads
keep <- rowSums(my_RawReads_bc2 > 5) >= 0.5 * ncol(my_RawReads_bc2)
my_RawReads_bc3 <- my_RawReads_bc2[keep, ] # now only 4439 genes (instead of 4499)

my_pipeSummary2 <- my_pipeSummary %>%
  filter(SampleID %in% colnames(my_RawReads_bc3)) %>%
  dplyr::select(SampleID, Replicate, PID, CASS_Pos, Cough_count, Cough_10, Cough_25, Cough_50, HHT, MainLineage) %>%
  filter(!is.na(Cough_25)) %>% # Need to remove any NA from the condition we will be comparing (cough25)
  column_to_rownames("SampleID")

# Now need to remove the NA Cough_25 from the RawReads
my_RawReads_bc4 <- my_RawReads_bc3 %>%
  dplyr::select(all_of(rownames(my_pipeSummary2)))


# Put the pipeSummary in the same order as the raw reads
my_pipeSummary3 <- my_pipeSummary2[colnames(my_RawReads_bc4), ]


################################################
############ GENERATE DESEQ2 OBJECT ############

# Ensure sample names line up
stopifnot(all(colnames(my_RawReads_bc4) == rownames(my_pipeSummary3)))
stopifnot(all(rownames(my_pipeSummary3) == colnames(my_RawReads_bc4)))

# Make the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = my_RawReads_bc4,
                              colData = my_pipeSummary3,
                              design = ~ Cough_25)
dds$Replicate
dds$PID # 93
colData(dds)

# Collapse the technical replicates (need to discuss if these are actually technical)
# https://bioinformatics.stackexchange.com/questions/18347/how-does-deseq2-collapsereplicates-function-work-on-read-counts-data
# DESeq2 is summing the technical replicates to collapse them I think...
ddsColl <- collapseReplicates(dds, dds$PID, renameCols = TRUE)
ddsColl$PID # 43
my_pipeSummary %>% filter(Type == "Clinical_Isolate") %>% filter(!is.na(Cough_25)) %>% pull(PID) %>% unique() # 43 here also!


################################################
################## RUN THE DEG #################

# Specify the reference level
ddsColl$Cough_25 <- relevel(ddsColl$Cough_25, ref = "25")

# Run the DEG
ddsColl_DEG <- DESeq(ddsColl)

# Get the results
res <- results(ddsColl_DEG)
res

# Optional: Log fold change shrinkage
# https://support.bioconductor.org/p/77461/
resultsNames(ddsColl_DEG)
# [1] "Intercept"         "Cough_25_75_vs_25"
resLFC <- lfcShrink(ddsColl_DEG, coef = "Cough_25_75_vs_25", type = "apeglm")
resLFC


################################################
############ SAVE AS NICE DATAFRAME ############

res_Cough75vs25_df <- as.data.frame(resLFC) %>% # Could use resLFC and see how it is different
  rownames_to_column("gene")

# Columns for Log2Fold > 1
res_Cough75vs25_df$DE1 <- ifelse(res_Cough75vs25_df$log2FoldChange < -1 & res_Cough75vs25_df$padj < 0.05, "significant down", ifelse(res_Cough75vs25_df$log2FoldChange > 1 & res_Cough75vs25_df$padj < 0.05, "significant up", "not significant"))
res_Cough75vs25_df$DE1 <- factor(res_Cough75vs25_df$DE1, levels = c("significant down", "not significant", "significant up"))
res_Cough75vs25_df$DE1_labels <- ifelse(res_Cough75vs25_df$DE1 != "not significant", res_Cough75vs25_df$gene, NA)

################################################
################# VOLCANO PLOT #################

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=10))

my_volcano <- res_Cough75vs25_df %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = DE1, label = DE1_labels, text = gene)) + # text is for plotly, could be GENE_ID
  geom_point(alpha = 0.7) + 
  labs(title = "DESeq2 Cough75 vs Cough25 Replicates") + 
  geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
  geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel 
  scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
plot_build <- ggplot_build(my_volcano)
y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
text_up <- res_Cough75vs25_df %>% filter(DE1 == "significant up") %>% nrow()
text_down <- res_Cough75vs25_df %>% filter(DE1 == "significant down") %>% nrow()
my_volcano_annotated <- my_volcano +
  annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
  annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
final_volcano <- my_volcano_annotated + my_plot_themes
final_volcano
# ggplotly(final_volcano)
ggsave(final_volcano,
       file = "Cough75vsCough25_DESeq2_v1.pdf",
       path = "Figures/Volcano/DESeq2",
       width = 7, height = 5, units = "in")



