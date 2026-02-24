# DESeq2 PTLD normal vs restrictive
# E. Lamont 
# 2.6.26

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
# library(GSEABase)
# library(tidyverse)
# library(pheatmap)
# library(dplyr)
# library(msigdbr)
# library(apeglm)


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
  dplyr::select(SampleID, Replicate, PID, PTLD_12_mo) %>%
  replace_na(list(PTLD_12_mo = "unknown")) %>% # Replace the NA because the dds will break with NA present
  column_to_rownames("SampleID")

# Put the pipeSummary in the same order as the raw reads
my_pipeSummary3 <- my_pipeSummary2[colnames(my_RawReads_bc3), ]


################################################
############ GENERATE DESEQ2 OBJECT ############

# Ensure sample names line up
stopifnot(all(colnames(my_RawReads_bc3) == rownames(my_pipeSummary3)))
stopifnot(all(rownames(my_pipeSummary3) == colnames(my_RawReads_bc3)))

# Make the DESeqDataSet object (using all the data)
dds <- DESeqDataSetFromMatrix(countData = my_RawReads_bc3,
                              colData = my_pipeSummary3,
                              design = ~ PTLD_12_mo)
dds$Replicate
dds$PID # 180
# colData(dds)
design(dds) # ~PTLD_12_mo
resultsNames(dds) # 0...?
levels(dds$PTLD_12_mo)
# [1] "Both"        "Normal"      "Obstructive" "Restrictive" "unknown"   

# Collapse the technical replicates (need to discuss if these are actually technical)
# https://bioinformatics.stackexchange.com/questions/18347/how-does-deseq2-collapsereplicates-function-work-on-read-counts-data
# DESeq2 is summing the technical replicates to collapse them I think...
ddsColl <- collapseReplicates(dds, dds$PID, renameCols = TRUE)
ddsColl$PID # 83
# my_pipeSummary %>% filter(Type == "Clinical_Isolate") %>% pull(PID) %>% unique() # 83 here also!


################################################
################## RUN THE DEG #################

# Specify the reference level
ddsColl$PTLD_12_mo <- relevel(ddsColl$PTLD_12_mo, ref = "Normal")

# Run the DEG
ddsColl_DEG <- DESeq(ddsColl)

# Look at what comparisons there are
resultsNames(ddsColl_DEG)
# [1] "Intercept"                        "PTLD_12_mo_Both_vs_Normal"        "PTLD_12_mo_Obstructive_vs_Normal"
# [4] "PTLD_12_mo_Restrictive_vs_Normal" "PTLD_12_mo_unknown_vs_Normal"   

# See how many samples are in each condition
table(ddsColl_DEG$PTLD_12_mo)
#  Normal        Both Obstructive Restrictive     unknown 
# 22           3           2          14          42 

# Get the results
res <- results(ddsColl_DEG, name = "PTLD_12_mo_Restrictive_vs_Normal")
res


# Optional: Log fold change shrinkage
# https://support.bioconductor.org/p/77461/
# resultsNames(ddsColl_DEG)
# # [1] "Intercept"       "CASS_Pos_1_vs_0"
# resLFC <- lfcShrink(ddsColl_DEG, coef = "CASS_Pos_1_vs_0", type = "apeglm")
# resLFC



################################################
############ SAVE AS NICE DATAFRAME ############

res_PTLD.NormalvsRestrictive_df <- as.data.frame(res) %>% # Could use res or resLFC and see how it is different
  rownames_to_column("gene")

add_DE_columns_function <- function(my_df, log2fold_cutoff = 1, padj_cutoff = 0.05) {
  
  ## Adds columns to DEG dataframe results for plotting with the volcano plot function
  
  # Determine column names
  DE_col_name <- paste0("DE", log2fold_cutoff)
  label_col_name <- paste0(DE_col_name, "_labels")
  
  # Make the df with new columns
  my_df[[DE_col_name]] <- ifelse(my_df$log2FoldChange < -log2fold_cutoff & my_df$padj < padj_cutoff, "significant down",
                                 ifelse(my_df$log2FoldChange > log2fold_cutoff & my_df$padj < padj_cutoff, "significant up", "not significant"))
  
  # Make the DE significance a factor
  my_df[[DE_col_name]] <- factor(my_df[[DE_col_name]], levels = c("significant down", "not significant", "significant up"))
  
  # Add gene name column
  my_df[[label_col_name]] <- ifelse(my_df[[DE_col_name]] != "not significant", my_df$gene, NA)
  
  return(my_df)
}

res_PTLD.NormalvsRestrictive_df <- add_DE_columns_function(res_PTLD.NormalvsRestrictive_df, log2fold_cutoff = 1, padj_cutoff = 0.05)

# Save as csv
# write.csv(res_PTLD.NormalvsRestrictive_df, "Data/DESeq2_Results/res_PTLD.NormalvsRestrictive_df.csv")

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

my_volcano <- res_PTLD.NormalvsRestrictive_df %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = DE1, label = DE1_labels, text = gene)) + # text is for plotly, could be GENE_ID
  geom_point(alpha = 0.7) + 
  labs(title = "DESeq2 PTLD Restrictive vs Normal Replicates collapsed") + 
  geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
  geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel 
  scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
plot_build <- ggplot_build(my_volcano)
y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
text_up <- res_PTLD.NormalvsRestrictive_df %>% filter(DE1 == "significant up") %>% nrow()
text_down <- res_PTLD.NormalvsRestrictive_df %>% filter(DE1 == "significant down") %>% nrow()
my_volcano_annotated <- my_volcano +
  annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
  annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
final_volcano <- my_volcano_annotated + my_plot_themes
final_volcano
# ggplotly(final_volcano)
# ggsave(final_volcano,
#        file = "PTLD.RestrictiveVsNormal_DESeq2_v1.pdf",
#        path = "Figures/Volcano/DESeq2",
#        width = 7, height = 5, units = "in")



