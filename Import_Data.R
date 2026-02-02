# Import data from All TBAT runs 
# E Lamont
# 2/2/26

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
library(ggprism) # for add_pvalue()
library(rstatix) # for adjust_pvalue
library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
library(pheatmap)
# library(dendextend) # May need this for looking at pheatmap clustering
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly
library(edgeR) # for cpm
library(sva) # For ComBat_seq batch correction
library(stringr)
library(readxl) # To import excel files as dataframes

# DuffyTools
# library(devtools)
# install_github("robertdouglasmorrison/DuffyTools")
# library(DuffyTools)
# install_github("robertdouglasmorrison/DuffyNGS")
# BiocManager::install("robertdouglasmorrison/DuffyTools")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biobase")


# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############

# TBAIT_Run_1
TBAIT_Run_1_pipeSummary <- read.csv("Data/TBAIT_Run_1/Pipeline.Summary.Details.csv") %>% select(-X) %>%
  mutate(Run = "TBAIT_Run_1")

# TBAIT_Seq
TBAIT_Seq_pipeSummary <- read.csv("Data/TBAIT_Seq/Pipeline.Summary.Details.csv") %>% select(-X) %>%
  mutate(Run = "TBAIT_Seq")

# 2025_10_LF_LA_EL
X2025_10_LF_LA_EL_pipeSummary <- read.csv("Data/2025_10_LF_LA_EL/Pipeline.Summary.Details.csv") %>% select(-X) %>%
  mutate(Run = "2025_10_LF_LA_EL")

# 20250508_LRF_pH_trim50
X20250508_LRF_pH_trim50_pipeSummary <- read.csv("Data/20250508_LRF_pH_trim50/Pipeline.Summary.Details.csv") %>% select(-X) %>%
  mutate(Run = "20250508_LRF_pH_trim50")

# Merge the pipeSummaries
All_pipeSummary <- merge(TBAIT_Run_1_pipeSummary, TBAIT_Seq_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, X2025_10_LF_LA_EL_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, X20250508_LRF_pH_trim50_pipeSummary, all = T)


###########################################################
###################### IMPORT METADATA ####################

All_metadata <- read.csv("Data/Metadata/TBAIT_seq_metadata_Clean.csv")

# Merge metadata with pipeSummary and only keep sampleIDs that are in metadata
my_pipeSummary <- All_pipeSummary %>% inner_join(All_metadata, by = "SampleID") # 227 samples

# Remove samples with low N_Genomic
my_pipeSummary <- my_pipeSummary %>% filter(CASS_Pos != "Low_N_Genomic") # 186 samples


###########################################################
##################### IMPORT RAW READS ####################

# TBAIT_Run_1
TBAIT_Run_1_RawReads <- read.csv("Data/TBAIT_Run_1/Mtb.Expression.Gene.Data.ReadsM.csv") %>% rename("Gene" = "X") %>%
  select(any_of(c("Gene", my_pipeSummary$SampleID)))

# TBAIT_Seq
TBAIT_Seq_RawReads <- read.csv("Data/TBAIT_Seq/Mtb.Expression.Gene.Data.ReadsM.csv") %>% rename("Gene" = "X") %>%
  select(any_of(c("Gene", my_pipeSummary$SampleID)))

# 2025_10_LF_LA_EL
X2025_10_LF_LA_EL_RawReads <- read.csv("Data/2025_10_LF_LA_EL/Mtb.Expression.Gene.Data.ReadsM.csv") %>% rename("Gene" = "X") %>%
  select(any_of(c("Gene", my_pipeSummary$SampleID)))

# 20250508_LRF_pH_trim50
X20250508_LRF_pH_trim50_RawReads <- read.csv("Data/20250508_LRF_pH_trim50/Mtb.Expression.Gene.Data.ReadsM.csv") %>% rename("Gene" = "X") %>%
  select(any_of(c("Gene", my_pipeSummary$SampleID)))

# Merge the Raw Reads
my_RawReads <- merge(TBAIT_Run_1_RawReads, TBAIT_Seq_RawReads, all = T)
my_RawReads <- merge(my_RawReads, X2025_10_LF_LA_EL_RawReads, all = T)
my_RawReads <- merge(my_RawReads, X20250508_LRF_pH_trim50_RawReads, all = T) # 187 variables


###########################################################
##################### IMPORT BOB'S TPM ####################

# TBAIT_Run_1
TBAIT_Run_1_TPM <- read.csv("Data/TBAIT_Run_1/Mtb.Expression.Gene.Data.TPM.csv") %>% rename("Gene" = "X") %>%
  select(any_of(c("Gene", my_pipeSummary$SampleID)))

# TBAIT_Seq
TBAIT_Seq_TPM <- read.csv("Data/TBAIT_Seq/Mtb.Expression.Gene.Data.TPM.csv") %>% rename("Gene" = "X") %>%
  select(any_of(c("Gene", my_pipeSummary$SampleID)))

# 2025_10_LF_LA_EL
X2025_10_LF_LA_EL_TPM <- read.csv("Data/2025_10_LF_LA_EL/Mtb.Expression.Gene.Data.TPM.csv") %>% rename("Gene" = "X") %>%
  select(any_of(c("Gene", my_pipeSummary$SampleID)))

# 20250508_LRF_pH_trim50
X20250508_LRF_pH_trim50_TPM <- read.csv("Data/20250508_LRF_pH_trim50/Mtb.Expression.Gene.Data.TPM.csv") %>% rename("Gene" = "X") %>%
  select(any_of(c("Gene", my_pipeSummary$SampleID)))

# Merge the Raw Reads
my_TPM <- merge(TBAIT_Run_1_TPM, TBAIT_Seq_TPM, all = T)
my_TPM <- merge(my_TPM, X2025_10_LF_LA_EL_TPM, all = T)
my_TPM <- merge(my_TPM, X20250508_LRF_pH_trim50_TPM, all = T)

my_TPM <- my_TPM %>% column_to_rownames("Gene")










