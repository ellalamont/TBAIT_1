# Make PCAs
# E. Lamont
# 2/2/26

source("BatchCorrection_CombatSeq.R") # my_tpm_bc

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9))



###########################################################
##################### PCA ALL TPM_BC ######################

# Keep genes with >5 tpm (should be using counts!) in at least 50% of samples 
# 2/2/26: Doing this, not sure if necessary
keep <- rowSums(my_tpm_bc > 5) >= 0.5 * ncol(my_tpm_bc)
my_tpm_bc2 <- my_tpm_bc[keep, ] # now only 4445 genes (instead of 4499)

# Transform the data
my_tpm_bc_t <- as.data.frame(t(my_tpm_bc2))

# Remove columns that are all zero so the scale works for prcomp
# 2/2/26: Not doing this because I removed everything with <5 counts above
# my_tpm_bc_t2 <- my_tpm_bc_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_bc_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 63.6% of variance
summary_PCA[2,1] # PC2 explains 4.3% of variance
summary_PCA[3,1] # PC3 explains 4.1% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_pipeSummary, by = "SampleID", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, shape = Run, color = CASS_Pos, text = Run)) + 
  geom_point(size = 5, alpha = 0.8, stroke = 0.8) +
  # scale_color_manual(values = c("1" = "#D81B60", "2" = "#1E88E5", "3" = "#FFC107", "4" = "#004D40")) +  
  # scale_shape_manual(values = c("0" = 21, "1" = 24, "H37Rv" = 23)) + 
  labs(title = "PCA All TBAIT Samples + H37Rv",
       subtitle = "RawReads -> Batch Correction -> TPM -> <5 removed",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig
ggplotly(PCA_fig)
# ggsave(PCA_fig,
#        file = paste0("GoodSamples_tpmf_txnCov60_v1.pdf"),
#        path = "Figures/PCA",
#        width = 10, height = 6, units = "in")



###########################################################
###################### PCA ALL TPM ########################

# Keep genes with >5 tpm (should be using counts!) in at least 50% of samples 
# 2/2/26: Doing this, not sure if necessary
keep <- rowSums(my_TPM > 5) >= 0.5 * ncol(my_TPM)
my_TPM2 <- my_TPM[keep, ] # now only 4445 genes (instead of 4499)

# Transform the data
my_tpm_t <- as.data.frame(t(my_TPM2))

# Remove columns that are all zero so the scale works for prcomp
# 2/2/26: Not doing this because I removed everything with <5 counts above
# my_tpm_bc_t2 <- my_tpm_bc_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 14.2% of variance
summary_PCA[2,1] # PC2 explains 11.0% of variance
summary_PCA[3,1] # PC3 explains 8.0% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_pipeSummary, by = "SampleID", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, shape = Run, color = CASS_Pos, text = Run)) + 
  geom_point(size = 5, alpha = 0.8, stroke = 0.8) +
  # scale_color_manual(values = c("1" = "#D81B60", "2" = "#1E88E5", "3" = "#FFC107", "4" = "#004D40")) +  
  # scale_shape_manual(values = c("0" = 21, "1" = 24, "H37Rv" = 23)) + 
  labs(title = "PCA All TBAIT Samples + H37Rv",
       subtitle = "TPM -> <5 removed",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig
ggplotly(PCA_fig)
# ggsave(PCA_fig,
#        file = paste0("GoodSamples_tpmf_txnCov60_v1.pdf"),
#        path = "Figures/PCA",
#        width = 10, height = 6, units = "in")


###########################################################
################ PCA ALL TPM TBAIT_Run_1 ##################

# Keep genes with >5 tpm (should be using counts!) in at least 50% of samples 
# 2/2/26: Doing this, not sure if necessary
keep <- rowSums(TBAIT_Run_1_TPM > 5) >= 0.5 * ncol(TBAIT_Run_1_TPM)
my_TPM2 <- TBAIT_Run_1_TPM[keep, ] # now only 4445 genes (instead of 4499)

my_TPM2 <- my_TPM2 %>% 
  tibble::remove_rownames() %>%
  column_to_rownames("Gene")

# Transform the data
my_tpm_t <- as.data.frame(t(my_TPM2))

# Remove columns that are all zero so the scale works for prcomp
# 2/2/26: Not doing this because I removed everything with <5 counts above
# my_tpm_bc_t2 <- my_tpm_bc_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] 
summary_PCA[2,1] 
summary_PCA[3,1] 

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_pipeSummary, by = "SampleID", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, shape = as.character(Replicate), fill = PID)) + 
  geom_point(size = 5, alpha = 0.8, stroke = 0.8) +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c("1" = 21, "2" = 24, "3" = 23)) + 
  labs(title = "PCA TBAIT_Run_1",
       subtitle = "TPM -> <5 removed",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  # guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  my_plot_themes + 
  guides(fill = guide_legend(override.aes = list(shape = 21, alpha = 1, stroke = 0))) # legend fix!!
PCA_fig
ggplotly(PCA_fig)
# ggsave(PCA_fig,
#        file = paste0("TBAIT_Run_1_TPM_v1.pdf"),
#        path = "Figures/PCA",
#        width = 8, height = 6, units = "in")


###########################################################
############ PCA ALL TPM 20250508_LRF_pH_trim50 ###########

# Keep genes with >5 tpm (should be using counts!) in at least 50% of samples 
# 2/2/26: Doing this, not sure if necessary
keep <- rowSums(X20250508_LRF_pH_trim50_TPM > 5) >= 0.5 * ncol(X20250508_LRF_pH_trim50_TPM)
my_TPM2 <- X20250508_LRF_pH_trim50_TPM[keep, ] # now only 4445 genes (instead of 4499)

my_TPM2 <- my_TPM2 %>% 
  tibble::remove_rownames() %>%
  column_to_rownames("Gene")

# Transform the data
my_tpm_t <- as.data.frame(t(my_TPM2))

# Remove columns that are all zero so the scale works for prcomp
# 2/2/26: Not doing this because I removed everything with <5 counts above
# my_tpm_bc_t2 <- my_tpm_bc_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] 
summary_PCA[2,1] 
summary_PCA[3,1] 

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_pipeSummary, by = "SampleID", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, shape = as.character(Replicate), fill = PID)) + 
  geom_point(size = 5, alpha = 0.8, stroke = 0.8) +
  # scale_fill_viridis_d() +
  scale_shape_manual(values = c("1" = 21, "2" = 24, "3" = 23)) + 
  labs(title = "PCA 20250508_LRF_pH_trim50_TPM",
       subtitle = "TPM -> <5 removed",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  # guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  my_plot_themes + 
  guides(fill = guide_legend(override.aes = list(shape = 21, alpha = 1, stroke = 0))) # legend fix!!
PCA_fig
# ggplotly(PCA_fig)
# ggsave(PCA_fig,
#        file = paste0("X20250508_LRF_pH_trim50_TPM_v1.pdf"),
#        path = "Figures/PCA",
#        width = 8, height = 6, units = "in")

###########################################################
############ PCA ALL TPM 20250508_LRF_pH_trim50 ###########

# Keep genes with >5 tpm (should be using counts!) in at least 50% of samples 
# 2/2/26: Doing this, not sure if necessary
keep <- rowSums(X2025_10_LF_LA_EL_TPM > 5) >= 0.5 * ncol(X2025_10_LF_LA_EL_TPM)
my_TPM2 <- X2025_10_LF_LA_EL_TPM[keep, ] # now only 4445 genes (instead of 4499)

my_TPM2 <- my_TPM2 %>% 
  tibble::remove_rownames() %>%
  column_to_rownames("Gene")

# Transform the data
my_tpm_t <- as.data.frame(t(my_TPM2))

# Remove columns that are all zero so the scale works for prcomp
# 2/2/26: Not doing this because I removed everything with <5 counts above
# my_tpm_bc_t2 <- my_tpm_bc_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] 
summary_PCA[2,1] 
summary_PCA[3,1] 

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_pipeSummary, by = "SampleID", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, shape = as.character(Replicate), fill = PID)) + 
  geom_point(size = 5, alpha = 0.8, stroke = 0.8) +
  # scale_fill_viridis_d() +
  scale_shape_manual(values = c("1" = 21, "2" = 24, "3" = 23)) + 
  labs(title = "PCA X2025_10_LF_LA_EL_TPM",
       subtitle = "TPM -> <5 removed",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  # guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  my_plot_themes + 
  guides(fill = guide_legend(override.aes = list(shape = 21, alpha = 1, stroke = 0))) # legend fix!!
PCA_fig
# ggplotly(PCA_fig)
ggsave(PCA_fig,
       file = paste0("X2025_10_LF_LA_EL_TPM_v1.pdf"),
       path = "Figures/PCA",
       width = 8, height = 6, units = "in")

