# Run Batch Correction with CombatSeq
# E. Lamont
# 2/2/26

source("Import_Data.R")
source("Function_CalculateTPM.R")

# https://github.com/zhangyuqing/ComBat-seq
# https://rnabio.org/module-03-expression/0003/06/02/Batch-Correction/




###########################################################
############# RUN COMBAT SEQ BATCH CORRECTION #############

count_matrix <- as.matrix(my_RawReads %>% column_to_rownames("Gene"))
# Ensure integer counts
mode(count_matrix) <- "integer"

# Reorder metadata to match column order in count_matrix
meta <- my_pipeSummary[match(colnames(count_matrix), my_pipeSummary$SampleID), ]
# Check alignment
all(meta$SampleID == colnames(count_matrix))  # should be TRUE

# # IF NOT TRUE RUN THESE:
# ## This will show the samples in count_matrix that don't match metadata
# # colnames(count_matrix)[!colnames(count_matrix) %in% my_pipeSummary$SampleID]
# ## And the opposite: metadata samples not in count_matrix
# # my_pipeSummary$SampleID[!my_pipeSummary$SampleID %in% colnames(count_matrix)]

# Extract batch (run) and condition (for checking later)
batch <- meta$Run
condition <- meta$CASS_Pos # Should condition be patient ID or CASS score?
# Tried condition being PID and it made the batch effect way worse on the PCA!



# Run ComBat-Seq
combat_counts <- ComBat_seq(
  count_matrix,
  batch = batch,
  group = NULL# , # optional, helps preserve biological signal
  # covar_mod = covar_mat # Would set group to NULL if using this
)


# # What if I want multiple biological variables?
# covar_mat <- cbind(meta$CASS_Pos, meta$PID)
# # Run ComBat-Seq
# combat_counts_MultVar <- ComBat_seq(
#   count_matrix,
#   batch = batch,
#   group = NULL, # optional, helps preserve biological signal
#   covar_mod = covar_mat # Would set group to NULL if using this
# )
# # Error in qr.default(design) : NA/NaN/Inf in foreign function call (arg 1)
# # Not sure why this isn't working

###########################################################
###################### CONVERT TO TPM #####################

my_RawReads_bc <- as.data.frame(combat_counts) %>% rownames_to_column("X")
my_tpm_bc <- CalculateTPM(my_RawReads_bc)
colSums(my_tpm_bc) # Just check it's all 1M










