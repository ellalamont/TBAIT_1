# Import Bob's differential gene expression data and the metagenesets data
# E. Lamont
# 2/2/26

source("Import_data.R")


###########################################################
################# IMPORT DESeq2 DATAFRAMES ################

parent_dir <- "Data/DESeq2_Results"
csv_files <- list.files(parent_dir, full.names = TRUE)

# Read all CSVs into a list of data.frames
DESeq2_dfs <- lapply(csv_files, read.csv)

# Add the names to the list of lists
names(DESeq2_dfs) <- tools::file_path_sans_ext(basename(csv_files))

# Also make a list of the names
DESeq2_df_names <- names(DESeq2_dfs)











# Data is coming from the Lenovo TBAIT_AllRuns

###########################################################
################### IMPORT BOB's DE DATA ##################

`Rep1_CASS1.ComparedTo.CASS0` <- read.delim("Data/DE_AllRuns/Rep1_only/Rep1_CASS0.vs.CASS1/CASS1_Rep1.MTb.Meta.JOINED.txt")
`Rep1_Cough50.ComparedTo.Cough49` <- read.delim("Data/DE_AllRuns/Rep1_only/Rep1_Cough50.vs.Cough49/C50_Rep1.MTb.Meta.JOINED.txt")
`Rep1_Cough75.ComparedTo.Cough25` <- read.delim("Data/DE_AllRuns/Rep1_only/Rep1_Cough75.vs.Cough25/C75_Rep1.MTb.Meta.JOINED.txt")
`Rep1_Cough90.ComparedTo.Cough10` <- read.delim("Data/DE_AllRuns/Rep1_only/Rep1_Cough90.vs.Cough10/C90_Rep1.MTb.Meta.JOINED.txt")

###########################################################
################ MAKE A LIST OF ALL DFs ###################
list_dfs <- list(`Rep1_CASS1.ComparedTo.CASS0`,
                 `Rep1_Cough50.ComparedTo.Cough49`, 
                 `Rep1_Cough75.ComparedTo.Cough25`, 
                 `Rep1_Cough90.ComparedTo.Cough10`)

# Make a list of the names
df_names <- c("Rep1_CASS1.ComparedTo.CASS0",
              "Rep1_Cough50.ComparedTo.Cough49", 
              "Rep1_Cough75.ComparedTo.Cough25", 
              "Rep1_Cough90.ComparedTo.Cough10")

# Give the df list the correct df names
names(list_dfs) <- df_names

###########################################################
############### ADD COLUMNS OF DE VALUES ##################

# 2/2/26: Using original p-values for all, and looking at log2fold 1 and 0.5

# Make a new list to hold dataframes with extra columns
list_dfs_2 <- list()

ordered_DE <- c("significant down", "not significant", "significant up")

# Add extra DE columns to each dataframe
for (i in 1:length(list_dfs)) {
  
  current_df <- list_dfs[[i]]
  current_df_name <- df_names[i]
  
  # Columns for Log2Fold > 1 and AVG_PVALUE < 0.05
  current_df$DE1 <- ifelse(current_df$LOG2FOLD < -1 & current_df$AVG_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 1 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE1 <- factor(current_df$DE1, levels = ordered_DE)
  current_df$DE1_labels <- ifelse(current_df$DE1 != "not significant", current_df$GENE_NAME, NA)
  
  # Columns for Log2Fold > 0.5 and AVG_PVALUE < 0.05
  current_df$DE0.5 <- ifelse(current_df$LOG2FOLD < -0.5 & current_df$AVG_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 0.5 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE0.5 <- factor(current_df$DE0.5, levels = ordered_DE)
  current_df$DE0.5_labels <- ifelse(current_df$DE0.5 != "not significant", current_df$GENE_NAME, NA)
  
  list_dfs_2[[current_df_name]] <- current_df
}

###########################################################
############# IMPORT BOB's METAGENESETS DATA ##############

`MetaGeneSets_Rep1_CASS1.ComparedTo.CASS0` <- read.delim("Data/DE_AllRuns/Rep1_only/Rep1_CASS0.vs.CASS1/CASS1_Rep1.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_Rep1_Cough90.ComparedTo.Cough10` <- read.delim("Data/DE_AllRuns/Rep1_only/Rep1_Cough90.vs.Cough10/C90_Rep1.MTb.MetaGeneSets.UP.txt")
