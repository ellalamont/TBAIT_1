# Use ggforestplot to make a forest plot but with the underlying data from Bob's Forest Plot
# E. Lamont
# 2/2/26

source("Import_DEG_sets.R")
source("Import_GeneSets.R")

# install.packages("tidytext")
library(tidytext)
devtools::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)

# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = "white", color = NA)) + 
  theme(legend.position = "right",legend.text=element_text(size=9),
        legend.title = element_text(size = 10),
        # legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=9, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=8), 
        plot.subtitle = element_text(size=10), 
        plot.margin = margin(10, 10, 10, 20)# ,
  )

facet_themes <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                      strip.text = element_text(size = 10))

###########################################################
########### SAVE DATA FROM BOB'S FOREST PLOTS #############

CASS1.vs.CASS0_ForestData <- plotGeneSetForest(file = list_dfs_2$Rep1_CASS1.ComparedTo.CASS0,
                                 geneSets = allGeneSetList$LipidBiosynthesisPathways_LRF_GeneSets,
                                 main = "CASS1 vs CASS0 (Replicate 1 only)",
                                 min.genes.per.set = 3,
                                 max.show = 80, # There will only be 71 because of the 3 threshold above
                                 text.cex = 1.1, pt.cex = 1.25, lwd = 3.5)
# write.csv(Sputum_data,
#           file = "Figures/ForestPlots/EllaNew/Sputum_iModulonsForest.csv")



# Split the N_mean_SD column 
CASS1.vs.CASS0_ForestData_2 <- CASS1.vs.CASS0_ForestData %>%
  separate(N_Mean_SD, into = c("N", "mean_sd"), sep = ", ") %>%
  mutate(Mean = str_extract(mean_sd, "-?\\d+\\.\\d+") %>% as.numeric(),
         SD = str_extract(mean_sd, "\\(.*?\\)") %>% str_remove_all("[()]") %>% as.numeric) %>% 
  select(-mean_sd)

# Get the single CI value for error bars
CASS1.vs.CASS0_ForestData_3 <- CASS1.vs.CASS0_ForestData_2 %>%
  mutate(
    CI_clean = str_remove_all(CI, "[()]"),
    CI_lower = as.numeric(str_extract(CI_clean, "^-?\\d*\\.?\\d+")),
    CI_upper = as.numeric(str_extract(CI_clean, "(?<=, )-?\\d*\\.?\\d+")),
    CI2 = abs(Mean - CI_lower)
  )


###########################################################
################### GGFORESTPLOT PACKAGE ##################

# https://nightingalehealth.github.io/ggforestplot/articles/ggforestplot.html


CASS1.vs.CASS0_Forest <- CASS1.vs.CASS0_ForestData_3 %>% 
  ggforestplot::forestplot(name = Label, estimate = Mean, pvalue = P_Value, 
                           se = CI2/1.96) + # Need to divide by 1.96 because I am giving it the true CI, not the SE.
  # scale_color_manual(values = my_fav_colors) + 
  labs(title = "CASS1 vs CASS0 (Replicate 1 only)", x = "Log2Fold Change") + 
  my_plot_themes 
# test_graph$layers[[1]]$aes_params$odd <- "#00000000" # Run this code to remove the stripes (https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes)
CASS1.vs.CASS0_Forest



