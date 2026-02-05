# Make graphs showing distribution of lineages
# E. Lamont 
# 2/4/26

source("Import_data.R")

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",legend.text=element_text(size=8),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"), legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=10))


###########################################################
############## CASS: MAKE A STACKED BARPLOT ###############

my_lineages_Percents <- my_pipeSummary %>%
  filter(CASS_Pos != "H37Rv") %>%
  filter(MixedLineages2 == "No") %>%
  filter(Replicate == "1") %>%
  mutate(MainLineage = paste0("L", sub("^Lineage ([0-9]+).*", "\\1", LineageID))) %>% 
  count(CASS_Pos, MainLineage) %>%
  group_by(CASS_Pos) %>%
  mutate(Percent = round(n/sum(n) * 100, 1)) %>%
  ungroup()
my_lineages_Percents


Lineage_plot <- my_lineages_Percents %>%
  ggplot(aes(x = CASS_Pos, y = Percent, fill = MainLineage)) + 
  geom_col(color = "black", alpha = 0.5) + 
  # geom_text(aes(label = scales::percent(after_stat(prop), accuracy = 1)), stat = "count", position = position_fill(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("L2" = "dodgerblue2", "L3" = "purple3", "L4" = "red3", "L1" = "magenta2", "mixed" = "grey")) + 
  scale_y_continuous(limits = c(0,100.1), breaks = seq(0, 100.1, 10), expand = expansion(mult = c(0, 0))) +
  annotate("text", x = "0", y = 85, label = "25.9% L2 (n=14)", size = 3) +
  annotate("text", x = "0", y = 60, label = "44.4% L3 (n=24)", size = 3) +
  annotate("text", x = "0", y = 15, label = "27.8% L4 (n=15)", size = 3) +
  annotate("text", x = "1", y = 85, label = "25.9% L2 (n=7)", size = 3) +
  annotate("text", x = "1", y = 60, label = "25.9% L3 (n=7)", size = 3) +
  annotate("text", x = "1", y = 15, label = "48.1% L4 (n=13)", size = 3) +
  labs(x = "CASS positivity", y = "% of samples", fill = "Lineage",
       title = NULL) + 
  my_plot_themes
Lineage_plot
ggsave(Lineage_plot,
       file = "All_Lineages_Bar1.pdf",
       path = "Figures/Lineages",
       width = 4, height = 4, units = "in")

my_lineages_Percents2 <- my_pipeSummary %>%
  filter(CASS_Pos != "H37Rv") %>%
  filter(MixedLineages2 == "No") %>%
  mutate(MainLineage = paste0("L", sub("^Lineage ([0-9]+).*", "\\1", LineageID))) %>% 
  count(CASS_Pos, MainLineage) %>%
  group_by(MainLineage) %>%
  mutate(Percent = round(n/sum(n) * 100, 1)) %>%
  ungroup()
my_lineages_Percents2
Lineage_plot2 <- my_lineages_Percents2 %>%
  ggplot(aes(x = MainLineage, y = Percent, fill = CASS_Pos)) + 
  geom_col(color = "black") + 
  # geom_text(aes(label = scales::percent(after_stat(prop), accuracy = 1)), stat = "count", position = position_fill(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("0" = "grey", "1" = "cyan3")) + 
  scale_y_continuous(limits = c(0,100.1), breaks = seq(0, 100.1, 10), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "% of samples", fill = "Lineage",
       title = NULL) + 
  my_plot_themes
Lineage_plot2
ggsave(Lineage_plot2,
       file = "All_Lineages_Bar2.pdf",
       path = "Figures/Lineages",
       width = 5, height = 4, units = "in")


###########################################################
############# COUGH75: MAKE A STACKED BARPLOT #############

my_lineages_Percents <- my_pipeSummary %>%
  filter(CASS_Pos != "H37Rv") %>%
  filter(MixedLineages2 == "No") %>%
  filter(Replicate == "1") %>%
  mutate(MainLineage = paste0("L", sub("^Lineage ([0-9]+).*", "\\1", LineageID))) %>% 
  count(Cough_25, MainLineage) %>%
  group_by(Cough_25) %>%
  mutate(Percent = round(n/sum(n) * 100, 1)) %>%
  ungroup() %>%
  filter(!is.na(Cough_25))
my_lineages_Percents

Cough75_Lineage_plot <- my_lineages_Percents %>%
  ggplot(aes(x = Cough_25, y = Percent, fill = MainLineage)) + 
  geom_col(color = "black", alpha = 0.5) + 
  # geom_text(aes(label = scales::percent(after_stat(prop), accuracy = 1)), stat = "count", position = position_fill(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("L2" = "dodgerblue2", "L3" = "purple3", "L4" = "red3", "L1" = "magenta2", "mixed" = "grey")) + 
  scale_y_continuous(limits = c(0,100.1), breaks = seq(0, 100.1, 10), expand = expansion(mult = c(0, 0))) +
  annotate("text", x = "25", y = 90, label = "37.5% L2 (n=6)", size = 3) +
  annotate("text", x = "25", y = 50, label = "31.2% L3 (n=5)", size = 3) +
  annotate("text", x = "25", y = 15, label = "25% L4 (n=4)", size = 3) +
  annotate("text", x = "75", y = 90, label = "19.2% L2 (n=5)", size = 3) +
  annotate("text", x = "75", y = 50, label = "42.3% L3 (n=11)", size = 3) +
  annotate("text", x = "75", y = 15, label = "38.5% L4 (n=10)", size = 3) +
  labs(x = "Cough", y = "% of samples", fill = "Lineage",
       title = NULL) + 
  my_plot_themes
Cough75_Lineage_plot
ggsave(Cough75_Lineage_plot,
       file = "Cough75_Lineage_Bar1.pdf",
       path = "Figures/Lineages",
       width = 4, height = 4, units = "in")




