# plotting HMP2 data combined (stacked bars)

# reading in HMP taxonomic data 
HMP_phylum <- readRDS("HMP2_Payami/Old/HMP2_phylum_relab.rds")
HMP_genus <- readRDS("HMP2_Payami/Old/HMP2_genus_relab.rds")
diagnosis <- readRDS("HMP2_Payami/HMP2 IBD diagnosis dataframe.rds")

# add in Sample column
HMP_phylum2 <- cbind(data.frame(Sample = rownames(HMP_phylum)), HMP_phylum)
rownames(HMP_phylum2) <- NULL
HMP_genus2 <- cbind(data.frame(Sample = rownames(HMP_genus)), HMP_genus)
rownames(HMP_genus2) <- NULL


# Step 2: Add the "Cohort" column from the diagnosis dataframe
HMP_phylum2$Cohort <- diagnosis$diagnosis2[match(HMP_phylum2$Sample, diagnosis$Sample)]
HMP_genus2$Cohort <- diagnosis$diagnosis2[match(HMP_genus2$Sample, diagnosis$Sample)]

# keeping only those sames that match the diagnosis dataframe (have been age filtered)
HMP_phylum2 <- HMP_phylum2[complete.cases(HMP_phylum2$Cohort), ]
HMP_genus2 <- HMP_genus2[complete.cases(HMP_genus2$Cohort), ]


# Convert to long format
phylum_long <- pivot_longer(HMP_phylum2, 
                            cols = -c(Sample, Cohort),
                            names_to = "Phylum",
                            values_to = "Abundance")

genus_long <- pivot_longer(HMP_genus2, 
                           cols = -c(Sample, Cohort),
                           names_to = "Genus",
                           values_to = "Abundance")

# convert to numeric
phylum_long <- phylum_long %>%
  mutate(Abundance = as.numeric(Abundance))

genus_long <- genus_long %>%
  mutate(Abundance = as.numeric(Abundance))


# Group by cohort and phylum, and calculate mean abundance
phylum_long_mean <- phylum_long %>%
  group_by(Cohort, Phylum) %>%
  summarize(mean_abundance = mean(Abundance))

genus_long_mean <- genus_long %>%
  group_by(Cohort, Genus) %>%
  summarize(mean_abundance = mean(Abundance))


library(RColorBrewer)
# Define the number of colors needed for each palette
pairedCount <- 8
set3Count <- 12

# Get the colors from the two palettes
pairedColors <- brewer.pal(pairedCount, "Paired")
set3Colors <- brewer.pal(set3Count, "Set3")

# Combine the colors into one vector
allColors <- c(set3Colors, pairedColors)

phylum_colors <- c("Acidobacteria" = "chocolate1", "Actinobacteria" = "cornflowerblue", 
                   "Bacteroidetes" = "aquamarine", "Deinococcus_Thermus" = "azure4", 
                   "Firmicutes" = "darkorchid", "Fusobacteria" = "darkgoldenrod", "Proteobacteria" =
                     "deeppink", "Spirochaetes" = "seagreen", "Synergistetes" = "lightpink", 
                   "Ternicutes" = "red1", "Verrucomicrobia" = "chartreuse1")

phylum_long_mean$Cohort <- factor(phylum_long_mean$Cohort, levels = c("nonIBD", "IBD"))

# Create stacked bar graph with the average cohort abundances 
ggplot(phylum_long_mean, aes(x = Cohort, y = mean_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "HMP2 Phylum-level Relative Abundances",
       x = "Cohort", y = "Relative Abundance (%)",
       fill = "Phylum") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    plot.title = element_text(size =21, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black")) +
  scale_fill_manual(values = phylum_colors) +
  scale_y_continuous(limits = c(0, 100))

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 phylum stacked bars.png", dpi = 600, units = "in",
       height = 7, width = 11)


# depict individual sample by sample 
phylum_long$Cohort <- factor(phylum_long$Cohort, levels = c("nonIBD", "IBD"))

# left to right 
ggplot(phylum_long, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ Cohort, scales = "free_x", space = "free_x") +
  labs(title = "HMP2 Phylum-level Relative Abundances",
       x = "Sample", y = "Relative Abundance (%)",
       fill = "Phylum") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
  theme(
    plot.title = element_text(size =21, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black"))+
  scale_fill_manual(values = phylum_colors) 
  #geom_segment(aes(x = 0, y = -10, xend = 1, yend = -10), color = "black", size = 1.2) +
  #annotate("text", x = 0.5, y = -15, label = "Cohorts", size = 4, fontface = "bold", color = "black")

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/phylum individual samples stacked bars.png", dpi = 600, units = "in",
       height = 6.5, width = 12)




# Select the top 10 genera by total abundance across all cohorts
top_genera2 <- genus_long_mean %>%
  group_by(Genus) %>%
  summarize(Total = sum(mean_abundance)) %>%
  top_n(10, Total) %>%
  pull(Genus)

# Filter the data to only include the top 20 genera
genus_long_mean_top20 <- genus_long_mean %>%
  filter(Genus %in% top_genera2)

phylum_colors3 <- c("Akkermansia" = "chocolate1", "Alistipes" = "gold", 
                    "Bacteroides" = "red1", "Barnesiella" = "darkorchid", 
                    "Bifidobacterium" = "indianred", "Blautia" = "darkgoldenrod", "Clostridium" =
                      "deeppink", "Dialister" = "seagreen", "Escherichia" = "lightpink", 
                    "Eubacterium" = "green3", "Faecalibacterium" = "chartreuse1", 
                    "Odoribacter" = "lightpink3", "Oscillibacter" = "cornflowerblue", 
                    "Parabacteroides" = "aquamarine", "Prevotella" = "azure4", 
                    "Roseburia" = "mediumorchid1", "Ruminococcus" = "moccasin", 
                    "Subdoligranulum" = "salmon", "Sutterella" = "darkblue", 
                    "Veillonella" = "darkcyan") 

genus_long_mean_top20$Cohort <- factor(genus_long_mean_top20$Cohort, levels = c("nonIBD", "IBD"))

# Plot the stacked bar chart for the top 20 genera
ggplot(genus_long_mean_top20, aes(x = Cohort, y = mean_abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "HMP2 Genus-level Relative Abundances (TOP 10)",
       x = "Cohort", y = "Relative Abundance (%)",
       fill = "Genus") +
  theme_classic() +
  theme(
    plot.title = element_text(size =20, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black"))+
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = phylum_colors3) 

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 genus stacked bars top 10.png", dpi = 600, units = "in",
       height = 7, width = 11)



# Create stacked bar graph with each individual sample depicted 
select_top_genera <- function(Genus, n = 10) {
  top_genera <- names(sort(colSums(Genus), decreasing = TRUE))[1:n]
  return(Genus[, top_genera2, drop = FALSE])
}

top_indiv_genera <- genus_long %>%
  group_by(Genus) %>%
  summarize(Total = sum(Abundance)) %>%
  top_n(10, Total) %>%
  pull(Genus)

# Filter the data to only include the top 10 genera
genus_cohorts_df_long_top10 <- genus_long %>%
  filter(Genus %in% top_indiv_genera)

# Calculate abundance for "Other" for each sample
other_abundance <- genus_cohorts_df_long_top10 %>%
  group_by(Sample, Cohort) %>%
  summarize(Total = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Genus = "Other", Abundance = 100 - Total)

# Add "Other" rows to the data frame
genus_cohorts_df_long_top10 <- bind_rows(genus_cohorts_df_long_top10, other_abundance)

genus_cohorts_df_long_top10$Genus <- factor(
  genus_cohorts_df_long_top10$Genus, 
  levels = c("Other", top_indiv_genera)
)

genus_cohorts_df_long_top10$Cohort <- factor(genus_cohorts_df_long_top10$Cohort, levels = c("nonIBD", "IBD"))

phylum_colors3 <- c("Other" = "seashell2", "Akkermansia" = "chocolate1", "Alistipes" = "gold", 
                    "Bacteroides" = "red1", "Barnesiella" = "darkorchid", 
                    "Bifidobacterium" = "indianred", "Blautia" = "darkgoldenrod", "Clostridium" =
                      "deeppink", "Dialister" = "seagreen", "Escherichia" = "dodgerblue", 
                    "Eubacterium" = "green3", "Faecalibacterium" = "chartreuse1", 
                    "Odoribacter" = "lightpink3", "Oscillibacter" = "cornflowerblue", 
                    "Parabacteroides" = "aquamarine", "Prevotella" = "azure4", 
                    "Roseburia" = "mediumorchid1", "Ruminococcus" = "moccasin", 
                    "Subdoligranulum" = "lightpink", "Sutterella" = "darkblue", 
                    "Veillonella" = "darkcyan") 

# left to right 
ggplot(genus_cohorts_df_long_top10, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ Cohort, scales = "free_x", space = "free_x") +
  labs(title = "HMP2 Genus-level Relative Abundances",
       x = "Sample", y = "Relative Abundance (%)",
       fill = "Genus") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black")
  ) +
  scale_fill_manual(values = phylum_colors3)

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/all genus individual samples stacked bars.png", dpi = 600, units = "in",
       height = 6.5, width = 12)



# ----------------------------------------------------------------------------
# reading in Payami PD+Control taxonomic data 
rel_ab_phylum <- readRDS("HMP2_Payami/Payami_phylum_relab.rds")
rel_ab_genus <- readRDS("HMP2_Payami/Payami_genus_relab.rds")

# adjusting RA values at each taxonomic level so they are all out of 100 
# Calculate the total for each row (excluding Sample column)
phylum_sum <- rowSums(rel_ab_phylum)
rel_ab_phylum <- rel_ab_phylum / phylum_sum * 100

genus_sum <- rowSums(rel_ab_genus)
rel_ab_genus <- rel_ab_genus / genus_sum * 100


# subset by cohort 
rel_ab_phylum2 <- cbind(
  data.frame(Sample_ID = rownames(rel_ab_phylum),
             cohort = ifelse(grepl("^DP\\d+|^SP\\d+", rownames(rel_ab_phylum)), "PD",
                             ifelse(grepl("^DC\\d+", rownames(rel_ab_phylum)), "Control", NA))),
  rel_ab_phylum
)

rel_ab_genus2 <- cbind(
  data.frame(Sample_ID = rownames(rel_ab_genus),
             cohort = ifelse(grepl("^DP\\d+|^SP\\d+", rownames(rel_ab_genus)), "PD",
                             ifelse(grepl("^DC\\d+", rownames(rel_ab_genus)), "Control", NA))),
  rel_ab_genus
)


# Convert to long format
phylum_cohorts_df_long <- pivot_longer(rel_ab_phylum2, 
                                       cols = -c(Sample_ID, cohort),
                                       names_to = "Phylum",
                                       values_to = "Abundance")

genus_cohorts_df_long <- pivot_longer(rel_ab_genus2, 
                                      cols = -c(Sample_ID, cohort),
                                      names_to = "Genus",
                                      values_to = "Abundance")

# Group by cohort and phylum, and calculate mean abundance
phylum_cohorts_df_long_mean <- phylum_cohorts_df_long %>%
  group_by(cohort, Phylum) %>%
  summarize(mean_abundance = mean(Abundance))

genus_cohorts_df_long_mean <- genus_cohorts_df_long %>%
  group_by(cohort, Genus) %>%
  summarize(mean_abundance = mean(Abundance))

# Define the number of colors needed for each palette
pairedCount <- 8
set3Count <- 12

# Get the colors from the two palettes
pairedColors <- brewer.pal(pairedCount, "Paired")
set3Colors <- brewer.pal(set3Count, "Set3")

# Combine the colors into one vector
allColors <- c(set3Colors, pairedColors)

phylum_colors2 <- c("Candidatus_Melainabacteria" = "chocolate1", "Actinobacteria" = "cornflowerblue", 
                    "Bacteroidetes" = "aquamarine", "Cyanobacteria" = "azure4", 
                    "Firmicutes" = "darkorchid", "Fusobacteria" = "darkgoldenrod", "Proteobacteria" =
                      "deeppink", "Synergistetes" = "lightpink", 
                    "Lentisphaerae" = "red1", "Verrucomicrobia" = "chartreuse1")

# Create stacked bar graph with the average cohort abundances 
ggplot(phylum_cohorts_df_long_mean, aes(x = cohort, y = mean_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "UAB Phylum-level Relative Abundances",
       x = "Cohort", y = "Relative Abundance (%)",
       fill = "Phylum") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    plot.title = element_text(size =20, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black")) +
  scale_fill_manual(values = phylum_colors2) 

ggsave("HMP2_Payami/Figures/UAB PDvsControl phylum stacked bars.png", dpi = 600, units = "in",
       height = 7, width = 11)



# depict individual sample by sample 
# left to right 
ggplot(phylum_cohorts_df_long, aes(x = Sample_ID, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ cohort, scales = "free_x", space = "free_x") +
  labs(title = "Wallen Phylum-level Relative Abundances",
       x = "Sample", y = "Relative Abundance (%)",
       fill = "Phylum") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
  theme(
    plot.title = element_text(size =21, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black"))+
  scale_fill_manual(values = phylum_colors) 
  #geom_segment(aes(x = 0, y = -10, xend = 1, yend = -10), color = "black", size = 1.2) +
  #annotate("text", x = 0.5, y = -15, label = "Cohorts", size = 4, fontface = "bold", color = "black")

ggsave("HMP2_Payami/Figures/Wallen phylum individual samples stacked bars.png", dpi = 600, units = "in",
       height = 6.5, width = 12)






# Select the top 20 genera by total abundance across all cohorts
top_genera <- genus_cohorts_df_long_mean %>%
  group_by(Genus) %>%
  summarize(Total = sum(mean_abundance)) %>%
  top_n(10, Total) %>%
  pull(Genus)

# Filter the data to only include the top 20 genera
genus_cohorts_df_long_mean_top20 <- genus_cohorts_df_long_mean %>%
  filter(Genus %in% top_genera)

phylum_colors4 <- c("Akkermansia" = "chocolate1", "Alistipes" = "gold", "Anaerostipes" = "hotpink",
                    "Bacteroides" = "red1", "Collinsella" = "lightgoldenrod1", "Dorea" = "darkorange4",  
                    "Bifidobacterium" = "darkgoldenrod", "Blautia" = "indianred",
                    "Escherichia" = "lightpink", "Eubacterium" = "green3",
                    "Faecalibacterium" = "chartreuse1", "Firmicutes_unclassified" = "forestgreen",
                    "Fusicatenibacter" = "darkturquoise", "Lachnospiraceae_unclassified" = "magenta", 
                    "Parabacteroides" = "aquamarine", "Prevotella" = "azure4", 
                    "Roseburia" = "mediumorchid1", "Ruminococcus" = "moccasin", 
                    "Ruminococcaceae_unclassified" = "honeydew", "Steptococcus" = "lightslateblue") 


# Plot the stacked bar chart for the top 20 genera
ggplot(genus_cohorts_df_long_mean_top20, aes(x = cohort, y = mean_abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Wallen Genus-level Relative Abundances (TOP 10)",
       x = "Cohort", y = "Relative Abundance (%)",
       fill = "Genus") +
  theme_classic() +
  theme(
    plot.title = element_text(size =20, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black"))+
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = phylum_colors4) 

ggsave("HMP2_Payami/Figures/UAB genus stacked bars top 10.png", dpi = 600, units = "in",
       height = 7, width = 11)


# Create stacked bar graph with each individual sample depicted 
select_top_genera <- function(Genus, n = 10) {
  top_genera <- names(sort(colSums(Genus), decreasing = TRUE))[1:n]
  return(Genus[, top_genera, drop = FALSE])
}

top_indiv_genera <- genus_cohorts_df_long %>%
  group_by(Genus) %>%
  summarize(Total = sum(Abundance)) %>%
  top_n(10, Total) %>%
  pull(Genus)

# Filter the data to only include the top 10 genera
genus_cohorts_df_long_top10 <- genus_cohorts_df_long %>%
  filter(Genus %in% top_indiv_genera)

# Calculate abundance for "Other" for each sample
other_abundance <- genus_cohorts_df_long_top10 %>%
  group_by(Sample_ID, cohort) %>%
  summarize(Total = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Genus = "Other", Abundance = 100 - Total)

# Add "Other" rows to the data frame
genus_cohorts_df_long_top10 <- bind_rows(genus_cohorts_df_long_top10, other_abundance)

genus_cohorts_df_long_top10$Genus <- factor(
  genus_cohorts_df_long_top10$Genus, 
  levels = c("Other", top_indiv_genera)
)

phylum_colors5 <- c("Other" = "seashell2", "Akkermansia" = "chocolate1", "Alistipes" = "gold", 
                    "Bacteroides" = "red1", "Barnesiella" = "darkorchid", 
                    "Bifidobacterium" = "dodgerblue", "Blautia" = "darkgoldenrod", "Clostridium" =
                      "deeppink", "Collinsella", "darkcyan", "Dialister" = "seagreen", "Escherichia" = "dodgerblue", 
                    "Eubacterium" = "green3", "Faecalibacterium" = "chartreuse1", 
                    "Lachnospiraceae_unclassified" = "hotpink",
                    "Odoribacter" = "lightpink3", "Oscillibacter" = "cornflowerblue", 
                    "Parabacteroides" = "aquamarine", "Prevotella" = "azure4", 
                    "Roseburia" = "mediumorchid1", "Ruminococcus" = "moccasin", 
                    "Subdoligranulum" = "lightpink", "Sutterella" = "darkblue") 

# left to right 
ggplot(genus_cohorts_df_long_top10, aes(x = Sample_ID, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ cohort, scales = "free_x", space = "free_x") +
  labs(title = "Wallen Genus-level Relative Abundances",
       x = "Sample", y = "Relative Abundance (%)",
       fill = "Genus") +
  theme_classic() +
  #scale_y_continuous(limits = c(0, 100)) +
  theme(
    plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black")
  ) +
  scale_fill_manual(values = phylum_colors5)

ggsave("HMP2_Payami/Figures/Wallen all genus individual samples stacked bars.png", dpi = 600, units = "in",
       height = 6.5, width = 12)
