# load packages to perform alpha diversity analysis 
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)

# reading .tsv file containing taxonomic classification relative abundances 
rel_ab <- readRDS("UFPF/Metaphlan output/rel_ab_cleaned.rds")

# filter at the phylum level 
rel_ab_phylum <- rel_ab[, grep("\\|p__[^|]*$", colnames(rel_ab))]
colnames(rel_ab_phylum) <- sub("^.*\\|p__", "", colnames(rel_ab_phylum))

# Subset the data by cohort 
PD_phylum <- rel_ab_phylum[grep("1N", rownames(rel_ab_phylum)), ]
Control_phylum <- rel_ab_phylum[grep("3N|3G", rownames(rel_ab_phylum)), ]
IBD_phylum <- rel_ab_phylum[grep("2G", rownames(rel_ab_phylum)), ]

PD_phylum_df <- data.frame(Sample = rownames(PD_phylum), PD_phylum)
Control_phylum_df <- data.frame(Sample = rownames(Control_phylum), Control_phylum)
IBD_phylum_df <- data.frame(Sample = rownames(IBD_phylum), IBD_phylum)
phylum_cohorts_df <- rbind(PD_phylum_df, Control_phylum_df, IBD_phylum_df)

# Convert to long format
phylum_cohorts_df_long <- pivot_longer(phylum_cohorts_df, 
                                       cols = -Sample,
                                       names_to = "Phylum",
                                       values_to = "Abundance")

# Extract cohort info from Sample column
phylum_cohorts_df_long$Cohort <- ifelse(grepl("1N", phylum_cohorts_df_long$Sample), "PD",
                                        ifelse(grepl("2G", phylum_cohorts_df_long$Sample), "IBD", "Control"))

# Group by cohort and phylum + calculate mean abundance
phylum_cohorts_df_long_mean <- phylum_cohorts_df_long %>%
  group_by(Cohort, Phylum) %>%
  summarize(mean_abundance = mean(Abundance))

# Define  number of colors needed for each palette
pairedCount <- 8
set3Count <- 12

# Get the colors from the two palettes
pairedColors <- brewer.pal(pairedCount, "Paired")
set3Colors <- brewer.pal(set3Count, "Set3")
allColors <- c(set3Colors, pairedColors)

# Create stacked bar graph with the average cohort abundances 
# HERE IT REMOVES ACTINOBACTERA FROM THE PD COHORT *************
ggplot(phylum_cohorts_df_long_mean, aes(x = Cohort, y = mean_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Phylum-level Relative Abundances by Cohort",
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
  scale_fill_manual(values = allColors) 

ggsave("UFPF/Figures/phylum stacked bars.png", dpi = 600, units = "in",
       height = 7, width = 11)


# have samples next to each other left to right 
# BUNCH OF WARNING MESSAGES FOR ROWS THAT IT REMOVES *********** 
ggplot(phylum_cohorts_df_long, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ Cohort, scales = "free_x", space = "free_x") +
  labs(title = "Phylum-level Relative Abundances",
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
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black"))+
  scale_fill_manual(values = allColors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) 
  #geom_segment(aes(x = 0, y = -10, xend = 1, yend = -10), color = "black", size = 1.2) +
  #annotate("text", x = 0.5, y = -15, label = "Cohorts", size = 4, fontface = "bold", color = "black")

ggsave("UFPF/Figures/phylum individual stacked bars.png", dpi = 600, units = "in",
       height = 6.5, width = 12)




# stacked bar charts at the genus level 
# --------------------------------------------------------------------------
# filter at the genus level 
rel_ab_genus <- rel_ab[, grep("\\|g__[^|]*$", colnames(rel_ab))]
colnames(rel_ab_genus) <- sub("^.*\\|g__", "", colnames(rel_ab_genus))

# Subset the data by cohort 
PD_genus <- rel_ab_genus[grep("1N", rownames(rel_ab_genus)), ]
Control_genus <- rel_ab_genus[grep("3N|3G", rownames(rel_ab_genus)), ]
IBD_genus <- rel_ab_genus[grep("2G", rownames(rel_ab_genus)), ]

# add a Sample column
PD_genus_df <- data.frame(Sample = rownames(PD_genus), PD_genus)
Control_genus_df <- data.frame(Sample = rownames(Control_genus), Control_genus)
IBD_genus_df <- data.frame(Sample = rownames(IBD_genus), IBD_genus)
genus_cohorts_df <- rbind(PD_genus_df, Control_genus_df, IBD_genus_df)

# Convert to long format
genus_cohorts_df_long <- pivot_longer(genus_cohorts_df, 
                                      cols = -Sample,
                                      names_to = "Genus",
                                      values_to = "Abundance")

# Extract cohort information from Sample column
genus_cohorts_df_long$Cohort <- ifelse(grepl("1N", genus_cohorts_df_long$Sample), "PD",
                                       ifelse(grepl("2G", genus_cohorts_df_long$Sample), "IBD", "Control"))

# Group by cohort and genus, and calculate mean abundance
genus_cohorts_df_long_mean <- genus_cohorts_df_long %>%
  group_by(Cohort, Genus) %>%
  summarize(mean_abundance = mean(Abundance))


# ---------------------------------------------------------------------------------------------------
# ONLY DEPICTING THE TOP  GENERA IN EACH COHORT 
# Select the top genera by total abundance across all cohorts
top_genera <- genus_cohorts_df_long_mean %>%
  group_by(Genus) %>%
  summarize(Total = sum(mean_abundance)) %>%
  top_n(10, Total) %>%
  pull(Genus)

# Filter the data to only include the top genera
genus_cohorts_df_long_mean_top10 <- genus_cohorts_df_long_mean %>%
  filter(Genus %in% top_genera)


# Identify non-top genera
non_top_genera <- genus_cohorts_df_long_mean_top10 %>%
  group_by(cohort) %>%
  summarize(mean_abundance = 100 - sum(mean_abundance)) %>%
  mutate(Genus = "Other")

# Combine the "Other" category with the existing data
genus_cohorts_df_long_mean_top10 <- rbind(genus_cohorts_df_long_mean_top10, non_top_genera)


# Plot the stacked bar chart for the top 20 genera
# No warnings I can see here ****
ggplot(genus_cohorts_df_long_mean_top10, aes(x = cohort, y = mean_abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Top Genus-level Relative Abundances by Cohort",
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
  scale_fill_manual(values = allColors) 

ggsave("UFPF/Figures/genus stacked bars.png", dpi = 600, units = "in",
       height = 7, width = 11)


# each individual sample depicted 
# Define a function to select the top  genera by abundance
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

# Define colors for each category
genus_colors <- c("Other" = "seashell2", "Porphyromonas" = "steelblue1", setNames(allColors, top_indiv_genera))


# left to right 
# doesn't appear to remove any *****
ggplot(genus_cohorts_df_long_top10, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ Cohort, scales = "free_x", space = "free_x") +
  labs(title = "Genus-level Relative Abundances",
       x = "Sample", y = "Relative Abundance (%)",
       fill = "Genus") +
  theme_classic() +
  #scale_y_continuous(limits = c(0, 100)) +
  theme(
    plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 11),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 7, angle =45, hjust = 1, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black")
  ) +
  scale_fill_manual(values = genus_colors)

ggsave("UFPF/Figures/genus individual stacked bars all.png", dpi = 600, units = "in",
       height = 6.5, width = 12)

