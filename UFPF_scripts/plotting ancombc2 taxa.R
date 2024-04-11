# plotting UFPF anombc2 data 
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)  
library(camcorder)

# running at the species level 
ancom <- readRDS("UFPF/ANCOMBC2/ancombc2 species.rds")

res_prim = ancom$res
res_pair = ancom$res_pair      # 233 species examined 

sig_taxa <- res_pair %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_"))))    # 1 sig taxa significant 

# plot significant taxon only 
standard_errors <- sig_taxa[, c(1, 5:7)]

sig_taxa <- sig_taxa[1:4]

# Renaming columns a less friendy way that works 
colnames(sig_taxa)[colnames(sig_taxa) == "lfc_Diagnosis2IBD"] <- "IBD vs Control"
colnames(sig_taxa)[colnames(sig_taxa) == "lfc_Diagnosis2PD"] <- "PD vs Control"
colnames(sig_taxa)[colnames(sig_taxa) == "lfc_Diagnosis2PD_Diagnosis2IBD"] <- "PD vs IBD"

colnames(standard_errors)[colnames(standard_errors) == "se_Diagnosis2IBD"] <- "IBD vs Control"
colnames(standard_errors)[colnames(standard_errors) == "se_Diagnosis2PD"] <- "PD vs Control"
colnames(standard_errors)[colnames(standard_errors) == "se_Diagnosis2PD_Diagnosis2IBD"] <- "PD vs IBD"

# Melt the data frame into long format for ggplot
sig_taxa_long <- melt(sig_taxa, id.vars = "taxon")
se_long <- melt(standard_errors, id.vars = "taxon")

sig_taxa_long <- sig_taxa_long %>%
  rename(LFC = value)

se_long <- se_long %>%
  rename(standard_error = value)

sig_taxa_long <- left_join(sig_taxa_long, se_long, by = c("taxon", "variable"))

# adding a column to denote significance 
sig_taxa_long <-  sig_taxa_long %>%
  mutate(significance = c(TRUE, FALSE, FALSE))

# for plot, pick colors for comparisons 
legend_colors <- list(
  "IBD vs Control" = "chartreuse3",
  "PD vs Control" = "dodgerblue3",
  "PD vs IBD" = "darkmagenta"
)

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 10, height = 6,
          dir = "UFPF/Figures/")

# plotting with significant stars
ggplot(sig_taxa_long, aes(x = LFC, y = taxon, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbarh(
    aes(xmin = LFC - standard_error, xmax = LFC + standard_error),
    position = position_dodge(0.9),
    height = 0.25,  # Adjust the height of the error bars
    size = 0.5,    # Adjust the size of the error bars
    color = "gray73"  # Set the color of the error bars
  ) +
  geom_text(
    aes(label = ifelse(significance, "*", "")),
    position = position_dodge(0.9),
    vjust = 0.4, # Adjust vertical position of stars
    size = 10      # Adjust size of stars
  ) +
  labs(
    title = "Species-Level Differential Abundance Analysis",
    x = "Log Fold Change",
    y = "Species",
    fill = "Comparison"
  ) +
  scale_fill_manual(
    values = legend_colors,
    name = "Pairwise Comparisons"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 18, color = "black", face = "italic"),
    axis.text.x = element_text(size = 16, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  )

ggsave("UFPF/Figures/ANCOMBC2 species.png", dpi = 600, units = "in",
       height = 8, width = 10)

gg_stop_recording()



# ---------------------------------------------------------------------
# plotting GENUS level data 
ancom <- readRDS("UFPF/ANCOMBC2/ancombc2 genus.rds")

res_prim = ancom$res
res_pair = ancom$res_pair      # 140 genera  examined 

sig_taxa <- res_pair %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_"))))    # 2 sig taxa

standard_errors <- sig_taxa[, c(1, 5:7)]

sig_taxa <- sig_taxa[1:4]

# Renaming columns a less friendy way that works 
colnames(sig_taxa)[colnames(sig_taxa) == "lfc_Diagnosis2IBD"] <- "IBD vs Control"
colnames(sig_taxa)[colnames(sig_taxa) == "lfc_Diagnosis2PD"] <- "PD vs Control"
colnames(sig_taxa)[colnames(sig_taxa) == "lfc_Diagnosis2PD_Diagnosis2IBD"] <- "PD vs IBD"

colnames(standard_errors)[colnames(standard_errors) == "se_Diagnosis2IBD"] <- "IBD vs Control"
colnames(standard_errors)[colnames(standard_errors) == "se_Diagnosis2PD"] <- "PD vs Control"
colnames(standard_errors)[colnames(standard_errors) == "se_Diagnosis2PD_Diagnosis2IBD"] <- "PD vs IBD"

# Melt the data frame into long format for ggplot
sig_taxa_long <- melt(sig_taxa, id.vars = "taxon")
se_long <- melt(standard_errors, id.vars = "taxon")

sig_taxa_long <- sig_taxa_long %>%
  rename(LFC = value)

se_long <- se_long %>%
  rename(standard_error = value)

sig_taxa_long <- left_join(sig_taxa_long, se_long, by = c("taxon", "variable"))

# adding a column to denote significance 
sig_taxa_long <-  sig_taxa_long %>%
  mutate(significance = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))

# for plot, pick colors for comparisons 
legend_colors <- list(
  "IBD vs Control" = "chartreuse3",
  "PD vs Control" = "dodgerblue3",
  "PD vs IBD" = "darkmagenta"
)

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 12, height = 8,
          dir = "UFPF/Figures/")

# plotting with significant stars
ggplot(sig_taxa_long, aes(x = LFC, y = taxon, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbarh(
    aes(xmin = LFC - standard_error, xmax = LFC + standard_error),
    position = position_dodge(0.9),
    height = 0.25,  # Adjust the height of the error bars
    size = 0.5,    # Adjust the size of the error bars
    color = "gray73"  # Set the color of the error bars
  ) +
  geom_text(
    aes(label = ifelse(significance, "*", "")),
    position = position_dodge(0.9),
    vjust = 0.4, # Adjust vertical position of stars
    size = 10      # Adjust size of stars
  ) +
  labs(
    title = "Genus-level Differential Abundance Analysis",
    x = "Log Fold Change",
    y = "Genus",
    fill = "Comparison"
  ) +
  scale_fill_manual(
    values = legend_colors,
    name = "Pairwise Comparisons"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black", face = "italic"),
    axis.text.x = element_text(size = 16, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  )

ggsave("UFPF/Figures/ANCOMBC2 genus.png", dpi = 600, units = "in",
       height = 8, width = 10)

gg_stop_recording()

