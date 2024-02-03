# plotting UFPF anombc2 data 
library(ggplot2)
library(reshape2)
library(paletteer)
library(dplyr)
library(tidyr)
library(openxlsx)
library(knitr)
library(kableExtra)
library(ggpubr)  
library(camcorder)


# running at the species level 
ancom <- readRDS("UFPF/ANCOMBC2/ancombc2 species.rds")

res_prim = ancom$res
res_pair = ancom$res_pair      # 233 species examined 

sig_taxa <- res_pair %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_"))))    # 1 sig taxa significant 



# creating a neater table with all results from ancombc2 
table <- res_pair[, c(1:7, 11:16)]

table <- table %>%
  rename(
    species = taxon,
    lfc_IBD = lfc_Diagnosis2IBD,
    lfc_PD = lfc_Diagnosis2PD,
    lfc_PDvsIBD = lfc_Diagnosis2PD_Diagnosis2IBD, 
    se_IBD = se_Diagnosis2IBD,
    se_PD = se_Diagnosis2PD,
    se_PDvsIBD = se_Diagnosis2PD_Diagnosis2IBD,
    pval_IBD = p_Diagnosis2IBD,
    pval_PD = p_Diagnosis2PD,
    pval_PDvsIBD = p_Diagnosis2PD_Diagnosis2IBD,
    adj_pval_IBD = q_Diagnosis2IBD,
    adj_pval_PD = q_Diagnosis2PD,
    adj_pval_PDvsIBD = q_Diagnosis2PD_Diagnosis2IBD
  )

# Create nice table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF ANCOMBC2 Species Taxonomic Output") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

save_kable(table_formatted, "UFPF/ANCOMBC2/table all species ancombc2 output.html")

# Save as Excel
write.xlsx(as.data.frame(table), "UFPF/ANCOMBC2/table all species ancombc2 output.xlsx")




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
          dir = "With Second Batch/Figures/")

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



# ---------------------------------------------------
# PLOT ALL 
# just keep the taxon and lfc columns
ancom_res_pair <- res_pair

# Renaming columns a less friendy way that works 
colnames(ancom_res_pair)[colnames(ancom_res_pair) == "lfc_Diagnosis2IBD"] <- "IBD vs Control"
colnames(ancom_res_pair)[colnames(ancom_res_pair) == "lfc_Diagnosis2PD"] <- "PD vs Control"
colnames(ancom_res_pair)[colnames(ancom_res_pair) == "lfc_Diagnosis2PD_Diagnosis2IBD"] <- "PD vs IBD"

colnames(ancom_res_pair)[colnames(ancom_res_pair) == "se_Diagnosis2IBD"] <- "IBD vs Control"
colnames(ancom_res_pair)[colnames(ancom_res_pair) == "se_Diagnosis2PD"] <- "PD vs Control"
colnames(ancom_res_pair)[colnames(ancom_res_pair) == "se_Diagnosis2PD_Diagnosis2IBD"] <- "PD vs IBD"

taxa <- ancom_res_pair[1:4]
se <- ancom_res_pair[, c(1, 5:7)]

# Melt the data frame into long format for ggplot
taxa_long <- melt(taxa, id.vars = "taxon")
se_long <- melt(se, id.vars = "taxon")

taxa_long <- taxa_long %>%
  rename(LFC = value)

se_long <- se_long %>%
  rename(standard_error = value)

taxa_long <- left_join(taxa_long, se_long, by = c("taxon", "variable"))

# for plot, pick colors for comparisons 
legend_colors <- list(
  "IBD vs Control" = "chartreuse3",
  "PD vs Control" = "dodgerblue3",
  "PD vs IBD" = "darkmagenta"
)

ggplot(taxa_long, aes(x = LFC, y = taxon, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbarh(
    aes(xmin = LFC - standard_error, xmax = LFC + standard_error),
    position = position_dodge(0.9),
    height = 0.25,  # Adjust the height of the error bars
    size = 0.5,    # Adjust the size of the error bars
    color = "gray73"  # Set the color of the error bars
  ) +
  labs(
    title = "Pairwise Differential Abundance Analysis",
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
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black", face = "italic"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  )

ggsave("UFPF/Figures/ANCOMBC2 all taxa.png", dpi = 600, units = "in",
       height = 28, width = 15)




# ----------------------------------


ancom_res_pair <- res_pair

colnames(ancom_res_pair)[colnames(ancom_res_pair) == "lfc_Diagnosis2IBD"] <- "IBD vs Control"
colnames(ancom_res_pair)[colnames(ancom_res_pair) == "lfc_Diagnosis2PD"] <- "PD vs Control"
colnames(ancom_res_pair)[colnames(ancom_res_pair) == "lfc_Diagnosis2PD_Diagnosis2IBD"] <- "PD vs IBD"

# saving enriched vs depleted species from ancombc2 regardless of significance 
# Create separate dataframes for each comparison
df_ibd_vs_control <- ancom_res_pair %>% select(taxon, 'IBD vs Control')
df_pd_vs_control <- ancom_res_pair %>% select(taxon, 'PD vs Control')
df_pd_vs_ibd <- ancom_res_pair %>% select(taxon, 'PD vs IBD')

# Create separate columns for positive and negative fold changes
df_ibd_vs_control_pos <- df_ibd_vs_control[df_ibd_vs_control$`IBD vs Control` > 0, ]
df_ibd_vs_control_neg <- df_ibd_vs_control[df_ibd_vs_control$`IBD vs Control` < 0, ]

df_pd_vs_control_pos <- df_pd_vs_control[df_pd_vs_control$`PD vs Control` > 0, ]
df_pd_vs_control_neg <- df_pd_vs_control[df_pd_vs_control$`PD vs Control` < 0, ]

df_pd_vs_ibd_pos <- df_pd_vs_ibd[df_pd_vs_ibd$`PD vs IBD` > 0, ]
df_pd_vs_ibd_neg <- df_pd_vs_ibd[df_pd_vs_ibd$`PD vs IBD` < 0, ]


df_ibd_vs_control_pos <- df_ibd_vs_control_pos %>%
  arrange(desc(abs(`IBD vs Control`)))
df_ibd_vs_control_neg <- df_ibd_vs_control_neg %>%
  arrange(desc(abs(`IBD vs Control`)))
df_pd_vs_control_pos <- df_pd_vs_control_pos %>%
  arrange(desc(abs(`PD vs Control`)))
df_pd_vs_control_neg <- df_pd_vs_control_neg %>%
  arrange(desc(abs(`PD vs Control`)))
df_pd_vs_ibd_pos <- df_pd_vs_ibd_pos %>%
  arrange(desc(abs(`PD vs IBD`)))
df_pd_vs_ibd_neg <- df_pd_vs_ibd_neg %>%
  arrange(desc(abs(`PD vs IBD`)))


write.csv(df_ibd_vs_control_pos, file = "UFPF/ANCOMBC2/All species ibd_vs_control_pos.csv", row.names = FALSE)
write.csv(df_ibd_vs_control_neg, file = "UFPF/ANCOMBC2/All species ibd_vs_control_neg.csv", row.names = FALSE)

write.csv(df_pd_vs_control_pos, file = "UFPF/ANCOMBC2/All species pd_vs_control_pos.csv", row.names = FALSE)
write.csv(df_pd_vs_control_neg, file = "UFPF/ANCOMBC2/All species pd_vs_control_neg.csv", row.names = FALSE)

write.csv(df_pd_vs_ibd_pos, file = "UFPF/ANCOMBC2/All species pd_vs_ibd_pos.csv", row.names = FALSE)
write.csv(df_pd_vs_ibd_neg, file = "UFPF/ANCOMBC2/All species pd_vs_ibd_neg.csv", row.names = FALSE)




# Making tables 
# running at the genus level 
ancom <- readRDS("UFPF/ANCOMBC2/ancombc2 genus.rds")
# fixed effects: age + sex + diagnosis 
res_pair = ancom$res_pair             # 140 genera 

# creating a neater table with all results from ancombc2 
table <- res_pair[, c(1:7, 11:16)]

table <- table %>%
  rename(
    genus = taxon,
    lfc_IBD = lfc_Diagnosis2IBD,
    lfc_PD = lfc_Diagnosis2PD,
    lfc_PDvsIBD = lfc_Diagnosis2PD_Diagnosis2IBD, 
    se_IBD = se_Diagnosis2IBD,
    se_PD = se_Diagnosis2PD,
    se_PDvsIBD = se_Diagnosis2PD_Diagnosis2IBD,
    pval_IBD = p_Diagnosis2IBD,
    pval_PD = p_Diagnosis2PD,
    pval_PDvsIBD = p_Diagnosis2PD_Diagnosis2IBD,
    adj_pval_IBD = q_Diagnosis2IBD,
    adj_pval_PD = q_Diagnosis2PD,
    adj_pval_PDvsIBD = q_Diagnosis2PD_Diagnosis2IBD
  )

# Create nice table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF ANCOMBC2 Genus Taxonomic Output") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

save_kable(table_formatted, "UFPF/ANCOMBC2/ANCOMBC2 tables/table genus taxonomic output.html")

# Save as Excel
write.xlsx(as.data.frame(table), "UFPF/ANCOMBC2/ANCOMBC2 tables/table genus taxonomic output.xlsx")



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
          dir = "With Second Batch/Figures/")

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



# creating tables of enriched vs depleted genera
ancom_res_pair <- res_pair

colnames(ancom_res_pair)[colnames(ancom_res_pair) == "lfc_Diagnosis2IBD"] <- "IBD vs Control"
colnames(ancom_res_pair)[colnames(ancom_res_pair) == "lfc_Diagnosis2PD"] <- "PD vs Control"
colnames(ancom_res_pair)[colnames(ancom_res_pair) == "lfc_Diagnosis2PD_Diagnosis2IBD"] <- "PD vs IBD"

# saving enriched vs depleted species from ancombc2 regardless of significance 
# Create separate dataframes for each comparison
df_ibd_vs_control <- ancom_res_pair %>% select(taxon, 'IBD vs Control')
df_pd_vs_control <- ancom_res_pair %>% select(taxon, 'PD vs Control')
df_pd_vs_ibd <- ancom_res_pair %>% select(taxon, 'PD vs IBD')

# Create separate columns for positive and negative fold changes
df_ibd_vs_control_pos <- df_ibd_vs_control[df_ibd_vs_control$`IBD vs Control` > 0, ]
df_ibd_vs_control_neg <- df_ibd_vs_control[df_ibd_vs_control$`IBD vs Control` < 0, ]

df_pd_vs_control_pos <- df_pd_vs_control[df_pd_vs_control$`PD vs Control` > 0, ]
df_pd_vs_control_neg <- df_pd_vs_control[df_pd_vs_control$`PD vs Control` < 0, ]

df_pd_vs_ibd_pos <- df_pd_vs_ibd[df_pd_vs_ibd$`PD vs IBD` > 0, ]
df_pd_vs_ibd_neg <- df_pd_vs_ibd[df_pd_vs_ibd$`PD vs IBD` < 0, ]


df_ibd_vs_control_pos <- df_ibd_vs_control_pos %>%
  arrange(desc(abs(`IBD vs Control`)))
df_ibd_vs_control_neg <- df_ibd_vs_control_neg %>%
  arrange(desc(abs(`IBD vs Control`)))
df_pd_vs_control_pos <- df_pd_vs_control_pos %>%
  arrange(desc(abs(`PD vs Control`)))
df_pd_vs_control_neg <- df_pd_vs_control_neg %>%
  arrange(desc(abs(`PD vs Control`)))
df_pd_vs_ibd_pos <- df_pd_vs_ibd_pos %>%
  arrange(desc(abs(`PD vs IBD`)))
df_pd_vs_ibd_neg <- df_pd_vs_ibd_neg %>%
  arrange(desc(abs(`PD vs IBD`)))


write.csv(df_ibd_vs_control_pos, file = "UFPF/ANCOMBC2/All genus ibd_vs_control_pos.csv", row.names = FALSE)
write.csv(df_ibd_vs_control_neg, file = "UFPF/ANCOMBC2/All genus ibd_vs_control_neg.csv", row.names = FALSE)

write.csv(df_pd_vs_control_pos, file = "UFPF/ANCOMBC2/All genus pd_vs_control_pos.csv", row.names = FALSE)
write.csv(df_pd_vs_control_neg, file = "UFPF/ANCOMBC2/All genus pd_vs_control_neg.csv", row.names = FALSE)

write.csv(df_pd_vs_ibd_pos, file = "UFPF/ANCOMBC2/All genus pd_vs_ibd_pos.csv", row.names = FALSE)
write.csv(df_pd_vs_ibd_neg, file = "UFPF/ANCOMBC2/All genus pd_vs_ibd_neg.csv", row.names = FALSE)
