library(ggplot2)
library(tidyverse)
library(tidyr)

# Plotting ANCOMBC2 data 
# SPECIES LEVEL DATA - Wallen PD 
ancom_PD1 <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD species ancombc2 output.rds")

res_prim_PD1 = ancom_PD1$res

sig_taxa_PD1 <- res_prim_PD1 %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))

sig_taxa_PD <- sig_taxa_PD1[, c(1, 3, 7, 19, 23)]

sig_taxa_PD_long <- sig_taxa_PD %>%
  pivot_longer(
    cols = starts_with("lfc_"),  
    values_to = "LFC"           
  )

colnames(sig_taxa_PD_long)[colnames(sig_taxa_PD_long) == "q_Case_statusPD"] <- "Adj P Value"
colnames(sig_taxa_PD_long)[colnames(sig_taxa_PD_long) == "diff_Case_statusPD"] <- "Significance"
colnames(sig_taxa_PD_long)[colnames(sig_taxa_PD_long) == "se_Case_statusPD"] <- "Standard_error"

sig_taxa_PD_long <- sig_taxa_PD_long[, -5]

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 7.5, height = 4,
          dir = "HMP2_Payami/Figures/for thesis")

ggplot(sig_taxa_PD_long, aes(x = LFC, y = taxon)) +
  geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue3") +
  geom_errorbarh(
    aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
    position = position_dodge(0.9),
    height = 0.25,  # Adjust the height of the error bars
    size = 0.5,    # Adjust the size of the error bars
    color = "gray73"  # Set the color of the error bars
  ) +
  labs(
    title = "Wallen PD-Associated Species",
    x = "Log Fold Change",
    y = "Species"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.position = "none", 
    legend.text = element_blank(), 
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 12, color = "black", face = "italic"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  )

ggsave("HMP2_Payami/Figures/PD significant species.png", dpi = 600, units = "in",
       height = 4.5, width = 8)

gg_stop_recording()



# plotting GENUS LEVEL DATA - WALLEN PD 
ancom_PD1 <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD genus ancombc2 output.rds")

res_prim_PD1 = ancom_PD1$res

sig_taxa_PD1 <- res_prim_PD1 %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))    # only 5 sigificant 

sig_taxa_PD <- sig_taxa_PD1[, c(1, 3, 7, 19, 23)]

sig_taxa_PD_long <- sig_taxa_PD %>%
  pivot_longer(
    cols = starts_with("lfc_"),  # Columns containing LFC values
    values_to = "LFC"           # New column name for LFC values
  )

colnames(sig_taxa_PD_long)[colnames(sig_taxa_PD_long) == "q_Case_statusPD"] <- "Adj P Value"
colnames(sig_taxa_PD_long)[colnames(sig_taxa_PD_long) == "diff_Case_statusPD"] <- "Significance"
colnames(sig_taxa_PD_long)[colnames(sig_taxa_PD_long) == "se_Case_statusPD"] <- "Standard_error"

sig_taxa_PD_long <- sig_taxa_PD_long[, -5]

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 6, height = 3,
          dir = "HMP2_Payami/Figures/for thesis")

ggplot(sig_taxa_PD_long, aes(x = LFC, y = taxon)) +
  geom_bar(stat = "identity", position = "dodge", fill = "skyblue") +
  geom_errorbarh(
    aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
    position = position_dodge(0.9),
    height = 0.25, 
    size = 0.5,    
    color = "gray73"  
  ) +
  labs(
    title = "Wallen PD-Associated Genera",
    x = "Log Fold Change",
    y = "Genus"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none", 
    legend.text = element_blank(), 
    axis.title.x = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 10, color = "black", face = "italic"),
    axis.text.x = element_text(size = 11, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  )

ggsave("HMP2_Payami/Figures/PD significant genera.png", dpi = 600, units = "in",
       height = 3, width = 6)

gg_stop_recording()


# --------------------------------------------------------------------------
# HMP2 IBD 
# plotting ancombc2 HMP2 data with CD UC combined into IBD 
# SPECIES LEVEL 
ancom <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD species ancombc2 output.rds")

res_prim = ancom$res

sig_taxa <- res_prim %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))    # 63 sig species  


sig_taxa_IBD <- sig_taxa[, c(1, 3, 7, 19, 23)]

sig_taxa_IBD_long <- sig_taxa_IBD %>%
  pivot_longer(
    cols = starts_with("lfc_"), 
    values_to = "LFC"          
  )

colnames(sig_taxa_IBD_long)[colnames(sig_taxa_IBD_long) == "q_diagnosis2IBD"] <- "Adj P Value"
colnames(sig_taxa_IBD_long)[colnames(sig_taxa_IBD_long) == "diff_diagnosis2IBD"] <- "Significance"
colnames(sig_taxa_IBD_long)[colnames(sig_taxa_IBD_long) == "se_diagnosis2IBD"] <- "Standard_error"

sig_taxa_IBD_long <- sig_taxa_IBD_long[,-5]

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 14, height = 13,
          dir = "HMP2_Payami/Figures/")

ggplot(sig_taxa_IBD_long, aes(x = LFC, y = taxon)) +
  geom_bar(stat = "identity", position = "dodge", fill = "darkseagreen4") +
  geom_errorbarh(
    aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
    position = position_dodge(0.9),
    height = 0.25,  
    size = 0.5,    
    color = "gray73"  
  ) +
  labs(
    title = "HMP2: IBD Associated Species",
    x = "Log Fold Change",
    y = "Species"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 21, face = "bold", hjust = 0.2),
    legend.position = "none", 
    legend.text = element_blank(), 
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 12, color = "black", face = "italic"),
    axis.text.x = element_text(size = 18, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  )

ggsave("HMP2_Payami/Figures/IBD significant species.png", dpi = 600, units = "in",
       height = 15, width = 12)

gg_stop_recording()



# GENUS LEVEL DATA - HMP2 IBD 
ancom_IBD <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD genus ancombc2 output.rds")

res_prim_IBD = ancom_IBD$res

sig_taxa_IBD <- res_prim_IBD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))    # 29 sig 

sig_taxa_IBD <- sig_taxa_IBD[, c(1, 3, 7, 19, 23)]

sig_taxa_IBD_long <- sig_taxa_IBD %>%
  pivot_longer(
    cols = starts_with("lfc_"),  
    values_to = "LFC"           
  )

colnames(sig_taxa_IBD_long)[colnames(sig_taxa_IBD_long) == "q_diagnosis2IBD"] <- "Adj P Value"
colnames(sig_taxa_IBD_long)[colnames(sig_taxa_IBD_long) == "diff_diagnosis2IBD"] <- "Significance"
colnames(sig_taxa_IBD_long)[colnames(sig_taxa_IBD_long) == "se_diagnosis2IBD"] <- "Standard_error"

sig_taxa_IBD_long <- sig_taxa_IBD_long[, -5]


gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 10, height = 8,
          dir = "HMP2_Payami/Figures/")

ggplot(sig_taxa_IBD_long, aes(x = LFC, y = taxon)) +
  geom_bar(stat = "identity", position = "dodge", fill = "darkseagreen") +
  geom_errorbarh(
    aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
    position = position_dodge(0.9),
    height = 0.25,
    size = 0.5,    
    color = "gray73"  
  ) +
  labs(
    title = "HMP2 IBD-Associated Genera",
    x = "Log Fold Change",
    y = "Genus"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.2),
    legend.position = "none", 
    legend.text = element_blank(), 
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 16, color = "black", face = "italic"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  )

ggsave("HMP2_Payami/Figures/IBD significant genera.png", dpi = 600, units = "in",
       height = 11, width = 10)

gg_stop_recording()
