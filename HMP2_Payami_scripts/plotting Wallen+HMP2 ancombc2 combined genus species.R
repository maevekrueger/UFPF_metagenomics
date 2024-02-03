depleted_IBD <- readRDS("HMP2_Payami/ANCOMBC2/depleted genus species HMP2 IBD.rds") 
enriched_IBD <-  readRDS("HMP2_Payami/ANCOMBC2/enriched genus species HMP2 IBD.rds")
depleted_PD <- readRDS("HMP2_Payami/ANCOMBC2/depleted genus species Wallen PD.rds") 
enriched_PD<- readRDS("HMP2_Payami/ANCOMBC2/enriched genus species Wallen PD.rds")

IBD <- rbind(depleted_IBD, enriched_IBD)    # 92
PD <- rbind(depleted_PD, enriched_PD)       # 16

# Wallen PD 
sig_taxa_PD <- PD[, c(1, 3, 7, 19, 23)]

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
          width = 10, height = 6,
          dir = "HMP2_Payami/Figures/")

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
    title = "Wallen: PD Associated Taxa",
    x = "Log Fold Change",
    y = "Taxa"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 21, face = "bold", hjust = 0.2),
    legend.position = "none", 
    legend.text = element_blank(), 
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black", face = "italic"),
    axis.text.x = element_text(size = 18, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  )

ggsave("HMP2_Payami/Figures/Wallen PD genus species ancombc2.png", dpi = 600, units = "in",
       height = 11, width = 10)

gg_stop_recording()


# plotting HMP2 IBD 
# combined into IBD 
sig_taxa_IBD <- IBD[, c(1, 3, 7, 19, 23)]

sig_taxa_IBD_long <- sig_taxa_IBD %>%
  pivot_longer(
    cols = starts_with("lfc_"),  # Columns containing LFC values
    values_to = "LFC"           # New column name for LFC values
  )

colnames(sig_taxa_IBD_long)[colnames(sig_taxa_IBD_long) == "q_diagnosis2IBD"] <- "Adj P Value"
colnames(sig_taxa_IBD_long)[colnames(sig_taxa_IBD_long) == "diff_diagnosis2IBD"] <- "Significance"
colnames(sig_taxa_IBD_long)[colnames(sig_taxa_IBD_long) == "se_diagnosis2IBD"] <- "Standard_error"

sig_taxa_IBD_long <- sig_taxa_IBD_long[,-5]

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 15, height = 18,
          dir = "HMP2_Payami/Figures/")

ggplot(sig_taxa_IBD_long, aes(x = LFC, y = taxon)) +
  geom_bar(stat = "identity", position = "dodge", fill = "darkseagreen4") +
  geom_errorbarh(
    aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
    position = position_dodge(0.9),
    height = 0.25,  # Adjust the height of the error bars
    size = 0.5,    # Adjust the size of the error bars
    color = "gray73"  # Set the color of the error bars
  ) +
  labs(
    title = "HMP2: IBD Associated Taxa",
    x = "Log Fold Change",
    y = "Taxa"
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

ggsave("HMP2_Payami/Figures/HMP2 IBD genus species ancombc2.png", dpi = 600, units = "in",
       height = 15, width = 12)

gg_stop_recording()

