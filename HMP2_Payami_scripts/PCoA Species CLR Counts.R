library(vegan)
library(ggplot2)
library(ape)
library(pairwiseAdonis)
library(stats)
library(dplyr)

# PCoA on CLR transformed Counts
PD <- readRDS("HMP2_Payami/Wallen_species_clr_counts.rds")
PD <- as.matrix(PD[, -(1:2)])

# Calculate Euclidean distance using the distance matrix "PD"
euclidean_dist_PD <- vegdist(PD, method = "euclidean")

# Run PCoA on aitchison distances (euclidean distances from clr-transformed counts)
pcoa_PD <- pcoa(euclidean_dist_PD, correction = "none")

# Access the eigenvalues and eigenvectors
eigenvalues_PD <- pcoa_PD$values$values
eigenvectors_PD <- pcoa_PD$vectors

pcoa_PD <- data.frame(Sample = PD$Sample,
                      diagnosis = PD$diagnosis,
                      PCoA1 = eigenvectors_PD[, 1],
                      PCoA2 = eigenvectors_PD[, 2])

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 8, height = 6,
          dir = "HMP2_Payami/Figures")

# PD plot 
ggplot(pcoa_PD, aes(x = PCoA1, y = PCoA2, color = diagnosis)) +
  geom_point(size = 2.0) +
  stat_ellipse() +
  scale_color_manual(values = c("mediumorchid1", "royalblue")) + 
  ylab("PCoA 2") +
  ggtitle("Wallen PCoA Aitchison Distances") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
  labs(color = NULL) + # Remove the legend title
  theme(
    plot.title = element_text(size =18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black"))

ggsave("HMP2_Payami/Figures/Payami PD PCoA aitch diss.png", dpi = 600, units = "in",
       height = 6, width = 8)

gg_stop_recording()

# PERMANOVA
permanova1 <- adonis2(euclidean_dist_PD ~ diagnosis, data = pcoa_PD)
permanova1    # significant effect of diagnosis 0.001

# Run PERMDISP
# calculates multivariate dispersions based on Aitchison distances
# assesses how spread out the samples are within each group
permdisp_result2 <- betadisper(euclidean_dist_PD, pcoa_PD$diagnosis)
permdisp_result2

# assess pairwise differences in multivariate dispersion between groups
permutest <- permutest(permdisp_result2, pairwise = TRUE, permutations = 999)
permutest

plot(permdisp_result2, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")



# HMP2 IBD 
# PCoA on CLR transformed Counts
IBD <- readRDS("HMP2_Payami/HMP2 IBD Age Filtered species clr counts.rds")

# change samples with CD/UC to "IBD"
IBD <- IBD %>%
  mutate(diagnosis = ifelse(diagnosis == "CD", "IBD", diagnosis))
IBD <- IBD %>%
  mutate(diagnosis = ifelse(diagnosis == "UC", "IBD", diagnosis))

# changing Control label to non-IBD 
IBD <- IBD %>%
  mutate(diagnosis = ifelse(diagnosis == "Control", "non-IBD", diagnosis))
IBD <- as.matrix(IBD[, -(1:2)])

# Calculate Euclidean distance using the distance matrix "IBD"
euclidean_dist_IBD <- vegdist(IBD, method = "euclidean")

# Run PCoA on aitchison distances
pcoa_IBD <- pcoa(euclidean_dist_IBD, correction = "none")

# Access the eigenvalues and eigenvectors 
eigenvalues_IBD <- pcoa_IBD$values$values
eigenvectors_IBD <- pcoa_IBD$vectors

pcoa_IBD <- data.frame(Sample = IBD$Sample,
                       diagnosis = IBD$diagnosis,
                       PCoA1 = eigenvectors_IBD[, 1],
                       PCoA2 = eigenvectors_IBD[, 2])

pcoa_IBD$diagnosis <- factor(pcoa_IBD$diagnosis, levels = c("non-IBD", "IBD"))

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 11, height = 5,
          dir = "HMP2_Payami/Figures")

# IBD plot 
ggplot(pcoa_IBD, aes(x = PCoA1, y = PCoA2, color = diagnosis)) +
  geom_point(size = 3.0) +
  stat_ellipse() +
  scale_color_manual(values = c("orange", "green3")) + 
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  ggtitle("HMP2 PCoA Aitchison Distances") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
  labs(color = NULL) +  # Remove the legend title
  theme(
    plot.title = element_text(size =21, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5)),
    strip.background = element_rect(fill = "black", color = "black"))

ggsave("HMP2_Payami/Figures/HMP2 IBD PCoA aitch diss.png", dpi = 600, units = "in",
       height = 6, width = 8)

gg_stop_recording()

# PERMANOVA
permanova_result2 <- adonis2(euclidean_dist_IBD ~ diagnosis, data = pcoa_IBD)
permanova_result2    # significant effect of diagnosis 0.001

# Run PERMDISP
# calculates multivariate dispersions based on Aitchison distances
# assesses how spread out the samples are within each group
permdisp_result2 <- betadisper(euclidean_dist_IBD, pcoa_IBD$diagnosis)
permdisp_result2

# assess pairwise differences in multivariate dispersion between groups
permutest <- permutest(permdisp_result2, pairwise = TRUE, permutations = 999)
permutest

plot(permdisp_result2, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
