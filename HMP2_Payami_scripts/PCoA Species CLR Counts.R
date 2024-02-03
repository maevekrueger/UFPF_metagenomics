library(vegan)
library(ggplot2)
library(ape)
library(pairwiseAdonis)
library(stats)
library(dplyr)

# PCoA on CLR transformed Counts
PD <- readRDS("HMP2_Payami/Payami PD species clr counts.rds")

# convert data to matrix, ignoring sample and diagnosis columns
PD1 <- as.matrix(PD[, -(1:2)])

# Calculate Euclidean distance using the distance matrix 
euclidean_dist_PD <- vegdist(PD1, method = "euclidean")

# Run PCoA on the euclidean distance matrix
pcoa_PD <- pcoa(euclidean_dist_PD, correction = "none")

# Access the eigenvalues and eigenvectors from the PCoA results
eigenvalues_PD <- pcoa_PD$values$values
eigenvectors_PD <- pcoa_PD$vectors

# Combine cohorts and produce PCoA data frames
pcoa_PD <- data.frame(Sample = PD$Sample,
                      diagnosis = PD$diagnosis,
                      PCoA1 = eigenvectors_PD[, 1],
                      PCoA2 = eigenvectors_PD[, 2])

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 8, height = 6,
          dir = "HMP2_Payami/Figures/for thesis")

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



# combined HMP2 IBD 
# PCoA on CLR transformed Counts
IBD <- readRDS("HMP2_Payami/HMP2 IBD Age Filtered species clr counts.rds")

# change samples with CD/UC to "IBD"
IBD <- IBD %>%
  mutate(diagnosis = ifelse(diagnosis == "CD", "IBD", diagnosis))
IBD <- IBD %>%
  mutate(diagnosis = ifelse(diagnosis == "UC", "IBD", diagnosis))

# changing Control label to nonIBD 
IBD <- IBD %>%
  mutate(diagnosis = ifelse(diagnosis == "Control", "nonIBD", diagnosis))

# convert data to matrix, ignoring sample and diagnosis columns
IBD1 <- as.matrix(IBD[, -(1:2)])

# Calculate Euclidean distance using the distance matrix 
euclidean_dist_IBD <- vegdist(IBD1, method = "euclidean")

# Run PCoA on the euclidean distance matrix
pcoa_IBD <- pcoa(euclidean_dist_IBD, correction = "none")

# Access the eigenvalues and eigenvectors from the PCoA results
eigenvalues_IBD <- pcoa_IBD$values$values
eigenvectors_IBD <- pcoa_IBD$vectors

# Combine cohorts and produce PCoA data frames
pcoa_IBD <- data.frame(Sample = IBD$Sample,
                       diagnosis = IBD$diagnosis,
                       PCoA1 = eigenvectors_IBD[, 1],
                       PCoA2 = eigenvectors_IBD[, 2])

pcoa_IBD$diagnosis <- factor(pcoa_IBD$diagnosis, levels = c("nonIBD", "IBD"))

gg_record(device = "png",
          dpi = 600,
          units = "in",
          width = 11, height = 5,
          dir = "HMP2_Payami/Figures/")

# IBD plot 
ggplot(pcoa_IBD, aes(x = PCoA1, y = PCoA2, color = diagnosis)) +
  geom_point(size = 3.0) +
  stat_ellipse() +
  scale_color_manual(values = c("orange", "green3")) + # Change the color scheme as desired
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

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 IBD PCoA aitch diss.png", dpi = 600, units = "in",
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
