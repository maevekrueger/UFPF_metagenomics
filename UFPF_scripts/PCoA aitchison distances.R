# PCoA on Aitchison distances 
library(compositions)
library(vegan)
library(dplyr)
library(ape)
library(ggplot2)
library(pairwiseAdonis)

# read in raw counts 
counts <- readRDS("UFPF/Metaphlan output/Counts/species raw counts.rds")

# add a small pseudocount to avoid any zeros in the df which clr doesn't like 
counts <- counts + 1e-10

# perform clr transformation
# make sure this is using the compositions package to perform clr 
clr_species <- compositions::clr(counts)

# subset cohorts 
PD_s <- clr_species[grep("1N", rownames(clr_species)), ]
Control_s <- clr_species[grep("3N|3G", rownames(clr_species)), ]
IBD_s <- clr_species[grep("2G", rownames(clr_species)), ]

# Combine the data into one data frame for each cohort
all_cohorts_s <- rbind(PD_s, Control_s, IBD_s)

# convert to data.frame 
all_cohorts_s <- as.data.frame(all_cohorts_s)

# these lines of code create a new data frame that has columns for sample ID and cohort 
all_cohorts_s <- cbind(
  data.frame(Sample_ID = rownames(all_cohorts_s),
             cohort = ifelse(grepl("1N", rownames(all_cohorts_s)), "PD",
                             ifelse(grepl("2G", rownames(all_cohorts_s)), "IBD",
                                    ifelse(grepl("3N|3G", rownames(all_cohorts_s)), "Control", NA)))),
  all_cohorts_s
)


# convert data to matrix, ignoring sample and cohort columns
all_cohorts_s2 <- all_cohorts_s[, -(1:2)]
matrix_s <- as.matrix(all_cohorts_s2)

# Calculate Euclidean distance using the distance matrix 
euclidean_dist_s <- vegdist(matrix_s, method = "euclidean")

# Run PCoA on the euclidean distance matrix
pcoa_s <- pcoa(euclidean_dist_s, correction = "none")

# Access the eigenvalues and eigenvectors from the PCoA results
eigenvalues_s <- pcoa_s$values$values
eigenvectors_s <- pcoa_s$vectors

# Combine cohorts and produce PCoA data frames
pcoa_s_df <- data.frame(Sample_ID = all_cohorts_s$Sample_ID,
                        cohort = all_cohorts_s$cohort,
                        PCoA1 = eigenvectors_s[, 1],
                        PCoA2 = eigenvectors_s[, 2])


ggplot(pcoa_s_df, aes(x = PCoA1, y = PCoA2, color = cohort)) +
  geom_point(size = 3.0) +
  stat_ellipse() +
  scale_color_manual(values = c("mediumorchid1", "limegreen", "royalblue")) + # Change the color scheme as desired
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  ggtitle("PCoA Aitchison Distances") +
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
    strip.background = element_rect(fill = "black", color = "black"))

ggsave("UFPF/Figures/PCoA species clr transformed counts.png", dpi = 600, units = "in",
       height = 6, width = 8)


# Run PERMANOVA
# Assuming 'cohort' is the grouping variable
permanova_result3 <- adonis2(euclidean_dist_s ~ cohort, data = pcoa_s_df)
permanova_result3     # 0.001*

# different way of calculating 
pairwise_result3 <- pairwise.adonis(
  euclidean_dist_s,
  pcoa_s_df$cohort,
  #sim.method = "euclidean", # don't think I need this b/c using a distance matrix 
  p.adjust.m = "bonferroni",
  reduce = NULL,
  perm = 999
)
pairwise_result3
# PD-Control  0.033*
# PD-IBD      0.003*
# IBD-Control 0.003*


