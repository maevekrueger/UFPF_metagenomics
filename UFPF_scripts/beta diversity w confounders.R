library(compositions)
library(vegan)
library(ggplot2)
library(ape)
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

all_cohorts_s <- rbind(PD_s, Control_s, IBD_s)
all_cohorts_s <- as.data.frame(all_cohorts_s)

# create a new data frame with columns for sample ID and cohort 
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


# BETA DIVERSITY ANALYSIS
# Calculate Euclidean distance using the distance matrix (again these are technically now considered Aitchison distances)
euclidean_dist_s <- vegdist(matrix_s, method = "euclidean")

# Run PCoA on the aitchison distance matrix
pcoa_s <- pcoa(euclidean_dist_s, correction = "none")

# Access the eigenvalues and eigenvectors from the PCoA results
eigenvalues_s <- pcoa_s$values$values
eigenvectors_s <- pcoa_s$vectors

# Combine cohorts and produce PCoA data frames
pcoa_s_df <- data.frame(Sample_ID = all_cohorts_s$Sample_ID,
                        cohort = all_cohorts_s$cohort,
                        PCoA1 = eigenvectors_s[, 1],
                        PCoA2 = eigenvectors_s[, 2])


# -------------------------------------------------------------------------------
# including confoundners 
# add in confounders from metadata file 
Metadata <- readRDS("UFPF/Metadata.rds")
colnames(Metadata)[c(2,6,12:53)] <- gsub("\\.", "_", colnames(Metadata)[c(2,6,12:53)])

# Scaling 
Metadata$Reads <- scale(Metadata$Reads)
Metadata$Age <- scale(Metadata$Age)

# pull out predictor variables according to the demographics/metadata table stats 
selected_columns <- Metadata %>%
  select(Sex, Age, Diagnosis2, Reads, Indigestion_meds, Anti_TNF, Anti_inflammatories__non_NSAID_, Depression_anxiety_meds, Iron_specific_supplement)
selected_columns <- selected_columns %>%
  rename(Anti_inflammatories = Anti_inflammatories__non_NSAID_)
selected_columns <- selected_columns %>%
  rename(Diagnosis = Diagnosis2)

selected_columns$Sample_ID <- rownames(selected_columns)
selected_columns <- selected_columns %>%
  select(Sample_ID, everything())

# Reorder selected_columns to match the order of Sample in pcoa_IBD
reordered <- selected_columns[match(pcoa_s_df$Sample, selected_columns$Sample), ]


# combine w/abundances 
combined <- cbind(pcoa_s_df[, 1:2], reordered, pcoa_s_df[, -c(1:2)])
combined <- combined[, -c(1, 2)]

combined$Sex <- factor(combined$Sex)
combined$Diagnosis <- factor(combined$Diagnosis)

convert_values <- function(x) {
  x <- ifelse(x %in% c("N", "N ", "never"), 0,
              ifelse(x %in% c("Y", "Y ", "yes"), 1, NA))
  return(x)
}
combined <- combined %>%
  mutate_at(vars(6:10), list(~convert_values(.)))

pcoa_s_df2 <- combined


# Run PERMANOVA (adonis2 from vegan package)
permanova_result3 <- adonis2(
  euclidean_dist_s ~ Sex + Age + Diagnosis + Reads + Indigestion_meds + 
    Anti_TNF + Anti_inflammatories + Depression_anxiety_meds + 
    Iron_specific_supplement, data = pcoa_s_df2
)

                             
permanova_result3     # significant: sex, age, diagnosis, read count


# Save
write.xlsx(permanova_result3, "UFPF/UFPF BETA DIVERSITY W confounders.xlsx")


