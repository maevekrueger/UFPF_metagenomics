library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)

# Wallen PD data 
# first we need to plot aitchison distances on filtered dataframe because some samples
# are missing metadata and we can't run the ADONIS analysis with confounders if there
# are missing variables
PD <- readRDS("HMP2_Payami/Wallen_species_clr_counts.rds")
PD2 <- as.matrix(PD[, -(1:2)])

# Calculate Euclidean distance using the distance matrix "PD"
euclidean_dist_PD <- vegdist(PD2, method = "euclidean")

# Run PCoA on aitchison distances (euclidean distances from clr-transformed counts)
pcoa_PD <- pcoa(euclidean_dist_PD, correction = "none")

# Access the eigenvalues and eigenvectors
eigenvalues_PD <- pcoa_PD$values$values
eigenvectors_PD <- pcoa_PD$vectors

pcoa_PD <- data.frame(Sample = PD$Sample,
                      diagnosis = PD$diagnosis,
                      PCoA1 = eigenvectors_PD[, 1],
                      PCoA2 = eigenvectors_PD[, 2])


# confounders in metadata
PD_metadata <- readRDS("HMP2_Payami/Wallen PD Metadata.rds")
PD_metadata <- PD_metadata %>%
  rename(Age = Age_at_collection, Diagnosis = Case_status)
PD_metadata <- PD_metadata %>%
  mutate(Sex = ifelse(Sex == "M", "Male", ifelse(Sex == "F", "Female", Sex)))

PD_metadata$total_sequences <- scale(PD_metadata$total_sequences)
PD_metadata$Age <- scale(PD_metadata$Age)

# for some reason "SC0637" "SC0558" are not present so I am removing these from the 
# metadata so it doesn't cause problems when we combine dataframes
PD_metadata <- PD_metadata[!rownames(PD_metadata) %in% c("SC0637", "SC0558"), ]

# pull out predictor variables according to the demographics/metadata table stats 
selected_columns <- PD_metadata %>%
  select(Diagnosis, Sex, Age, total_sequences, Laxatives, Pain_med, Depression_anxiety_mood_med, 
         Birth_control_or_estrogen, Antihistamines, Probiotic, Sleep_aid)

selected_columns <- rownames_to_column(selected_columns, var = "Sample")

# Reorder selected_columns to match the order of Sample in pcoa_IBD
selected_columns <- selected_columns[match(pcoa_PD$Sample, selected_columns$Sample), ]


# combine 
combined <- cbind(pcoa_PD[, 1:2], selected_columns, pcoa_PD[, -c(1:2)])
combined <- combined[, -c(1, 2)]

combined$Sex <- factor(combined$Sex)
combined$Diagnosis <- factor(combined$Diagnosis)

convert_values <- function(x) {
  x <- ifelse(x %in% c("N", "N ", "never"), 0,
              ifelse(x %in% c("Y", "Y ", "yes"), 1, NA))
  return(x)
}
combined <- combined %>%
  mutate_at(vars(6:12), list(~convert_values(.)))


# remove any samples with missing metadata
combined <- combined[complete.cases(combined[, 6:12]), ]
pcoa_PD_df2 <- combined



# have to filter out these samples with missing metadata from the pcoa
PD <- readRDS("HMP2_Payami/Wallen_species_clr_counts.rds")

# filter out missing samples 
# Extract the Sample names from the combined dataframe
sample_names_to_keep <- combined$Sample

# match the sample names so PD is filtered properly 
PD_filtered <- PD[PD$Sample %in% sample_names_to_keep, ]

PD_filtered2 <- as.matrix(PD_filtered[, -(1:2)])

# Calculate Euclidean distance using the distance matrix "PD"
euclidean_dist_PD <- vegdist(PD_filtered2, method = "euclidean")

# Run PCoA on aitchison distances (euclidean distances from clr-transformed counts)
pcoa_PD <- pcoa(euclidean_dist_PD, correction = "none")

# Access the eigenvalues and eigenvectors
eigenvalues_PD <- pcoa_PD$values$values
eigenvectors_PD <- pcoa_PD$vectors

pcoa_PD <- data.frame(Sample = PD_filtered$Sample,
                      diagnosis = PD_filtered$diagnosis,
                      PCoA1 = eigenvectors_PD[, 1],
                      PCoA2 = eigenvectors_PD[, 2])



permanova_result <- adonis2(
  euclidean_dist_PD ~ Diagnosis + Sex + Age + total_sequences + Laxatives + 
    Pain_med + Depression_anxiety_mood_med + Birth_control_or_estrogen + 
    Antihistamines + Probiotic + Sleep_aid, data = pcoa_PD_df2
)

permanova_result  # significant: sex, age, total_sequences, laxatives, pain meds, depression/anxiety meds


# Save
write.xlsx(permanova_result, "HMP2_Payami/Wallen BETA DIVERSITY W confounders.xlsx")





# ---------------------------------------------------------------------------------
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
IBD2 <- as.matrix(IBD[, -(1:2)])

# Calculate Euclidean distance using the distance matrix "IBD"
euclidean_dist_IBD <- vegdist(IBD2, method = "euclidean")

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


# including confoundners 
# add in confounders from metadata file 
IBD_metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")
IBD_metadata <- column_to_rownames(IBD_metadata, var = "Sample")

# removing some columns I don't need
columns_to_remove <- c(14:19)
IBD_metadata <- IBD_metadata[, -columns_to_remove]

IBD_metadata <- IBD_metadata %>%
  rename(Age = consent_age, Diagnosis = diagnosis2, Sex = sex, Reads = reads_filtered,
         Immunosuppressants = Immunosuppressants..e.g..oral.corticosteroids., 
         Diagnosis_Age = Age.at.diagnosis, Site = site_name, 
         Colonoscopy_past_2_weeks = X2..In.the.past.2.weeks..have.you.undergone.a.colonoscopy.or.other.procedure,
         Diarrhea_past_2_weeks = X4..In.the.past.2.weeks..have.you.had.diarrhea., 
         Bowel_surgery = X6..Have.you.ever.had.bowel.surgery.)

# label probiotic column N or Y for probiotics taken within 7 days of collection
IBD_metadata <- IBD_metadata %>%
  mutate(Probiotic = case_when(
    grepl("^No", Probiotic) ~ "N",
    Probiotic == "" ~ NA_character_,
    TRUE ~ "Y"
  ))

order <- c("nonIBD", "IBD", "Overall")
IBD_metadata$Diagnosis <- factor(IBD_metadata$Diagnosis, levels = order)

IBD_metadata <- IBD_metadata %>%
  mutate(Age = as.numeric(Age)) 

IBD_metadata$Reads <- scale(IBD_metadata$Reads)
IBD_metadata$Age <- scale(IBD_metadata$Age)


# pull out predictor variables according to the demographics/metadata table stats 
selected_columns <- IBD_metadata %>%
  select(Sex, Age, Diagnosis, Reads, Site, Probiotic, Immunosuppressants, Antibiotics) %>%
  rownames_to_column(var = "Sample")


# Reorder selected_columns to match the order of Sample in pcoa_IBD
selected_columns <- selected_columns[match(pcoa_IBD$Sample, selected_columns$Sample), ]


# combine w/abundances 
combined <- cbind(pcoa_IBD[, 1:2], selected_columns, pcoa_IBD[, -c(1:2)])
combined <- combined[, -c(1, 2)]

combined$Sex <- factor(combined$Sex)
combined$Diagnosis <- factor(combined$Diagnosis)
combined$Site <- factor(combined$Site)

convert_values <- function(x) {
  x <- ifelse(x %in% c("N", "N ", "No"), 0,
              ifelse(x %in% c("Y", "Y ", "Yes"), 1, NA))
  return(x)
}
combined <- combined %>%
  mutate_at(vars(7:9), list(~convert_values(.)))

# remove any samples with missing metadata
combined <- combined[complete.cases(combined[, 6:9]), ]
pcoa_IBD_df2 <- combined



# have to filter out these samples with missing metadata from the pcoa
IBD <- readRDS("HMP2_Payami/HMP2 IBD Age Filtered species clr counts.rds")
# change samples with CD/UC to "IBD"
IBD <- IBD %>%
  mutate(diagnosis = ifelse(diagnosis == "CD", "IBD", diagnosis))
IBD <- IBD %>%
  mutate(diagnosis = ifelse(diagnosis == "UC", "IBD", diagnosis))

# changing Control label to non-IBD 
IBD <- IBD %>%
  mutate(diagnosis = ifelse(diagnosis == "Control", "non-IBD", diagnosis))

# filter out missing samples 
# Extract the Sample names from the combined dataframe
sample_names_to_keep <- combined$Sample

# match the sample names so PD is filtered properly 
IBD_filtered <- IBD[IBD$Sample %in% sample_names_to_keep, ]
IBD_filtered2 <- as.matrix(IBD_filtered[, -(1:2)])

# Calculate Euclidean distance using the distance matrix "IBD"
euclidean_dist_IBD <- vegdist(IBD_filtered2, method = "euclidean")

# Run PCoA on aitchison distances
pcoa_IBD <- pcoa(euclidean_dist_IBD, correction = "none")

# Access the eigenvalues and eigenvectors 
eigenvalues_IBD <- pcoa_IBD$values$values
eigenvectors_IBD <- pcoa_IBD$vectors

pcoa_IBD <- data.frame(Sample = IBD_filtered$Sample,
                       diagnosis = IBD_filtered$diagnosis,
                       PCoA1 = eigenvectors_IBD[, 1],
                       PCoA2 = eigenvectors_IBD[, 2])

pcoa_IBD$diagnosis <- factor(pcoa_IBD$diagnosis, levels = c("non-IBD", "IBD"))



# Run PERMANOVA (adonis2 from vegan package)
permanova_result2 <- adonis2(
  euclidean_dist_IBD ~ Sex + Age + Diagnosis + Reads + Site + Probiotic +
    Immunosuppressants + Antibiotics, data = pcoa_IBD_df2
)

permanova_result2  #everything is significant  0.001

# Save
write.xlsx(permanova_result2, "HMP2_Payami/HMP2 BETA DIVERSITY W confounders.xlsx")

