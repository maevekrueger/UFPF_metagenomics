# creating phyloseq objects for ANCOMBC2
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

# (1) CREATING PHYLOSEQ FROM RAW COUNTS WITH UNCLASSIFIED 
# MAKE TAXONOMY TABLE 
# read in metaphlan data 
s_abund2 <- readRDS("UFPF/Metaphlan output/Counts/Counts w Unclassified.rds")
s_abund2 <- column_to_rownames(s_abund2, var = "Sample")
s_abund2 <- s_abund2[, -c(1:2)]

s_abund2 <- t(s_abund2)
rownames(s_abund2) <- gsub("\\|", ".", rownames(s_abund2))

# Function to extract the last level of classification from the row names
extract_last_level <- function(row_names) {
  last_levels <- sub('.*\\|', '', row_names)
  last_levels
}

# Function to extract the classification level from the row names
extract_classification_level <- function(row_names) {
  levels <- sub('k__|p__|c__|o__|f__|g__|s__', '', row_names)
  levels
}

# Create a new data frame to store the transformed data
transformed_data <- data.frame(
  Taxonomy = extract_last_level(rownames(s_abund2)),
  Classification_Level = extract_classification_level(rownames(s_abund2)),
  s_abund2
)

# isolate species names from row names - 1067 species 
transformed_data <- transformed_data[grep("\\.s__[^.]*$", rownames(transformed_data)), ]

transformed_data_levels <- sub("^.*\\.s__(.+)$", "\\1", rownames(transformed_data))
rownames(transformed_data) <- transformed_data_levels


# filtering out species that aren't present in at least 10% of samples 
# Calculate the number of samples
total_samples <- ncol(transformed_data) - 2                # 96 total samples 
# Calculate the threshold count level (10% of total samples)
threshold_count <- 0.10 * total_samples                    # threshold is 9.6 samples     
# Filter rows based on the threshold count
filtered_data <- transformed_data[rowSums(transformed_data[, -c(1, 2)] != 0) >= threshold_count, ]
transformed_data <- filtered_data                          # filtered from 1167 to 233 species 

# move row names to a new column 
transformed_data <- transformed_data %>%
  rownames_to_column(var = "Species")

# remove unnecessary column 
transformed_data <- transformed_data[, -3]

# Separating taxonomic levels into separate columns 
separated_data <- transformed_data %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "\\.")

# clean up 
separated_data <- apply(separated_data, 2, function(col) gsub(".*__", "", col))
separated_data <- as.data.frame(separated_data)
separated_data <- separated_data[, c(1:7)]

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(separated_data)))
# Assign labels to row names
rownames(separated_data) <- labels

#----------------------------------------------------------------------------------------------
# MAKE OTU TABLE 
# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data)))
# Assign labels to row names
rownames(transformed_data) <- labels

# remove columns we don't need in the OTU table 
transformed_data <- transformed_data[, -(1:2)]

# -------------------------------------------------------------------------------------------
# creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
Metadata <- readRDS("UFPF/Metadata.rds")

convert_values <- function(x) {
  x <- ifelse(x %in% c("N", "N ", "never"), 0,
              ifelse(x %in% c("Y", "Y ", "yes"), 1, NA))
  return(x)
}
Metadata <- Metadata %>%
  mutate_at(vars(13:31), list(~convert_values(.)))

# Scaling Reads using Scale function (to be used as fixed effect later)
Metadata$Reads <- scale(Metadata$Reads)

# transform OTU and taxonomy table into matrices
separated_data <- as.matrix(separated_data)
transformed_data <- as.matrix(transformed_data)

# create phyloseq object
OTU = otu_table(transformed_data, taxa_are_rows = TRUE)
TAX = tax_table(separated_data)
samples = sample_data(Metadata)

phyloseq_object <- phyloseq(OTU, TAX, samples)

# save object 
saveRDS(phyloseq_object, "UFPF/Phyloseq Objects/phyloseq object species effect of meds.rds")


# running in hipergator 
phyloseq_object <- readRDS("UFPF/Phyloseq Objects/phyloseq object species effect of meds.rds")

# running ancombc 
ancom <- ancombc2(phyloseq_object, 
                  fix_formula = "Laxatives + Indigestion.meds + Anti.inflammatories..non.NSAID. + Anti.TNF + Cholesterol.meds + Antihistamines + Diabetes.meds + Depression.anxiety.meds + NSAIDs",
                  tax_level = "Species",
                  p_adj_method = "BH",
                  group = "Diagnosis2",
                  lib_cut=0,
                  struc_zero=FALSE,
                  neg_lb=FALSE,
                  alpha=0.05,
                  global = FALSE, 
                  pairwise = FALSE)

saveRDS(ancom, "UFPF/Effect of Meds/ancombc2 species effect of meds.rds")

ancom <- readRDS("UFPF/Effect of Meds/ancombc2 species effect of meds.rds")
res_prim = ancom$res

sig_taxa <- res_prim %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_"))))    # 0 significant species 

write.csv(sig_taxa, "UFPF/Effect of Meds/medication effects on species.csv", row.names = FALSE)




# redoing at the genus level 
# (1) CREATING PHYLOSEQ FROM RAW COUNTS WITH UNCLASSIFIED 
# MAKE TAXONOMY TABLE 
# read in metaphlan data 
s_abund2 <- readRDS("UFPF/Metaphlan output/Counts/Counts w Unclassified.rds")
s_abund2 <- column_to_rownames(s_abund2, var = "Sample")
s_abund2 <- s_abund2[, -c(1:2)]

s_abund2 <- t(s_abund2)
rownames(s_abund2) <- gsub("\\|", ".", rownames(s_abund2))

# Function to extract the last level of classification from the row names
extract_last_level <- function(row_names) {
  last_levels <- sub('.*\\|', '', row_names)
  last_levels
}

# Function to extract the classification level from the row names
extract_classification_level <- function(row_names) {
  levels <- sub('k__|p__|c__|o__|f__|g__|', '', row_names)
  levels
}

# Create a new data frame to store the transformed data
transformed_data <- data.frame(
  Taxonomy = extract_last_level(rownames(s_abund2)),
  Classification_Level = extract_classification_level(rownames(s_abund2)),
  s_abund2
)

# isolate species names from row names - 595 genera 
transformed_data <- transformed_data[grep("\\.g__[^.]*$", rownames(transformed_data)), ]

transformed_data_levels <- sub("^.*\\.g__(.+)$", "\\1", rownames(transformed_data))
rownames(transformed_data) <- transformed_data_levels


# filtering out species that aren't present in at least 10% of samples 
# Calculate the number of samples
total_samples <- ncol(transformed_data) - 2                # 96 total samples 
# Calculate the threshold count level (10% of total samples)
threshold_count <- 0.10 * total_samples                    # threshold is 9.6 samples     
# Filter rows based on the threshold count
filtered_data <- transformed_data[rowSums(transformed_data[, -c(1, 2)] != 0) >= threshold_count, ]
transformed_data <- filtered_data                          # filtered from 650 to 140 genera

# move row names to a new column 
transformed_data <- transformed_data %>%
  rownames_to_column(var = "Genus")

# remove unnecessary column 
transformed_data <- transformed_data[, -3]

# Separating taxonomic levels into separate columns 
separated_data <- transformed_data %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = "\\.")

# clean up 
separated_data <- apply(separated_data, 2, function(col) gsub(".*__", "", col))
separated_data <- as.data.frame(separated_data)
separated_data <- separated_data[, c(1:7)]

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(separated_data)))
# Assign labels to row names
rownames(separated_data) <- labels

#----------------------------------------------------------------------------------------------
# MAKE OTU TABLE 

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data)))
# Assign labels to row names
rownames(transformed_data) <- labels

# remove columns we don't need in the OTU table 
transformed_data <- transformed_data[, -(1:2)]

# -------------------------------------------------------------------------------------------
# creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
Metadata <- readRDS("UFPF/Metadata.rds")

convert_values <- function(x) {
  x <- ifelse(x %in% c("N", "N ", "never"), 0,
              ifelse(x %in% c("Y", "Y ", "yes"), 1, NA))
  return(x)
}

Metadata <- Metadata %>%
  mutate_at(vars(13:31), list(~convert_values(.)))

# Scaling Reads using Scale function (to be used as fixed effect later)
Metadata$Reads <- scale(Metadata$Reads)

# transform OTU and taxonomy table into matrices
separated_data <- as.matrix(separated_data)
transformed_data <- as.matrix(transformed_data)

# create phyloseq object
OTU = otu_table(transformed_data, taxa_are_rows = TRUE)
TAX = tax_table(separated_data)
samples = sample_data(Metadata)

phyloseq_object_g <- phyloseq(OTU, TAX, samples)

# save object 
saveRDS(phyloseq_object_g, "UFPF/Phyloseq Objects/phyloseq object genus effect of meds.rds")


phyloseq_object_g <- readRDS("UFPF/Phyloseq Objects/phyloseq object genus effect of meds.rds")

# running ancombc 
ancom <- ancombc2(phyloseq_object_g, 
                  fix_formula = "Laxatives + Indigestion.meds + Anti.inflammatories..non.NSAID. + Anti.TNF + Cholesterol.meds + Antihistamines + Diabetes.meds + Depression.anxiety.meds + NSAIDs",
                  tax_level = "Genus",
                  p_adj_method = "BH",
                  group = "Diagnosis2",
                  lib_cut=0,
                  struc_zero=FALSE,
                  neg_lb=FALSE,
                  alpha=0.05,
                  global = FALSE, 
                  pairwise = FALSE)

saveRDS(ancom, "UFPF/Effect of Meds/ancombc2 genus effect of meds.rds")

ancom <- readRDS("UFPF/Effect of Meds/ancombc2 genus effect of meds.rds")

res_prim = ancom$res

sig_taxa <- res_prim %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_"))))    #  1 significant - Faecalimonas 

write.csv(sig_taxa, "UFPF/Effect of Meds/medication effects on genus.csv", row.names = FALSE)
