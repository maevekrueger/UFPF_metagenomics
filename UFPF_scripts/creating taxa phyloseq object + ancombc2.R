# creating phyloseq objects for ANCOMBC2
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

# (1) CREATING PHYLOSEQ FROM RAW COUNTS WITH UNCLASSIFIED ESTIMATION 
# MAKE TAXONOMY TABLE ---  read in metaphlan data 
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

transformed_data <- transformed_data[, -3]  # remove column we don't need

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
# (2) MAKE OTU TABLE 
# we obviously don't have "OTUs" - this is just based on the tutorial, and OTU is just used here as the stand-in for the species name 

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data)))
# Assign labels to row names
rownames(transformed_data) <- labels

# remove columns we don't need in the OTU table 
transformed_data <- transformed_data[, -(1:2)]

# -------------------------------------------------------------------------------------------
# (3) creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
Metadata <- readRDS("UFPF/Metadata.rds")

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
saveRDS(phyloseq_object, "UFPF/phyloseq object species.rds")


# running in hipergator 
phyloseq_object <- readRDS("UFPF/phyloseq object species.rds")

# running ancombc 
ancom <- ancombc2(phyloseq_object, 
                  fix_formula = "Diagnosis2 + Reads",
                  tax_level = "Species",
                  p_adj_method = "BH",
                  group = "Diagnosis2",
                  lib_cut=0,
                  struc_zero=FALSE,
                  neg_lb=FALSE,
                  alpha=0.05,
                  global = FALSE, 
                  pairwise = TRUE)

saveRDS(ancom, "UFPF/ANCOMBC2/UFPF ancombc2 species no sex age.rds")

res_prim = ancom$res
res_pair = ancom$res_pair     

sig_taxa <- res_pair %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_"))))   # 1 significant: Faecalimonas_umbilicata

write.csv(sig_taxa, file = "UFPF/ANCOMBC2/ANCOMBC2 tables/No sex or age/enriched species IBD.csv", row.names = FALSE)



# GENUS-LEVEL
# a separate phyloseq object was created because some genera were not detected at the species-level
# so they would have been filtered out in the previous species-level phyloseq object creation 

# (1) CREATING PHYLOSEQ FROM RAW COUNTS WITH UNCLASSIFIED ESTIMATION
# MAKE TAXONOMY TABLE -- read in metaphlan data 
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
# (2) MAKE OTU TABLE 

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data)))
# Assign labels to row names
rownames(transformed_data) <- labels

# remove columns we don't need in the OTU table 
transformed_data <- transformed_data[, -(1:2)]

# -------------------------------------------------------------------------------------------
# (3) creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
Metadata <- readRDS("UFPF/Metadata.rds")

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
saveRDS(phyloseq_object_g, "UFPF/phyloseq object genus.rds")


phyloseq_object_g <- readRDS("UFPF/phyloseq object genus.rds")

# running ancombc 
ancom <- ancombc2(phyloseq_object_g, 
                  fix_formula = "Diagnosis2 + Reads",
                  tax_level = "Genus",
                  p_adj_method = "BH",
                  group = "Diagnosis2",
                  lib_cut=0,
                  struc_zero=FALSE,
                  neg_lb=FALSE,
                  alpha=0.05,
                  global = FALSE, 
                  pairwise = TRUE)

saveRDS(ancom, "UFPF/ANCOMBC2/UFPF ancombc2 genus.rds")

res_prim = ancom$res
res_pair = ancom$res_pair      

sig_taxa <- res_pair %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_"))))   # 2 significant: Klebsiella + Faecalimonas

write.csv(sig_taxa, file = "UFPF/ANCOMBC2/ANCOMBC2 tables/enriched genera IBD.csv", row.names = FALSE)
