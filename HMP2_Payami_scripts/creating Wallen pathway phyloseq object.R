library(readxl)
library(tidyverse)
library(phyloseq)

# creating humann pathway abundance phyloseq object 
# Payami = Wallen PD dataset
payami_pathways <- read_excel("HMP2_Payami/PAYAMI DATA- Source_Data_24Oct2022.xlsx", sheet = 6)

payami_pathways <- column_to_rownames(payami_pathways, var = "Pathway")

# Remove '[' and everything after it including "
rownames(payami_pathways) <- sub("\\[.*", "", rownames(payami_pathways))
rownames(payami_pathways) <- gsub('"', '', rownames(payami_pathways))

# filtering out pathways that aren't present in at least 25% of samples 
# Calculate the number of samples
total_samples <- ncol(payami_pathways) 
# Calculate the threshold count level (25% of total samples)
threshold_count <- 0.25 * total_samples
# Filter rows based on the threshold count
filtered_data <- payami_pathways[rowSums(payami_pathways != 0) >= threshold_count, ]
transformed_data <- filtered_data
transformed_data <- as.data.frame(transformed_data)

# move row names to a new column 
transformed_data <- transformed_data %>%
  rownames_to_column(var = "Genus")   

# Create a new separated data frame with the specified column names
separated_data <- data.frame(
  Kingdom = character(0),
  Phylum = character(0),
  Class = character(0),
  Order = character(0),
  Family = character(0),
  Genus = character(0)
)

separated_data <- separated_data[1:nrow(transformed_data), ]
separated_data$Genus <- transformed_data$Genus
separated_data[, -which(names(separated_data) == "Genus")] <- "NA"

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(separated_data)))
# Assign labels to row names
rownames(separated_data) <- labels

#----------------------------------------------------------------------------------------------
# MAKE OTU TABLE 
# we clearly don't have OTUs - this is just to keep the terminology consistent with the
#phyloseq tutorial we used. OTU here is just a stand-in for the pathway name 
# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data)))
# Assign labels to row names
rownames(transformed_data) <- labels

# remove column we don't need in the OTU table 
transformed_data <- transformed_data[, -1]

# -------------------------------------------------------------------------------------------
# creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
Metadata <- readRDS("HMP2_Payami/Wallen PD Metadata.rds")

# remove unnecessary column
Metadata <- Metadata[, -which(names(Metadata) == "Sample.1")]

# OPTIONAL scaling Reads using Scale function (to be used as fixed effect later)
Metadata$total_sequences <- scale(Metadata$total_sequences)

# transform OTU and taxonomy table into matrices
separated_data <- as.matrix(separated_data)
transformed_data <- as.matrix(transformed_data)

# create phyloseq object
OTU = otu_table(transformed_data, taxa_are_rows = TRUE)
TAX = tax_table(separated_data)
samples = sample_data(Metadata)

humann_phyloseq_object2 <- phyloseq(OTU, TAX, samples)

# save object 
saveRDS(humann_phyloseq_object2, "HMP2_Payami/Phyloseq Objects/Wallen PD pathways phyloseq object.rds")

