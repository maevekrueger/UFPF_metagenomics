# HUMANN FUNCTIONAL ANALYSIS 

# create Phyloseq object 
#  NEED TO DO THIS WITH CONSOLIDATED DATA ****************
# filtering out genes that aren't present in at least 25% of samples 
load("HMP2_Payami/HMP2 Consolidated and Age Filtered KO Groups.RData")
IBD_KO <- t(df)

# Calculate the number of samples
total_samples <- ncol(IBD_KO)             # 332 
# Calculate the threshold count level (25% of total samples)
threshold_count <- 0.25 * total_samples                        # 83
# Filter rows based on the threshold count
filtered_data <- IBD_KO[rowSums(IBD_KO != 0) >= threshold_count, ]
transformed_data <- filtered_data
transformed_data <- as.data.frame(transformed_data)            # 7462 filtered to 2767 KO groups 

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
# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data)))
# Assign labels to row names
rownames(transformed_data) <- labels

# remove column we don't need in the OTU table 
transformed_data <- transformed_data[, -1]


# -------------------------------------------------------------------------------------------
# creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
Metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")

# remove samples from Metadata that are missing in KO group file 
missing_samples <- setdiff(Metadata$Sample, colnames(IBD_KO))
missing_samples  # ^ "CSM67UDR_TR" "CSM79HLA_TR" "CSM79HPA_TR" "MSM79HF9_TR" "MSM9VZNH_TR"
# Remove the missing samples from Metadata
Metadata <- Metadata[!(Metadata$Sample %in% missing_samples), ]

# OPTIONAL scaling Reads using Scale function (to be used as fixed effect later)
Metadata$reads_filtered <- scale(Metadata$reads_filtered)

rownames(Metadata) <- Metadata$Sample
Metadata <- Metadata[, -1]

# transform OTU and taxonomy table into matrices
separated_data <- as.matrix(separated_data)
transformed_data <- as.matrix(transformed_data)

# create phyloseq object
OTU = otu_table(transformed_data, taxa_are_rows = TRUE)
TAX = tax_table(separated_data)
samples = sample_data(Metadata)

humann_phyloseq_object <- phyloseq(OTU, TAX, samples)

# save object 
saveRDS(humann_phyloseq_object, "HMP2_Payami/Phyloseq Objects/HMP2 KO Groups phyloseq object.rds")





# create Phyloseq object - Pathways 
save(paths, file = "Hmp2_Payami/HMP2 Consolidated and Age Filtered Pathways.RData")
IBD_paths <- t(paths)

# filtering out genes that aren't present in at least 25% of samples 
# Calculate the number of samples
total_samples <- ncol(IBD_paths)                              # 332
# Calculate the threshold count level (25% of total samples)
threshold_count <- 0.25 * total_samples                       # 83
# Filter rows based on the threshold count
filtered_data <- IBD_paths[rowSums(IBD_paths != 0) >= threshold_count, ]
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
# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data)))
# Assign labels to row names
rownames(transformed_data) <- labels

# remove column we don't need in the OTU table 
transformed_data <- transformed_data[, -1]

# -------------------------------------------------------------------------------------------
# creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
Metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")

# remove samples from Metadata that are missing in original file 
missing_samples <- setdiff(Metadata$Sample, colnames(IBD_paths))
missing_samples  # ^ "CSM67UDR_TR" "CSM79HLA_TR" "CSM79HPA_TR" "MSM79HF9_TR" "MSM9VZNH_TR"
# Remove the missing samples from Metadata
Metadata <- Metadata[!(Metadata$Sample %in% missing_samples), ]

# OPTIONAL scaling Reads using Scale function (to be used as fixed effect later)
Metadata$reads_filtered <- scale(Metadata$reads_filtered)

rownames(Metadata) <- Metadata$Sample
Metadata <- Metadata[, -1]

# transform OTU and taxonomy table into matrices
separated_data <- as.matrix(separated_data)
transformed_data <- as.matrix(transformed_data)

# create phyloseq object
OTU = otu_table(transformed_data, taxa_are_rows = TRUE)
TAX = tax_table(separated_data)
samples = sample_data(Metadata)

humann_phyloseq_object <- phyloseq(OTU, TAX, samples)

# save object 
saveRDS(humann_phyloseq_object, "HMP2_Payami/Phyloseq Objects/HMP2 pathway phyloseq object.rds")
