# GENUS LEVEL 
# Making phyloseq objects with scaled reads in metadata file 
IBD <- readRDS("HMP2_Payami/HMP2 IBD Age Filtered All Levels Counts.rds")
IBD_t <- as.data.frame(t(IBD))

IBD_t <- IBD_t %>%
  slice(2:n()) %>%
  setNames(unlist(IBD_t[1, ]))
IBD_t <- IBD_t[-1, ]

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
transformed_data2 <- data.frame(
  Taxonomy = extract_last_level(rownames(IBD_t)),
  Classification_Level = extract_classification_level(rownames(IBD_t)),
  IBD_t
)

# isolate genus names from row names 
transformed_data2 <- transformed_data2[grep("\\|g__[^|]*$", rownames(transformed_data2)), ]

transformed_data_levels2 <- sub("^.+\\|g__(.+)$", "\\1", rownames(transformed_data2))
rownames(transformed_data2) <- transformed_data_levels2

# filtering out genera that aren't present in at least 5% of samples 
# Calculate the number of samples
total_samples2 <- ncol(transformed_data2) - 2  # 337 samples 
# Calculate the threshold count level (5% of total samples)
threshold_count2 <- 0.05 * total_samples2      # threshold = 16.85 samples 
# Filter rows based on the threshold count
filtered_data2 <- transformed_data2[rowSums(transformed_data2[, -c(1, 2)] != 0) >= threshold_count2, ]
transformed_data2 <- filtered_data2            # 142 --> now filtered to 106 genera 

# move row names to a new column 
transformed_data2 <- transformed_data2 %>%
  rownames_to_column(var = "Genus")

# remove unnecessary column 
transformed_data2 <- transformed_data2[, -2]

# Separating taxonomic levels into separate columns 
separated_data2 <- transformed_data2 %>%
  separate(Classification_Level, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = "\\|")

# clean up 
separated_data2 <- apply(separated_data2, 2, function(col) gsub(".*__", "", col))
separated_data2 <- as.data.frame(separated_data2)
separated_data2 <- separated_data2[, c(1:6)]

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(separated_data2)))
# Assign labels to row names
rownames(separated_data2) <- labels

# MAKE OTU TABLE 
# remove unnecessary column 
transformed_data2 <- transformed_data2[, -2]

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data2)))
# Assign labels to row names
rownames(transformed_data2) <- labels

# remove columns we don't need in the OTU table 
transformed_data2 <- transformed_data2[, -1]


# -------------------------------------------------------------------------------------------
# creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
IBD_metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")    # metadata age filtered

IBD_metadata <- IBD_metadata %>%
  column_to_rownames(var = "Sample")

# scaling total_sequences using Scale function (to be used as fixed effect later)
IBD_metadata$reads_filtered <- as.vector(scale(IBD_metadata$reads_filtered))

# transform OTU and taxonomy table into matrices
separated_data2 <- as.matrix(separated_data2)
transformed_data2 <- as.matrix(transformed_data2)

transformed_data2 <- as.data.frame(as.matrix(transformed_data2)) %>%
  mutate_all(as.numeric)

# create phyloseq object
OTU2 = otu_table(transformed_data2, taxa_are_rows = TRUE)
TAX2 = tax_table(separated_data2)
samples2 = sample_data(IBD_metadata)

phyloseq_object <- phyloseq(OTU2, TAX2, samples2)

# save object 
saveRDS(phyloseq_object, "HMP2_Payami/HMP2 IBD genus phyloseq object.rds")




# Payami PD Data 
PD <- readRDS("HMP2_Payami/Payami all levels Counts.rds")
# move taxa into rownames 
PD_t <- t(PD)

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
  Taxonomy = extract_last_level(rownames(PD_t)),
  Classification_Level = extract_classification_level(rownames(PD_t)),
  PD_t
)

# isolate species names from row names 
transformed_data <- transformed_data[grep("\\|g__[^|]*$", rownames(transformed_data)), ]

transformed_data_levels <- sub("^.+\\|g__(.+)$", "\\1", rownames(transformed_data))
rownames(transformed_data) <- transformed_data_levels

# filtering out genera that aren't present in at least 5% of samples 
# Calculate the number of samples
total_samples <- ncol(transformed_data) - 2  # 722 samples 
# Calculate the threshold count level (5% of total samples)
threshold_count <- 0.05 * total_samples      # threshold = 36.1 samples 
# Filter rows based on the threshold count
filtered_data <- transformed_data[rowSums(transformed_data[, -c(1, 2)] != 0) >= threshold_count, ]
transformed_data <- filtered_data            # 216 genera --> now filtered to 105 

# move row names to a new column 
transformed_data <- transformed_data %>%
  rownames_to_column(var = "Genus")

# remove unnecessary column 
transformed_data <- transformed_data[, -2]

# Separating taxonomic levels into separate columns 
separated_data <- transformed_data %>%
  separate(Classification_Level, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = "\\|")

# clean up 
separated_data <- apply(separated_data, 2, function(col) gsub(".*__", "", col))
separated_data <- as.data.frame(separated_data)
separated_data <- separated_data[, c(1:6)]

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(separated_data)))
# Assign labels to row names
rownames(separated_data) <- labels


# MAKE OTU TABLE 
# remove unnecessary column 
transformed_data <- transformed_data[, -2]

# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data)))
# Assign labels to row names
rownames(transformed_data) <- labels

# remove columns we don't need in the OTU table 
transformed_data <- transformed_data[, -1]

# -------------------------------------------------------------------------------------------
# creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
PD_metadata <- readRDS("HMP2_Payami/Payami PD Metadata.rds")

# remove unnecessary column
PD_metadata <- PD_metadata[, -which(names(PD_metadata) == "Sample.1")]

# scaling total_sequences using Scale function (to be used as fixed effect later)
PD_metadata$total_sequences <- as.vector(scale(PD_metadata$total_sequences))

# transform OTU and taxonomy table into matrices
separated_data <- as.matrix(separated_data)
transformed_data <- as.matrix(transformed_data)

# create phyloseq object
OTU = otu_table(transformed_data, taxa_are_rows = TRUE)
TAX = tax_table(separated_data)
samples = sample_data(PD_metadata)

phyloseq_object2 <- phyloseq(OTU, TAX, samples)

# save object 
saveRDS(phyloseq_object2, "HMP2_Payami/Phyloseq Objects/Wallen PD genus phyloseq object.rds")

