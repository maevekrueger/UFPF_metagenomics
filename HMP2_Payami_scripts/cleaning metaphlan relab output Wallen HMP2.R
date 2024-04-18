library(readxl)
library(tidyverse)

# cleaning Wallen metaphlan relative abundance output 
# (Wallen might also be referred to as Payami)
Wallen_counts <- read_excel("HMP2_Payami/PAYAMI DATA- Source_Data_24Oct2022.xlsx", sheet = 3)
Wallen_taxa <- read_excel("HMP2_Payami/PAYAMI DATA- Source_Data_24Oct2022.xlsx", sheet = 4)

# Extract the "clade_name" column
clade_names <- Wallen_taxa$clade_name
Wallen_taxa <- Wallen_taxa[, -1]
rownames(Wallen_taxa) <- clade_names
Wallen_taxa <- t(Wallen_taxa)

# subset data by cohort
Wallen_PD <- Wallen_taxa[grep("^DP\\d+|^SP\\d+", rownames(Wallen_taxa)), ]
Wallen_Control <- Wallen_taxa[grep("^DC\\d+", rownames(Wallen_taxa)), ]

Wallen_all <- rbind(Wallen_PD, Wallen_Control)

# Subsetting Bacteria 
bacteria <- Wallen_all[, grepl("\\_Bacteria", colnames(Wallen_all))]

# filter out columns based on phylum level 
rel_ab_phylum <- bacteria[, grep("\\|p__[^|]*$", colnames(bacteria))]
colnames(rel_ab_phylum) <- sub("^.*\\|p__", "", colnames(rel_ab_phylum))
# filter out columns based on genus level 
rel_ab_genus <- bacteria[, grep("\\|g__[^|]*$", colnames(bacteria))]
colnames(rel_ab_genus) <- sub("^.*\\|g__", "", colnames(rel_ab_genus))
# filter out columns based on species level 
rel_ab_species <- bacteria[, grep("\\|s__[^|]*$", colnames(bacteria))]
colnames(rel_ab_species) <- sub("^.*\\|s__", "", colnames(rel_ab_species))

rel_ab_phylum <- as.data.frame(rel_ab_phylum)
rel_ab_genus <- as.data.frame(rel_ab_genus)
rel_ab_species <- as.data.frame(rel_ab_species)

# adjusting RA values at each taxonomic level so they are all out of 100 
# Calculate the total for each row (excluding Sample column)
phylum_sum <- rowSums(rel_ab_phylum)
rel_ab_phylum <- rel_ab_phylum / phylum_sum * 100

# CEHCK  *good to go 
row_totals <- rowSums(rel_ab_phylum)
all_totals_100 <- all(abs(row_totals - 100) < 1e-6)
# Print the result
if (all_totals_100) {
  print("All rows sum to 100 after rescaling.")
} else {
  print("Not all rows sum to 100 after rescaling.")
}

genus_sum <- rowSums(rel_ab_genus)
rel_ab_genus <- rel_ab_genus / genus_sum * 100

species_sum <- rowSums(rel_ab_species)
rel_ab_species <- rel_ab_species / species_sum * 100


# SAVE
saveRDS(rel_ab_phylum, "HMP2_Payami/Wallen_relab_phylum.rds")
saveRDS(rel_ab_genus, "HMP2_Payami/Wallen_relab_genus.rds")
saveRDS(rel_ab_species, "HMP2_Payami/Wallen_relab_species.rds")



# ------------------------------------------------------------------
# HMP2 IBD dataset 
# getting relative abundances at phylum, genus, species-level
IBD_filtered <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")

# loading in IBD taxanomic data (have to decompress file)
path <- "HMP2_Payami/taxonomic_profiles.tsv.gz"
gzipped_file <- gzfile(path, "rt")
# reading in the column names was being funky, but this way seems to work: 
first_row <- readLines(gzipped_file, n = 1)
# Split the first row into column names
column_names <- unlist(strsplit(first_row, "\t"))
IBD_taxa <- read.table(gzipped_file, header = FALSE, sep = "\t", col.names = column_names)
close(gzipped_file)

# clean up column names 
header <- names(IBD_taxa)
new_header <- gsub("_taxonomic_profile", "", header)
names(IBD_taxa) <- new_header
colnames(IBD_taxa)[1] <- "Sample"

# there are 1339 columns- just checking to ensure no duplicates -- *none found  
duplicate_columns <- colnames(IBD_taxa)[duplicated(colnames(IBD_taxa))]

if (length(duplicate_columns) > 0) {
  cat("Duplicate column names found:", paste(duplicate_columns, collapse = ", "), "\n")
} else {
  cat("No duplicate column names found.\n")
}

# cleaning 
IBD_taxa_t <- as.data.frame(t(IBD_taxa))
colnames(IBD_taxa_t) <- unlist(IBD_taxa_t[1, ])
IBD_taxa_t <- IBD_taxa_t[-1, ]

IBD_taxa_t <- cbind("Sample" = rownames(IBD_taxa_t), IBD_taxa_t)
rownames(IBD_taxa_t) <- NULL

# filter the taxonomic data to match the filtered metadata file 
filtered_IBD_samples <- unique(IBD_filtered$Sample)
# Filter IBD taxa based on matching sample names
filtered_IBD_age <- IBD_taxa_t %>%
  filter(Sample %in% filtered_IBD_samples)

# Extract "Sample" and "diagnosis" and "diagnosis2" columns for later 
diagnosis <- IBD_filtered %>% select(Sample, diagnosis, diagnosis2)

saveRDS(diagnosis, "HMP2_Payami/HMP2 IBD diagnosis dataframe.rds")

filtered_IBD_age <- filtered_IBD_age %>%
  column_to_rownames(var = "Sample")

# Subsetting Bacteria 
bacteria <- filtered_IBD_age[, grepl("\\_Bacteria", colnames(filtered_IBD_age))]

# saving relative abundace output 
saveRDS(bacteria, "HMP2_Payami/HMP2 IBD Age Filtered All Levels relab.rds")

# filter out columns based on phylum level 
phylum <- bacteria[, grep("\\|p__[^|]*$", colnames(bacteria))]
colnames(phylum) <- sub("^.*\\|p__", "", colnames(phylum))
# filter out columns based on genus level 
genus <- bacteria[, grep("\\|g__[^|]*$", colnames(bacteria))]
colnames(genus) <- sub("^.*\\|g__", "", colnames(genus))
# filter out columns based on species level 
species <- bacteria[, grep("\\|s__[^|]*$", colnames(bacteria))]
colnames(species) <- sub("^.*\\|s__", "", colnames(species))

saveRDS(phylum, "HMP2_Payami/HMP2 IBD Age Filtered phylum relab.rds") 
saveRDS(genus, "HMP2_Payami/HMP2 IBD Age Filtered genus relab.rds")
saveRDS(species, "HMP2_Payami/HMP2 IBD Age Filtered species relab.rds")