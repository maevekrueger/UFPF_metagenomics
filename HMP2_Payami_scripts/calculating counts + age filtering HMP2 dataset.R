library(dplyr)
# converting relab to counts 

# Aging Filtering the HMP2 IBD dataset to include subjects 40 and older 
IBD_metadata <- readRDS("HMP2_Payami/old NOT age filtered HMP2 IBD metadata.rds")
# calculate average age 
mean(IBD_metadata$consent_age, na.rm = TRUE)       # 27.47816 years old 

# determine how many samples are at least 40 years old 
sum(IBD_metadata$consent_age >= 40, na.rm = TRUE)    # 337 samples 

# filter to only include samples 40 and over (337 total samples)
IBD_filtered <- IBD_metadata %>%
  filter(consent_age >= 40)

diagnosis_counts <- table(IBD_filtered$diagnosis)
diagnosis_counts    # CD = 101   UC = 97   Control = 139

saveRDS(IBD_filtered, "HMP2_Payami/IBD Metadata Age Filtered.rds")

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

# extract reads_filtered from IBD metadata file 
IBD_reads <- IBD_filtered %>%
  select(Sample, reads_filtered)

bacteria <- bacteria %>%
  rownames_to_column(var = "Sample")
IBD_bacteria_reads <- left_join(IBD_reads, bacteria, by = c("Sample"))

# Create a new data frame to store the updated values
IBD_counts <- IBD_bacteria_reads
rownames(IBD_counts) <- IBD_counts$Sample
IBD_counts$Sample <- NULL

# determining which samples contain all zeros (want to remove them b/c it messes up calculations later)  
row_sums <- rowSums(IBD_counts[, -1])
row_sums <- data.frame(Sample = rownames(IBD_counts), Sum = row_sums)
rownames(row_sums) <- NULL

IBD_counts[, -1] <- sapply(IBD_counts[, -1], as.numeric)
IBD_bacteria_reads[, -1] <- sapply(IBD_bacteria_reads[, -1], as.numeric)


IBD_counts <- IBD_bacteria_reads
# convert relative abundances into counts
# Iterate over each sample in the data frame
for (i in 1:nrow(IBD_bacteria_reads)) {
  # Get the reference value for the current sample
  total <- IBD_bacteria_reads$"reads_filtered"[i]
  
  # Divide the relative abundance values (excluding columns 1 and 2) by 100 and multiply by the reference value
  IBD_counts[i, -(1:2)] <- IBD_bacteria_reads[i, -(1:2)] / 100 * total
}

# new normalized count table
saveRDS(IBD_counts, "HMP2_Payami/HMP2 IBD Age Filtered All Levels Counts.rds")


rownames(IBD_counts) <- IBD_counts$Sample
IBD_counts$Sample <- NULL
IBD_counts <- IBD_counts[, -1]

# filter out columns based on phylum level 
phylum <- IBD_counts[, grep("\\|p__[^|]*$", colnames(IBD_counts))]
colnames(phylum) <- sub("^.*\\|p__", "", colnames(phylum))
# filter out columns based on genus level 
genus <- IBD_counts[, grep("\\|g__[^|]*$", colnames(IBD_counts))]
colnames(genus) <- sub("^.*\\|g__", "", colnames(genus))
# filter out columns based on species level 
species <- IBD_counts[, grep("\\|s__[^|]*$", colnames(IBD_counts))]
colnames(species) <- sub("^.*\\|s__", "", colnames(species))

saveRDS(phylum, "HMP2_Payami/HMP2 IBD Age Filtered phylum Counts.rds")
saveRDS(genus, "HMP2_Payami/HMP2 IBD Age Filtered genus Counts.rds")
saveRDS(species, "HMP2_Payami/HMP2 IBD Age Filtered species Counts.rds")


# -------------------------------------------------------------------------------
# Wallen PD DATA 
PD_counts <- read_excel("HMP2_Payami/PAYAMI DATA- Source_Data_24Oct2022.xlsx", sheet = 3)

# Extract the "clade_name" column
clade_names <- PD_counts$clade_name
PD_counts <- PD_counts[, -1]
rownames(PD_counts) <- clade_names
PD_counts <- t(PD_counts)

# Subsetting Bacteria 
PD_counts <- PD_counts[, grepl("\\_Bacteria", colnames(PD_counts))]

# subset data by cohort
Wallen_PD <- PD_counts[grep("^DP\\d+|^SP\\d+", rownames(PD_counts)), ]
Wallen_Control <- PD_counts[grep("^DC\\d+", rownames(PD_counts)), ]
Wallen_all <- rbind(Wallen_PD, Wallen_Control)

saveRDS(Wallen_all, "HMP2_Payami/Wallen Counts.rds")


# filter out columns based on phylum level 
phylum <- Wallen_all[, grep("\\|p__[^|]*$", colnames(Wallen_all))]
colnames(phylum) <- sub("^.*\\|p__", "", colnames(phylum))
# filter out columns based on genus level 
genus <- Wallen_all[, grep("\\|g__[^|]*$", colnames(Wallen_all))]
colnames(genus) <- sub("^.*\\|g__", "", colnames(genus))
# filter out columns based on species level 
species <- Wallen_all[, grep("\\|s__[^|]*$", colnames(Wallen_all))]
colnames(species) <- sub("^.*\\|s__", "", colnames(species))

saveRDS(phylum, "HMP2_Payami/Wallen_counts_phylum.rds")
saveRDS(genus, "HMP2_Payami/Wallen_counts_genus.rds")
saveRDS(species, "HMP2_Payami/Wallen_counts_species.rds")
