library(readxl)
library(tidyr)
library(dplyr)
library(tibble)

# cleaning Wallen (payami) metaphlan relative abundance output 
payami_counts <- read_excel("HMP2_Payami/PAYAMI DATA- Source_Data_24Oct2022.xlsx", sheet = 3)
payami_taxa <- read_excel("HMP2_Payami/PAYAMI DATA- Source_Data_24Oct2022.xlsx", sheet = 4)

# Extract the "clade_name" column
clade_names <- payami_taxa$clade_name
payami_taxa <- payami_taxa[, -1]
rownames(payami_taxa) <- clade_names
payami_taxa <- t(payami_taxa)

# subset data by cohort
uab_PD <- payami_taxa[grep("^DP\\d+|^SP\\d+", rownames(payami_taxa)), ]
uab_Control <- payami_taxa[grep("^DC\\d+", rownames(payami_taxa)), ]

payami_all <- rbind(uab_PD, uab_Control)

# Subsetting Bacteria 
bacteria <- payami_all[, grepl("\\_Bacteria", colnames(payami_all))]

# filter out columns based on phylum level 
rel_ab_phylum <- bacteria[, grep("\\|p__[^|]*$", colnames(bacteria))]
colnames(rel_ab_phylum) <- sub("^.*\\|p__", "", colnames(rel_ab_phylum))
# filter out columns based on genus level 
rel_ab_genus <- bacteria[, grep("\\|g__[^|]*$", colnames(bacteria))]
colnames(rel_ab_genus) <- sub("^.*\\|g__", "", colnames(rel_ab_genus))
# filter out columns based on species level 
rel_ab_species <- bacteria[, grep("\\|s__[^|]*$", colnames(bacteria))]
colnames(rel_ab_species) <- sub("^.*\\|s__", "", colnames(rel_ab_species))

# make it a data frame 
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
saveRDS(rel_ab_phylum, "HMP2_Payami/uAB_phylum.rds")
saveRDS(rel_ab_genus, "HMP2_Payami/uAB_ab_genus.rds")
saveRDS(rel_ab_species, "HMP2_Payami/uAB_species.rds")
