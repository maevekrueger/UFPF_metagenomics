library(tidyverse)

# computing individual counts from relative abundance data generated from MetaPhlAn
# we primarily relied on using calculated counts that included the unclassified estimations in MetaPhlAn

# getting counts from RAs with the unclassified flag added 
rel_ab_unkn <- read.table("UFPF/Metaphlan output/metaphlan_unknown_all_combined.tsv", header = TRUE, row.names = 1, sep = "\t")

# removing "_metaphlan_rel_ab" from the end of all the sample IDs in the table
header <- names(rel_ab_unkn)
new_header <- gsub("_metaphlan_rel_ab_unknown", "", header)
names(rel_ab_unkn) <- new_header

# removing sample 2020-027-1N that had zero reads and 105 with no Metadata
# also removing samples with duplicate sample collection (only using the first sample a patient submitted)
rel_ab_unkn <- rel_ab_unkn %>%
  select(-c("UF.PF.2020.027.1N", "UF.PF.2022.105.1N", "UF.PF.2022.084.1N2", 
            "UF.PF.2022.085.3N2"))
rel_ab_unkn <- as.data.frame(t(rel_ab_unkn))

rel_ab_unkn <- rel_ab_unkn %>%
  rownames_to_column(var = "Sample")
# fixing mislabeled sample: changing 2021-023-1N to 2020-023-1N 
rel_ab_unkn <- rel_ab_unkn %>%
  mutate(
    Sample = case_when(
      Sample == "UF.PF.2021.023.1N" ~ "UF.PF.2020.023.1N",
      TRUE ~ Sample
    )
  ) %>%
  arrange(as.numeric(gsub("\\D", "", Sample)))

rel_ab_unkn <- column_to_rownames(rel_ab_unkn, var = "Sample")

saveRDS(rel_ab_unkn, "UFPF/Metaphlan output/rel_ab_unkn_cleaned.rds")



rel_ab_unkn <- readRDS("UFPF/Metaphlan output/rel_ab_unkn_cleaned.rds")
rel_ab_unkn <- rownames_to_column(rel_ab_unkn, "Sample")

# add in total read count that's located in the Metadata file 
Metadata<- readRDS("UFPF/Metadata.rds")
nreads <- data.frame(Sample = rownames(Metadata), Reads = Metadata$Reads)

rel_ab_t_nreads <- left_join(nreads, rel_ab_unkn, by = c("Sample"))

# Create a new data frame to store the updated values
counts_data2 <- rel_ab_t_nreads

# create a function to compute the individual counts from the RA values
for (i in 1:nrow(rel_ab_t_nreads)) {
  # Get the Total Read Count for the current sample
  total <- rel_ab_t_nreads$"Reads"[i]
  
  # Divide the RA values (excluding columns 1 and 2) by 100 and multiply by the Total Read Count for that sample
  counts_data2[i, -(1:2)] <- rel_ab_t_nreads[i, -(1:2)] / 100 * total
}

# rename nreads column
colnames(counts_data2)[colnames(counts_data2) == "Reads"] <- "Total Reads"

saveRDS(counts_data2, "UFPF/Metaphlan output/Counts/Counts w Unclassified.rds")



# -----------------------------------------------------------------------------------------
# without the unclassified flag added 
# reading .tsv file containing taxonomic classification relative abundances from MetaPhlAn
rel_ab <- read.table("UFPF/Metaphlan output/metaphlan_all_combined.tsv", header = TRUE, row.names = 1, sep = "\t")

# removing "_metaphlan_rel_ab" from the end of all the sample IDs in the table
header <- names(rel_ab)
new_header <- gsub("_metaphlan_rel_ab", "", header)
names(rel_ab) <- new_header

# removing sample 2020-027-1N that had zero reads 
rel_ab <- rel_ab %>%
  select(-c("UF.PF.2020.027.1N", "UF.PF.2022.105.1N", "UF.PF.2022.084.1N2", 
            "UF.PF.2022.085.3N2"))

rel_ab_transposed <- as.data.frame(t(rel_ab))
rel_ab_transposed <- rel_ab_transposed %>%
  rownames_to_column(var = "Sample")

# changing 2021-023-1N to 2020-023-1N 
rel_ab_transposed <- rel_ab_transposed%>%
  mutate(Sample = case_when(
    Sample == "UF.PF.2021.023.1N" ~ "UF.PF.2020.023.1N",
    TRUE ~ Sample))

rel_ab_transposed <- column_to_rownames(rel_ab_transposed, var = "Sample")
saveRDS(rel_ab_transposed, "UFPF/Metaphlan output/rel_ab_cleaned.rds")


# Read in Metadata with total read count for each sample 
Metadata<- readRDS("UFPF/Metadata.rds")
nreads <- data.frame(Sample = rownames(Metadata), Reads = Metadata$Reads)

rel_ab_transposed <- rownames_to_column(rel_ab_transposed, "Sample")
rel_ab_t_nreads <- left_join(nreads, rel_ab_transposed, by = c("Sample"))

# Create a new data frame to store the updated values
counts_data <- rel_ab_t_nreads

# Create a function to calculate the individual counts from RA values 
for (i in 1:nrow(rel_ab_t_nreads)) {
  # Get the Total Read Count for the current sample
  total <- rel_ab_t_nreads$"Reads"[i]
  
  # Divide the relative abundance values (excluding columns 1 and 2) by 100 and multiply by the total read count 
  counts_data[i, -(1:2)] <- rel_ab_t_nreads[i, -(1:2)] / 100 * total
}

# rename nreads column
colnames(counts_data)[colnames(counts_data) == "Reads"] <- "Total Reads"

saveRDS(counts_data, "UFPF/Metaphlan output/Counts/Counts WITHOUT unclassified.rds")
