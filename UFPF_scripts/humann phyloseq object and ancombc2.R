# creating phyloseq object for pathway files 
# DOING THIS ON RPK OUTPUT NOT RELATIVE ABUNDANCE 
library(phyloseq)
library(ANCOMBC)
library(tidyverse)

# ------------------------------------------------------------------------------------------
# humann pathway abundance phyloseq object 
# MAKE TAXONOMY TABLE (which will contain pathways in this case)
load(file = "UFPF/Humann output/consolidated pathways.RData")
pathways <- df2
pathways <- t(pathways)

# filtering out pathways that aren't present in at least 25% of samples 
# Calculate the number of samples
total_samples <- ncol(pathways) 
# Calculate the threshold count level (25% of total samples)
threshold_count <- 0.25 * total_samples
# Filter rows based on this threshold count
filtered_data <- pathways[rowSums(pathways != 0) >= threshold_count, ]
transformed_data <- filtered_data
transformed_data <- as.data.frame(transformed_data)    # filtered from 521 paths to 350

# move row names to a new column 
transformed_data <- transformed_data %>%
  rownames_to_column(var = "Genus")

# Create a new separated data frame with the specified column names. 
# taxonomy is irrelevant for this analysis, but it keeps how I'm creating these phyoloseq objects consistent 
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
# "OTU" is based on the tutorial for 16S data- OTU here is just used as a stand-in for our pathway names 
# adding OTU row names 
labels <- paste0("OTU", seq_len(nrow(transformed_data)))
# Assign labels to row names
rownames(transformed_data) <- labels

# remove column we don't need in the OTU table 
transformed_data <- transformed_data[, -1]

# -------------------------------------------------------------------------------------------
# creating a phyloseq object from the OTU table, Taxonomy table, Metadata table 
Metadata <- readRDS("UFPF/Metadata.rds")

# Scaling Reads using Scale function (to be used as fixed effect later)
Metadata$Reads <- scale(Metadata$Reads)

# transform OTU and taxonomy table into matrices
separated_data <- as.matrix(separated_data)
transformed_data <- as.matrix(transformed_data)

# create the phyloseq object
OTU = otu_table(transformed_data, taxa_are_rows = TRUE)
TAX = tax_table(separated_data)
samples = sample_data(Metadata)

humann_phyloseq_object <- phyloseq(OTU, TAX, samples)

# save object 
saveRDS(humann_phyloseq_object, "UFPF/Phyloseq Objects/humann pathways phyloseq object.rds")


# ran in hipergator 
# running ancombc at the "genus" level -- aka where the pathway info is 
ancom_pathways <- ancombc2(humann_phyloseq_object, 
                           fix_formula = "Sex + Age + Diagnosis2 + Reads",
                           tax_level = "Genus",
                           p_adj_method = "BH",
                           group = "Diagnosis2",
                           lib_cut=0,
                           struc_zero=FALSE,
                           neg_lb=FALSE,
                           alpha=0.05,
                           global = FALSE, 
                           pairwise = TRUE)

ancom_pathways <- readRDS("UFPF/ANCOMBC2/ancombc2 pathways output.rds")

res_pair_pathways = ancom_pathways$res_pair    # 140 paths investigated
res_prim_pathways = ancom_pathways$res

sig_path <- res_pair_pathways %>%               # 132 significant pathways 
  rowwise() %>%
  filter(any(c_across(starts_with("diff_"))))



# -----------------------------------------------------------------------------------------
# creating tables humann ancombc2 data - MetaCyc pathway abundance 
sig_path <- sig_path[, c(1:7, 14:19)]

sig_taxa_long <- sig_path %>%
  pivot_longer(
    cols = starts_with("lfc_"),  
    names_to = "Comparison",    
    values_to = "LFC"           
  )

sig_taxa_long <- sig_taxa_long %>%
  select(-2:-10)

# create separate data frame with significance info 
diff_columns <- sig_path %>%
  select(starts_with("taxon"), starts_with("diff_"))

# Convert to long format
diff_columns <- diff_columns %>%
  pivot_longer(
    cols = starts_with("diff_"),  
    names_to = "Comparison",    
    values_to = "Significance"  
  )

# create separate data frame with q values  
q_columns <- sig_path %>%
  select(starts_with("taxon"), starts_with("q_"))

# Convert to long format
q_columns <- q_columns %>%
  pivot_longer(
    cols = starts_with("q_"),  
    names_to = "Comparison",    
    values_to = "Adj P Value"  
  )

# create separate data frame with standard error values  
se_columns <- sig_path %>%
  select(starts_with("taxon"), starts_with("se_"))

# Convert to long format
se_columns <- se_columns %>%
  pivot_longer(
    cols = starts_with("se_"),  
    names_to = "Comparison",    
    values_to = "Standard_error"  
  )

# adding TRUE/FALSE info 
sig_taxa_long$'Adj P Value' <- q_columns$'Adj P Value'
sig_taxa_long$Significance <- diff_columns$Significance
sig_taxa_long$Standard_error <- se_columns$Standard_error


# Rename groups in the "Comparison" column
sig_taxa_long <- sig_taxa_long %>%
  mutate(Comparison = case_when(
    Comparison == "lfc_Diagnosis2IBD" ~ "IBD vs Control",
    Comparison == "lfc_Diagnosis2PD" ~ "PD vs Control",
    Comparison == "lfc_Diagnosis2PD_Diagnosis2IBD" ~ "PD vs IBD",
    TRUE ~ Comparison  # Keep other values as they are
  ))

# Create three separate data frames for each group- IBD, PD, and PD compared to IBD 
# in this case there are no significant PD pathways 
IBD_data <- sig_taxa_long %>%
  filter(Comparison == "IBD vs Control", Significance == TRUE) %>%    # 93 sig paths 
  arrange(as.numeric(`Adj P Value`))

PD_data <- sig_taxa_long %>%
  filter(Comparison == "PD vs Control", Significance == TRUE) %>%     # 1 sig path
  arrange(as.numeric(`Adj P Value`))

PD_vs_IBD_data <- sig_taxa_long %>%
  filter(Comparison == "PD vs IBD", Significance == TRUE) %>%         # 125 sig paths
  arrange(as.numeric(`Adj P Value`))


# for later 
# create tables sig pathways 
depleted_IBD <- IBD_data[IBD_data$LFC < 0, ]
enriched_IBD <- IBD_data[IBD_data$LFC > 0, ]
depleted_PD <- PD_data[PD_data$LFC < 0, ]
enriched_PD <- PD_data[PD_data$LFC > 0, ]
depleted_PDvIBD <- PD_vs_IBD_data[PD_vs_IBD_data$LFC < 0, ]
enriched_PDvIBD <- PD_vs_IBD_data[PD_vs_IBD_data$LFC > 0, ]

saveRDS(depleted_IBD, "UFPF/ANCOMBC2/depleted pathways IBD.rds") 
saveRDS(enriched_IBD, "UFPF/ANCOMBC2/enriched pathways IBD.rds")
saveRDS(depleted_PD, "UFPF/ANCOMBC2/depleted pathways PD.rds")
saveRDS(enriched_PD, "UFPF/ANCOMBC2/enriched pathways PD.rds")
saveRDS(depleted_PDvIBD, "UFPF/ANCOMBC2/depleted pathways PDvIBD.rds")
saveRDS(enriched_PDvIBD, "UFPF/ANCOMBC2/enriched pathways PDvIBD.rds")
write.csv(depleted_IBD, file = "UFPF/ANCOMBC2/depleted pathways IBD.csv", row.names = FALSE)
write.csv(enriched_IBD, file = "UFPF/ANCOMBC2/enriched pathways IBD.csv", row.names = FALSE)
write.csv(depleted_PD, file = "UFPF/ANCOMBC2/depleted pathways PD.csv", row.names = FALSE)
write.csv(depleted_PDvIBD, file = "UFPF/ANCOMBC2/depleted pathways PDvIBD.csv", row.names = FALSE)
write.csv(enriched_PDvIBD, file = "UFPF/ANCOMBC2/enriched pathways PDvIBD.csv", row.names = FALSE)
# ______________________________________________________________________________

