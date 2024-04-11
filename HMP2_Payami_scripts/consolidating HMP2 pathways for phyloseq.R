# cleaning up HMP2 functional data for ancombc2 
library(tidyverse)

# read in pathways file 
df <- read_tsv("HMP2_Payami/pathabundances.tsv.gz") 

df$`# Pathway`[3703] #this index was previously a problematic one - now see that the length is fine 

df <- df %>%
  dplyr::rename("Pathway" = "# Pathway")

colnames(df) <- str_remove_all(colnames(df), "_Abundance")

df <- df %>%
  column_to_rownames(var = "Pathway")

df <- t(df) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")

Metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")
samples <- Metadata$Sample
filtered_path <- df[df$Sample %in% samples, ]
IBD_path_filtered <- filtered_path

saveRDS(IBD_path_filtered, "HMP2_Payami/HMP2 Pathways Age Filtered.rds")

# Want to consolidate identical pathways and sum their abundances 
# I previously cleaned this up and age filtered so it only includes the subjects I want 
IBD_path_filtered <- readRDS("HMP2_Payami/HMP2 Pathways Age Filtered.rds")

rownames(IBD_path_filtered) <- NULL
IBD_path_filtered <- column_to_rownames(IBD_path_filtered, var = "Sample")

colnames(IBD_path_filtered) <- sub("\\|.*", "", colnames(IBD_path_filtered))
colnames(IBD_path_filtered) %>%
  unique()

str(IBD_path_filtered)

#  consolidate the identical pathway names 
paths <- sapply(split.default(IBD_path_filtered, names(IBD_path_filtered)), rowSums, na.rm = T) %>%
  as.data.frame() 

names <- colnames(paths)
names <- sub("\\:.*", "", names)

length(unique(names)) == length(names)       # <- TRUE

any(duplicated(rownames(names)))             # <- FALSE 


save(paths, file = "Hmp2_Payami/HMP2 Consolidated and Age Filtered Pathways.RData")

