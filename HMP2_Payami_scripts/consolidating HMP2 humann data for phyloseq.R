# cleaning up HMP2 functional data for ancombc2 
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr) 
library(stringr) 

# ran this in hipergator 
path1 <- "HMP2_Payami/genefamilies.tsv.gz"
gzipped1 <- gzfile(path1, "rt")
# reading in the column names was being funky, but this way seems to work: 
first_row1 <- readLines(gzipped1, n = 1)
# Split the first row into column names
column_names1 <- unlist(strsplit(first_row1, "\t"))
IBD_gene <- read.table(gzipped1, header = FALSE, sep = "\t", col.names = column_names1)
close(gzipped1)


# genefamillies file was then further proccessed in hipergator 
# 1- regrouped to KO groups 
# 2- then renamed KO groups 

IBD_KO <- read_tsv("HMP2_Payami/hmp2_KO_groups_renamed.tsv")

# clean up column names 
colnames(IBD_KO)[1] <- "KO Groups"
header <- names(IBD_KO)
new_header <- gsub("_Abundance.RPKs", "", header)
names(IBD_KO) <- new_header


# filter so only contains age filtered IBD samples 
Metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")
KO_t <- t(IBD_KO)
KO_t <- as.data.frame(`colnames<-`(KO_t[-1, ], KO_t[1, ]))
KO_t <- rownames_to_column(KO_t, var = "Sample")
samples <- Metadata$Sample
filtered_KO <- KO_t[KO_t$Sample %in% samples, ]

# Note: Metadata contains 337 samples but KO groups only contains 332
# Identify samples in Metadata that are missing from filtered_KO 
missing_samples <- setdiff(Metadata$Sample, filtered_KO$Sample)
missing_samples 
# ^ "CSM67UDR_TR" "CSM79HLA_TR" "CSM79HPA_TR" "MSM79HF9_TR" "MSM9VZNH_TR"

saveRDS(filtered_KO, "HMP2_Payami/HMP2 KO Groups Age Filtered.rds")


# transpose so samples are row names
#genefamilies_ko <- t(genefamilies_ko)
#genefamilies_ko <- as.data.frame(genefamilies_ko)
#colnames(genefamilies_ko) <- genefamilies_ko[1, ] 
#genefamilies_ko <- genefamilies_ko[-1, ]

rownames(filtered_KO) <- NULL
filtered_KO <- column_to_rownames(filtered_KO, var = "Sample")

#~~~~ I use dplyr::select to remove columns whose names contain "UNGROUPED" + UNMAPPED \
filtered_KO <- filtered_KO %>%
  dplyr::select(-matches("UNGROUPED"), -UNMAPPED) %>%
  mutate_all(as.numeric)

colnames(filtered_KO) <- sub("\\[EC.*", "", colnames(filtered_KO))
colnames(filtered_KO) <- sub("\\|.*", "", colnames(filtered_KO))
colnames(filtered_KO)

str(filtered_KO)

df <- sapply(split.default(filtered_KO, names(filtered_KO)), rowSums, na.rm = T) %>%
  as.data.frame() 

names <- colnames(df)
names <- sub("\\:.*", "", names)

length(unique(names)) == length(names)
#The above line should yield a TRUE in the console
#This means that the number of unique column names in the data frame is equal to 
#the total number of columns
#In other words, each column name is unique

save(df, file = "HMP2_Payami/HMP2 Consolidated and Age Filtered KO Groups.RData")



# -------------------------------------------------------------------------------
# pathways 
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

`#  consolidate the identical pathway names 
paths <- sapply(split.default(IBD_path_filtered, names(IBD_path_filtered)), rowSums, na.rm = T) %>%
  as.data.frame() 

names <- colnames(paths)
names <- sub("\\:.*", "", names)

length(unique(names)) == length(names)       # <- TRUE

any(duplicated(rownames(names)))             # <- FALSE 


save(paths, file = "Hmp2_Payami/HMP2 Consolidated and Age Filtered Pathways.RData")
