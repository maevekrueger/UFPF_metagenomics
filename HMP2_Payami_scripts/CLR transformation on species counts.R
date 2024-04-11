# CLR transformation on raw counts 
library(compositions)
library(tidyverse)

# read in species-level counts 
IBD <- readRDS("HMP2_Payami/HMP2 IBD Age Filtered species Count.rds")
PD <- readRDS("HMP2_Payami/Wallen_counts_species.rds")

# add a small pseudocount to avoid any zeros in the df which clr doesn't like 
IBD <- IBD + 1e-10
PD <- PD + 1e-10

# perform clr transformation with compositions package 
clr_IBD <- as.data.frame(clr(IBD))
clr_PD <- as.data.frame(clr(PD))

# Wallen PD DATA
payami_all <- cbind(
  data.frame(Sample = rownames(clr_PD),
             diagnosis = ifelse(grepl("^DP\\d+|^SP\\d+", rownames(clr_PD)), "PD",
                                ifelse(grepl("^DC\\d+", rownames(clr_PD)), "Control", NA))),
  clr_PD
)
rownames(payami_all) <- NULL

saveRDS(payami_all, "HMP2_Payami/Wallen_species_clr_counts.rds")



# HMP2 IBD DATA
IBD_metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")
# Extract "Sample" and "diagnosis" 
diagnosis2 <- IBD_metadata %>% select(Sample, diagnosis2)

clr_IBD$Sample <- rownames(clr_IBD)
clr_IBD <- clr_IBD[, c("Sample", setdiff(names(clr_IBD), "Sample"))]
rownames(clr_IBD) <- NULL

IBD <- merge(clr_IBD, diagnosis2, by = "Sample")
column_order <- c("Sample", "diagnosis2", setdiff(names(IBD), c("Sample", "diagnosis2")))
IBD <- IBD[, column_order]

ibd <- IBD[IBD$diagnosis2 == "IBD", ]               # 198 IBD 
nonIBD <- IBD[IBD$diagnosis2 == "nonIBD", ]        # 139 controls (nonIBD)
all <- rbind(ibd, nonIBD)

saveRDS(all, "HMP2_Payami/HMP2 IBD Combined Age Filtered species clr counts.rds")
