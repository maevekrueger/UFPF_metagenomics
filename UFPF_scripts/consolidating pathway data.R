library(readr) 
library(tidyverse)

pathway_abundance <- read_tsv("UFPF/Humann output/pathway_abundance_all.tsv") 

# removing extra from the end of all the sample IDs in the table
header <- names(pathway_abundance)
new_header <- str_remove_all(header, "_all_Abundance") %>%
  str_replace_all("-", ".")
names(pathway_abundance) <- new_header

# removing samples with few reads (027), duplicate sample collections, and samples with no metadata(105)
pathway_abundance <- subset(pathway_abundance, select = -c(
  UF.PF.2020.027.1N,
  UF.PF.2022.069.1N_Abundance,
  UF.PF.2022.070.3N_Abundance,
  UF.PF.2022.075.1N_Abundance,
  UF.PF.2022.076.1N_Abundance,
  UF.PF.2022.084.1N2_Abundance,
  UF.PF.2022.085.3N2_Abundance,
  UF.PF.2022.105.1N_Abundance
))

# removing other extra from the end of all the new sample IDs in the table
header <- names(pathway_abundance)
new_header <- str_remove_all(header, "_Abundance") %>%
  str_replace_all("-", ".")
names(pathway_abundance) <- new_header

# transpose so samples are row names
pathway_abundance <- t(pathway_abundance)
pathway_abundance <- as.data.frame(pathway_abundance)
colnames(pathway_abundance) <- pathway_abundance[1, ] 
pathway_abundance <- pathway_abundance[-1, ]

# group into cohorts 
PDp <- pathway_abundance[grep("1N", rownames(pathway_abundance)), ]
Controlp <- pathway_abundance[grep("3N|3G", rownames(pathway_abundance)), ]
IBDp <- pathway_abundance[grep("2G", rownames(pathway_abundance)), ]

all_cohorts_p <- rbind(PDp, Controlp, IBDp)


#~~~~ I use dplyr::select to remove columns whose names contain "UNINTEGRATED" + UNMAPPED \
all_cohorts_p <- all_cohorts_p %>%
  dplyr::select(-matches("UNINTEGRATED"), -UNMAPPED) %>%
  mutate_all(as.numeric)

colnames(all_cohorts_p) <- sub("\\[EC.*", "", colnames(all_cohorts_p))
colnames(all_cohorts_p) <- sub("\\|.*", "", colnames(all_cohorts_p))
colnames(all_cohorts_p)

str(all_cohorts_p)

df2 <- sapply(split.default(all_cohorts_p, names(all_cohorts_p)), rowSums, na.rm = T) %>%
  as.data.frame() 
# %>% mutate(UNMAPPED = all_cohorts_g$UNMAPPED)

names <- colnames(df2)
names <- sub("\\:.*", "", names)

length(unique(names)) == length(names)
#The above line should yield a TRUE in the console
#This means that the number of unique column names in the data frame is equal to 
#the total number of columns
#In other words, each column name is unique

save(df2, file = "UFPF/Humann output/consolidated pathways.RData")

