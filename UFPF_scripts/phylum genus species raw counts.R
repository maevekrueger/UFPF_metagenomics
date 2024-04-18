library(tidyverse)

# generating count tables (counts computed w/the unclassified estimation) 
# at the phylum, genus, and species levels 
# read in raw counts 
counts  <- readRDS("UFPF/Metaphlan output/Counts/Counts w Unclassified.rds")
counts <- column_to_rownames(counts, var = "Sample")

# remove unnecessary column Total Reads
counts <- counts[, -1]
counts <- t(counts)
rownames(counts) <- gsub("\\|", ".", rownames(counts))

# filter out rows based on phylum level 
phylum <- counts[grep("\\.p__[^.]*$", rownames(counts)), ]

phylum_levels <- sub("^.*\\.p__(.+)$", "\\1", rownames(phylum))
rownames(phylum) <- phylum_levels

# filter out rows based on genus level 
genus <- counts[grep("\\.g__[^.]*$", rownames(counts)), ]

genus_levels <- sub("^.*\\.g__(.+)$", "\\1", rownames(genus))
rownames(genus) <- genus_levels

# filter out rows based on species level 
species <- counts[grep("\\.s__[^.]*$", rownames(counts)), ]

species_levels <- sub("^.*\\.s__(.+)$", "\\1", rownames(species))
rownames(species) <- species_levels

# transpose so samples are in row names like originally 
phylym <- t(phylum)
genus <- t(genus)
species <- t(species)

saveRDS(phylum, "UFPF/Metaphlan output/Counts/phylum raw counts.rds")
saveRDS(genus, "UFPF/Metaphlan output/Counts/genus raw counts.rds")
saveRDS(species, "UFPF/Metaphlan output/Counts/species raw counts.rds")

