library(ggplot2)
library(tibble)
library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)

# HMP2 IBD data IBD combined 
ancom_IBD <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD species ancombc2 output.rds")

res_prim_IBD = ancom_IBD$res
sig_taxa_IBD <- res_prim_IBD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))  # 63 sig species 

sig_taxa_IBD <- sig_taxa_IBD[, c(1, 3, 19, 23)]


# Payami PD data 
ancom_PD <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD species ancombc2 output.rds")

res_prim_PD = ancom_PD$res
sig_taxa_PD <- res_prim_PD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))  # 11 sig species

sig_taxa_PD <- sig_taxa_PD[, c(1, 3, 19, 23)]


# subsetting enriched species 
enriched_IBD <- sig_taxa_IBD[sig_taxa_IBD$lfc_diagnosis2IBD > 0, ]
enriched_PD <- sig_taxa_PD[sig_taxa_PD$lfc_Case_statusPD > 0, ]

# subsetting depleted species 
depleted_IBD <- sig_taxa_IBD[sig_taxa_IBD$lfc_diagnosis2IBD < 0, ]
depleted_PD <- sig_taxa_PD[sig_taxa_PD$lfc_Case_statusPD < 0, ]

# Extract enriched species names 
enriched_species_IBD <- unique(enriched_IBD$taxon)
enriched_species_PD <- unique(enriched_PD$taxon)

# Extract depleted species names 
depleted_species_IBD <- unique(depleted_IBD$taxon)
depleted_species_PD <- unique(depleted_PD$taxon)


# using ggVennDiagram 
x <- list(
  "HMP2 IBD" = depleted_species_IBD,
  "Wallen PD" = depleted_species_PD
)

venn <- Venn(x)
data <- process_data(venn)

# CUSTOMIZE GROUP COLORS
colorGroups <- c("HMP2 IBD"= "green", "Wallen PD" = "blue")  

colfunc <- colorRampPalette(colorGroups)
col <- colfunc(3)

ggplot() +
  geom_sf(aes(fill = name), data = venn_region(data), show.legend = FALSE) +
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(aes(label = name), size = 16, data = venn_setlabel(data)) +
  geom_sf_text(aes(label = count), size = 20, data = venn_region(data)) +
  scale_fill_manual(values = alpha(col, .25)) +
  scale_color_manual(values = col) +
  theme_void()

ggsave("HMP2_Payami/Figures/HMP2 Wallen Venn Species Depleted.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")



# ENRICHED
# using ggVennDiagram 
x <- list(
  "HMP2 IBD" = enriched_species_IBD,
  "Wallen PD" = enriched_species_PD
)

venn <- Venn(x)
data <- process_data(venn)

# CUSTOMIZE GROUP COLORS
colorGroups <- c("UFPF IBD"= "green", "UFPF PD" = "blue")  

colfunc <- colorRampPalette(colorGroups)
col <- colfunc(3)

ggplot() +
  geom_sf(aes(fill = name), data = venn_region(data), show.legend = FALSE) +
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(aes(label = name), size = 16, data = venn_setlabel(data)) +
  geom_sf_text(aes(label = count), size = 20, data = venn_region(data)) +
  scale_fill_manual(values = alpha(col, .25)) +
  scale_color_manual(values = col) +
  theme_void()

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 Wallen Venn Species Enriched.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")



# PATHWAYS ----------------------------------------------------------
# HMP2 IBD 
paths <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD combined pathway ancombc2 output.rds")
results <- paths$res     

sig_path <- results %>%               # 158 significant pathways 
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))

# Create a data frame for IBD significant pathways
ibd_sig_pathways <- sig_path[, c(1, 3, 7, 11, 15, 19)]
ibd_sig_pathways <- ibd_sig_pathways[order(ibd_sig_pathways$q_diagnosis2IBD), ]

ibd_enriched <- ibd_sig_pathways[ibd_sig_pathways$lfc_diagnosis2IBD > 0, ]   # 96
ibd_depleted <- ibd_sig_pathways[ibd_sig_pathways$lfc_diagnosis2IBD < 0, ]   # 62


# Wallen PD paths 
paths <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD ancombc2 metacyc pathways.rds")
results <- paths$res     

sig_path <- results %>%               # 198 significant pathways 
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))

# Create a data frame for IBD significant pathways
pd_sig_pathways <- sig_path[, c(1, 3, 7, 11, 15, 19)]
pd_sig_pathways <- pd_sig_pathways[order(pd_sig_pathways$q_Case_statusPD), ]

pd_enriched <- pd_sig_pathways[pd_sig_pathways$lfc_Case_statusPD > 0, ]   # 46
pd_depleted <- pd_sig_pathways[pd_sig_pathways$lfc_Case_statusPD < 0, ]   # 152


# identifying shared depleted pathways for later 
common_depleted_paths <- merge(ibd_depleted, pd_depleted, by = "taxon")   # 26
# identifying shared enriched pathways 
common_enriched_paths <- merge(ibd_enriched, pd_enriched, by = "taxon")   # 8 
# pathways they shared but different in enriched/depleted
common_paths <- merge(ibd_sig_pathways, pd_sig_pathways, by = "taxon")   # 53-26-8 =(19 paths)


# remove unmapped or unintegrated 
ibd_depleted <- ibd_depleted %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED", taxon))
ibd_enriched <- ibd_enriched %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED", taxon))
pd_depleted <- pd_depleted %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED", taxon))
pd_enriched <- pd_enriched %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED", taxon))


# Extract enriched pathway names 
enriched_paths_IBD <- unique(ibd_enriched$taxon)
enriched_paths_PD <- unique(pd_enriched$taxon)

# Extract depleted pathway names 
depleted_paths_IBD <- unique(ibd_depleted$taxon)
depleted_paths_PD <- unique(pd_depleted$taxon)


# using ggVennDiagram 
x <- list(
  "HMP2 IBD" = depleted_paths_IBD,
  "Wallen PD" = depleted_paths_PD
)

venn <- Venn(x)
data <- process_data(venn)

# CUSTOMIZE GROUP COLORS
colorGroups <- c("HMP2 IBD"= "green", "Wallen PD" = "blue")  

colfunc <- colorRampPalette(colorGroups)
col <- colfunc(3)

ggplot() +
  geom_sf(aes(fill = name), data = venn_region(data), show.legend = FALSE) +
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(aes(label = name), size = 16, data = venn_setlabel(data)) +
  geom_sf_text(aes(label = count), size = 20, data = venn_region(data)) +
  scale_fill_manual(values = alpha(col, .25)) +
  scale_color_manual(values = col) +
  theme_void()

ggsave("HMP2_Payami/Figures/HMP2 Wallen Venn Pathways Depleted.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")


# ENRICHED 
# using ggVennDiagram 
x <- list(
  "HMP2 IBD" = enriched_paths_IBD,
  "Wallen PD" = enriched_paths_PD
)

venn <- Venn(x)
data <- process_data(venn)

# CUSTOMIZE GROUP COLORS
colorGroups <- c("HMP2 IBD"= "green", "Wallen PD" = "blue")  

colfunc <- colorRampPalette(colorGroups)
col <- colfunc(3)

ggplot() +
  geom_sf(aes(fill = name), data = venn_region(data), show.legend = FALSE) +
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(aes(label = name), size = 16, data = venn_setlabel(data)) +
  geom_sf_text(aes(label = count), size = 20, data = venn_region(data)) +
  scale_fill_manual(values = alpha(col, .25)) +
  scale_color_manual(values = col) +
  theme_void()

ggsave("HMP2_Payami/Figures/HMP2 Wallen Venn Pathways Enriched.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")
