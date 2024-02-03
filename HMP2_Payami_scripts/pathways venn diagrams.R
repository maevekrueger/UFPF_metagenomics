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

write.csv(ibd_enriched, file = "HMP2_Payami/ANCOMBC2/enriched HMP2 IBD paths.csv", row.names = FALSE)
write.csv(ibd_depleted, file = "HMP2_Payami/ANCOMBC2/depleted HMP2 IBD paths.csv", row.names = FALSE)


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

write.csv(pd_enriched, file = "HMP2_Payami/ANCOMBC2/enriched Wallen PD paths.csv", row.names = FALSE)
write.csv(pd_depleted, file = "HMP2_Payami/ANCOMBC2/depleted Wallen PD paths.csv", row.names = FALSE)


# identifying shared depleted pathways 
common_depleted_paths <- merge(ibd_depleted, pd_depleted, by = "taxon")   # 26
# identifying shared enriched pathways 
common_enriched_paths <- merge(ibd_enriched, pd_enriched, by = "taxon")   # 8 
# pathways they shared but different in enriched/depleted
common_paths <- merge(ibd_sig_pathways, pd_sig_pathways, by = "taxon")   # 53-26-8 =(19 paths)

write.csv(common_depleted_paths, file = "HMP2_Payami/ANCOMBC2/paths depleted in HMP2 and Wallen.csv", row.names = FALSE)
write.csv(common_enriched_paths, file = "HMP2_Payami/ANCOMBC2/paths enriched in HMP2 and Wallen.csv", row.names = FALSE)


# remove unmapped or unintegrated 
ibd_depleted <- ibd_depleted %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED", taxon))
ibd_enriched <- ibd_enriched %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED", taxon))
pd_depleted <- pd_depleted %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED", taxon))
pd_enriched <- pd_enriched %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED", taxon))

# creating venn diagrams 
library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)

# Extract enriched species names 
enriched_species_IBD <- unique(ibd_enriched$taxon)
enriched_species_PD <- unique(pd_enriched$taxon)

# Extract depleted species names 
depleted_species_IBD <- unique(ibd_depleted$taxon)
depleted_species_PD <- unique(pd_depleted$taxon)


# using ggVennDiagram 
x <- list(
  IBD = depleted_species_IBD,
  PD = depleted_species_PD
)

ggVennDiagram(x[1:3], label_alpha = 0)

ggVennDiagram(
  x[1:3],
  category.names = names(x),
  show_intersect = FALSE,
  set_color = "black",
  set_size = 10,
  label = "both",    # or can use just "count" 
  label_alpha = 0,
  label_geom = "label",
  label_color = "black",
  label_size = 8,  # Adjust label size
  label_percent_digit = 1,  # Adjust decimal places for percentages
  label_txtWidth = 25,  # Adjust text width for labels
  edge_lty = "solid",  # Change edge line type
  edge_size = 1.5,  # Change edge size
  fill = custom_colors,
  title = "Depleted Pathways" # Add a title
) + scale_fill_gradient(low = "lightpink", high = "red")

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 Payami Venn Diagram Pathways Depleted.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")


# using ggVennDiagram 
x <- list(
  IBD = enriched_species_IBD,
  PD = enriched_species_PD
)

ggVennDiagram(x[1:3], label_alpha = 0)

ggVennDiagram(
  x[1:3],
  category.names = names(x),
  show_intersect = FALSE,
  set_color = "black",
  set_size = 10,
  label = "both",    # or can use just "count" 
  label_alpha = 0,
  label_geom = "label",
  label_color = "black",
  label_size = 8,  
  label_percent_digit = 1,  # Adjust decimal places for percentages
  label_txtWidth = 25,  # Adjust text width for labels
  edge_lty = "solid",  # Change edge line type
  edge_size = 1.5,  # Change edge size
  fill = custom_colors,
  title = "Enriched Pathways" # Add a title
) + scale_fill_gradient(low = "lightblue", high = "dodgerblue")

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 Payami Venn Diagram Pathways Enriched.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")
