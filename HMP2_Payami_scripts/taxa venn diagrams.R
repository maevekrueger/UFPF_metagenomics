library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)

# ----------------------------------------------------------------------------------------------
# plotting with HMP2 IBD combined 
# HMP2 IBD data 
ancom_IBD <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD species ancombc2 output.rds")

res_prim_IBD = ancom_IBD$res
sig_taxa_IBD <- res_prim_IBD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))  # 63 sig species 

sig_taxa_IBD <- sig_taxa_IBD[, c(1, 3, 19, 23)]


# Payami PD data 
ancom_PD1 <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD species ancombc2 output.rds")

res_prim_PD1 = ancom_PD1$res
sig_taxa_PD1 <- res_prim_PD1 %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))

sig_taxa_PD <- sig_taxa_PD1[, c(1, 3, 19, 23)]


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
  title = "Depleted Species" # Add a title
) + scale_fill_gradient(low = "lightpink", high = "red")

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 Payami Venn Diagram Depleted.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")

# Extract enriched species names 
enriched_species_IBD <- unique(enriched_IBD$taxon)
enriched_species_PD <- unique(enriched_PD$taxon)


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
  title = "Enriched Species" # Add a title
) + scale_fill_gradient(low = "lightblue", high = "dodgerblue")

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 Payami Venn Diagram Enriched.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")


# ----------------------------------------------------------------------------------
# genus and species level info combined 
depleted_IBD <- readRDS("HMP2_Payami/ANCOMBC2/depleted genus species HMP2 IBD.rds") 
enriched_IBD <- readRDS("HMP2_Payami/ANCOMBC2/enriched genus species HMP2 IBD.rds")
depleted_PD <- readRDS("HMP2_Payami/ANCOMBC2/depleted genus species Wallen PD.rds") 
enriched_PD <- readRDS("HMP2_Payami/ANCOMBC2/enriched genus species Wallen PD.rds")

# Extract depleted names 
depleted_IBD_names <- unique(depleted_IBD$taxon)
depleted_PD_names <- unique(depleted_PD$taxon)

# using ggVennDiagram 
x <- list(
  IBD = depleted_IBD_names,
  PD = depleted_PD_names
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
  title = "Depleted Taxa" # Add a title
) + scale_fill_gradient(low = "lightpink", high = "red")

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 Payami Venn Diagram Depleted Species Genus.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")

# Extract enriched names 
enriched_IBD_names <- unique(enriched_IBD$taxon)
enriched_PD_names <- unique(enriched_PD$taxon)

# using ggVennDiagram 
x <- list(
  IBD = enriched_IBD_names,
  PD = enriched_PD_names
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
  title = "Enriched Taxa" # Add a title
) + scale_fill_gradient(low = "lightblue", high = "dodgerblue")

ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 Payami Venn Diagram Enriched Species Genus.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")
