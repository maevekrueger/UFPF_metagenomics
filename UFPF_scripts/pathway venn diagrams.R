# enriched and depleted IBD and PD pathways 
# creating venn diagrams 
library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)

ibd_depleted <- readRDS("UFPF/ANCOMBC2/ANCOMBC2 Tables/Enriched vs Depleted/depleted pathways IBD.rds") 
ibd_enriched <- readRDS("UFPF/ANCOMBC2/ANCOMBC2 Tables/Enriched vs Depleted/enriched pathways IBD.rds")
pd_depleted <- readRDS("UFPF/ANCOMBC2/ANCOMBC2 Tables/Enriched vs Depleted/depleted pathways PD.rds")
pd_enriched <- readRDS("UFPF/ANCOMBC2/ANCOMBC2 Tables/Enriched vs Depleted/enriched pathways PD.rds")

# Extract enriched names 
enriched_species_IBD <- unique(ibd_enriched$taxon)
enriched_species_PD <- unique(pd_enriched$taxon)

# Extract depleted names 
depleted_species_IBD <- unique(ibd_depleted$taxon)
depleted_species_PD <- unique(pd_depleted$taxon)


# using ggVennDiagram 
x <- list(
  "UFPF IBD" = depleted_species_IBD,
  "UFPF PD" = depleted_species_PD
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

ggsave("UFPF/Figures/Venn Diagram Pathways Depleted.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")


# using ggVennDiagram 
x <- list(
  "UFPF IBD" = enriched_species_IBD,
  "UFPF PD" = enriched_species_PD
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

ggsave("UFPF/Figures/Venn Diagram Pathways Enriched.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")




# creating a nice tables 
# IBD DEPLETED
table <- ibd_depleted
table <- table %>%
  rename(pathway = taxon)

# using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF IBD Depleted Pathways") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

save_kable(table_formatted, "UFPF/ANCOMBC2/table IBD depleted paths.html")
write.xlsx(as.data.frame(table), "UFPF/ANCOMBC2/table IBD depleted paths.xlsx")

#IBD ENRICHED
table <- ibd_enriched
table <- table %>%
  rename(pathway = taxon)

# Create nice table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF IBD Enriched Pathways") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

save_kable(table_formatted, "UFPF/ANCOMBC2/table IBD enriched paths.html")
write.xlsx(as.data.frame(table), "UFPF/ANCOMBC2/table IBD enriched paths.xlsx")

# PD DEPLETED 
table <- pd_depleted
table <- table %>%
  rename(pathway = taxon)

# Create nice table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF PD Depleted Pathways") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

save_kable(table_formatted, "UFPF/ANCOMBC2/table PD depleted paths.html")
write.xlsx(as.data.frame(table), "UFPF/ANCOMBC2/table PD depleted paths.xlsx")

