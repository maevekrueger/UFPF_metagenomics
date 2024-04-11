library(ggplot2)
library(ggVennDiagram)

# enriched and depleted IBD and PD pathways 
ibd_depleted <- readRDS("UFPF/ANCOMBC2/ANCOMBC2 Tables/Enriched vs Depleted/Pathways/depleted pathways IBD.rds") 
ibd_enriched <- readRDS("UFPF/ANCOMBC2/ANCOMBC2 Tables/Enriched vs Depleted/Pathways/enriched pathways IBD.rds")
pd_depleted <- readRDS("UFPF/ANCOMBC2/ANCOMBC2 Tables/Enriched vs Depleted/Pathways/depleted pathways PD.rds")
pd_enriched <- readRDS("UFPF/ANCOMBC2/ANCOMBC2 Tables/Enriched vs Depleted/Pathways/enriched pathways PD.rds")

# Extract enriched names 
enriched_paths_IBD <- unique(ibd_enriched$taxon)
enriched_paths_PD <- unique(pd_enriched$taxon)

# Extract depleted names 
depleted_paths_IBD <- unique(ibd_depleted$taxon)
depleted_paths_PD <- unique(pd_depleted$taxon)

# DEPLETED
x <- list(
  "UFPF IBD" = depleted_paths_IBD,
  "UFPF PD" = depleted_paths_PD
)

venn <- Venn(x)
data <- process_data(venn)


# CUSTOMIZE GROUP COLORS
# vector for colors - optinos- pick the one thats the least hideous 
colorGroups <- c("UFPF IBD"= "green", "UFPF PD" = "blue")  

# use colorRampPalette to create function that interpolates colors 
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

ggsave("UFPF/Figures/Venn Pathways Depleted.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")


# ----------------------------------------------------------
# ENRICHED
x <- list(
  "UFPF IBD" = enriched_paths_IBD,
  "UFPF PD" = enriched_paths_PD
)

venn <- Venn(x)
data <- process_data(venn)

# CUSTOMIZE GROUP COLORS
# vector for colors - optinos- pick the one thats the least hideous 
colorGroups <- c("UFPF IBD"= "green", "UFPF PD" = "blue")  

# use colorRampPalette to create function that interpolates colors 
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

ggsave("UFPF/Figures/Venn Pathways Enriched.png", dpi = 600, width = 10, height = 8, units = "in", bg = "white")


