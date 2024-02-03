library(tidyverse)
library(WGCNA)
library(afex)
library(emmeans)
library(multcomp)
library(ggplot2)
library(ggbeeswarm)
library(paletteer)
library(camcorder)
library(ggpubr)

# Pathways 
#Loading RDS data ----
hmp2.meta <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")
uf.meta <- readRDS("UFPF/Metadata.rds")
wallen.meta <- readRDS("HMP2_Payami/Payami PD Metadata.rds")

hmp2.stats <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD combined pathway ancombc2 output.rds")
uf.stats <- readRDS("UFPF/ANCOMBC2/ancombc2 pathways output.rds")
wallen.stats <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD ancombc2 metacyc pathways.rds")

hmp2.counts <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD pathway bias corrected abund.rds")
uf.counts <- readRDS("UFPF/ANCOMBC2/bias corrected abund pathways.rds")
wallen.counts <- readRDS("HMP2_Payami/ANCOMBC2/Wallen paths bias corrected abund ancombc2.rds")

#Pull out taxa that are significantly associated with different diseases from external datasets according to lfc and pval ----
hmp2.stats$res %>%
  filter(diff_diagnosis2IBD == T) %>% 
  filter(lfc_diagnosis2IBD > 0) %>%
  pull(taxon) -> hmp2.features.ibd.up

hmp2.stats$res %>%
  filter(diff_diagnosis2IBD == T) %>% 
  filter(lfc_diagnosis2IBD < 0) %>%
  pull(taxon) -> hmp2.features.ibd.down

wallen.stats$res %>%
  filter(diff_Case_statusPD == T) %>%
  filter(lfc_Case_statusPD > 0) %>%
  pull(taxon) -> wallen.features.up

wallen.stats$res %>%
  filter(diff_Case_statusPD == T) %>%
  filter(lfc_Case_statusPD < 0) %>%
  pull(taxon) -> wallen.features.down

# combining taxa between UC/CD in the enriched and depleted modules 
#hmp2.features.down <- union(hmp2.features.ibd.down)
#hmp2.features.up <- union(hmp2.features.ibd.up)


#This function generates the average abundance dataframe in a dataset of interest (counts_of_interest)
#based on over- and under-abundant taxonomic features (up_features and down_features, respectively)
#Those feature lists were generated in the previous section ---
make_coabundance_module <- function(counts_of_interest, up_features, down_features){
  l <- ncol(counts_of_interest)
  
  modules <- vector(length = l)
  
  taxa <- colnames(counts_of_interest)
  
  up_idx <- taxa %in% up_features
  down_idx <- taxa %in% down_features
  
  modules[up_idx] <- "up"
  modules[down_idx] <- "down"
  modules[modules == "FALSE"] <- "ns"
  
  MEs <- moduleEigengenes(counts_of_interest,
                          modules)
  
  data <- MEs$averageExpr
  colnames(data) <- str_remove_all(colnames(data), "AE")
  
  return(data)
}

#Get module scores of PD, UC, and CD in UF data ----
uf_pd_module <- make_coabundance_module(uf.counts,
                                        wallen.features.up,
                                        wallen.features.down)
uf_ibd_module <- make_coabundance_module(uf.counts,
                                         hmp2.features.ibd.up,
                                         hmp2.features.ibd.down)


#Get module scores in Payami's data of PD, UC, and CD associated taxa ----
#In this case, we are both creating and testing the PD module in the PD dataset
#The results from this should be obvious - this is more of a quality control step on the coding part.
#i.e., the module made from species enriched in the PD microbiome should definitely be UP in the PD samples,
#and if we don't get this result then something is wrong with the code I've written
wallen_pd_module <- make_coabundance_module(wallen.counts,
                                            wallen.features.up,
                                            wallen.features.down)
wallen_ibd_module <- make_coabundance_module(wallen.counts,
                                             hmp2.features.ibd.up,
                                             hmp2.features.ibd.down)


#Get module scores in the HMP2 data of PD, UC, and CD associated taxa ----
#As above, this is a bit of a circular analysis with the UC and CD associated taxa, 
#but it lets me make sure that the code is running right 
hmp2_pd_module <- make_coabundance_module(hmp2.counts,
                                          wallen.features.up,
                                          wallen.features.down)
hmp2_ibd_module <- make_coabundance_module(hmp2.counts,
                                           hmp2.features.ibd.up,
                                           hmp2.features.ibd.down)


#These functions runs stats on modules of interest in the data of interest and draws the plot ----
#The `expr` variable is the module dataframe generated above using the make_coabundance_module function
#The `module` variable is a string of "up" or "down," referring to what features you want to see plotted
#Basically, pick the external features you want to see in the UF data and whether those features are up or down originally
#Using uf_pd_module means you're looking at the PD-associated taxa from Wallen et al
#Then, specifying "up" means you want to know about the abundance of the PD-enriched taxa in the UF data
#The logic is the same for the other functions I've created below, just focusing on a different test dataset 
uf_module_plot <- function(expr, module){
  data <- expr %>%
    rownames_to_column(var = "id") %>%
    left_join(uf.meta %>%
                rownames_to_column(var = "id"),
              by = "id")
  
  lm <- aov_ez(id = "id", #aov_ez from the afex package
               dv = module, #can include covariates here if you wanted - can discuss if needed
               between = "Diagnosis2",
               data = data)
  
  print(lm) #this will print the ANOVA table in the console so you can see p, F values for yourself
  
  emm <- emmeans(lm, specs = pairwise ~ Diagnosis2) #pairwise comparisons using the emmeans package
  
  cld <- multcomp::cld(emm, Letters = letters) #LETTERS!!!!!!! from the multcomp package
  cld <- cld %>%
    mutate(.group = str_remove_all(.group, " "))
  
  module2 <- ensym(module) #don't ask me about this... idk why it works
  #I just needed a way to use the direction specified by the `module` variable to also call a column in ggplot's aes(),
  #which doesn't work with character strings. This fixed it for some reason
  
  data %>% 
    left_join(cld, by = "Diagnosis2") %>%
    ggplot(aes(x = Diagnosis2, y = !!module2)) + #again, don't ask me about the !! ... s/o to Stack Overflow 
    geom_quasirandom(aes(fill = Diagnosis2, shape = Diagnosis2),
                     size = 6,
                     show.legend = F,
                     alpha = 0.7) +
    stat_summary(color = "black",
                 fun = "mean",
                 geom = "crossbar",
                 show.legend = F) +
    stat_summary(color = "black",
                 fun.data = mean_se,
                 geom = "errorbar",
                 show.legend = F,
                 width = 0.6,
                 linewidth = 1.2) +
    geom_text(aes(color = Diagnosis2, label = .group,
                  y = max(!!module2 + 0.3 * sd(!!module2, na.rm = T))),
              size = 8,
              show.legend = F) + 
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_paletteer_d(palette = "tvthemes::kimPossible") + #change colors and shapes however you'd like 
    scale_color_paletteer_d(palette = "tvthemes::kimPossible") +
    labs(y = "Average z-scored bias-corrected abundaces") +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_blank())
  
}
gg_record(device= "png",
          units = "in",
          dpi = 600,
          height = 6,
          width = 6) 

uf_module_plot(uf_pd_module, "down") #Try this guy out 
#Ask for a different feature dataframe and different module and see what you get 

gg_stop_recording()

wallen_feature_plot <- function(expr, module){
  data <- expr %>%
    rownames_to_column(var = "id") %>%
    left_join(wallen.meta %>%
                rownames_to_column(var = "id"),
              by = "id")
  
  lm <- aov_ez(id = "id", 
               dv = module, 
               between = "Case_status",
               data = data)
  
  emm <- emmeans(lm, specs = pairwise ~ Case_status) 
  
  cld <- multcomp::cld(emm, Letters = letters) 
  cld <- cld %>%
    mutate(.group = str_remove_all(.group, " "))
  
  module2 <- ensym(module) 
  
  data %>% 
    left_join(cld, by = "Case_status") %>%
    ggplot(aes(x = Case_status, y = !!module2)) + 
    geom_quasirandom(aes(fill = Case_status, shape = Case_status),
                     size = 4,
                     show.legend = F,
                     alpha = 0.7) +
    stat_summary(color = "black",
                 fun = "mean",
                 geom = "crossbar",
                 show.legend = F) +
    stat_summary(color = "black",
                 fun.data = mean_se,
                 geom = "errorbar",
                 show.legend = F,
                 width = 0.6,
                 linewidth = 1.2) +
    geom_text(aes(color = Case_status, label = .group,
                  y = max(!!module2 + 0.3 * sd(!!module2, na.rm = T))),
              size = 8,
              show.legend = F) + 
    scale_shape_manual(values = c(21, 22)) +
    scale_fill_paletteer_d(palette = "ggthemes::stata_s1color") + #change colors and shapes however you'd like 
    scale_color_paletteer_d(palette = "ggthemes::stata_s1color") +
    labs(y = "Average z-scored bias-corrected abundaces") +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_blank())
}
gg_record(device= "png",
          units = "in",
          dpi = 600,
          height = 6,
          width = 4.5) 

wallen_feature_plot(wallen_pd_module, "up")

gg_stop_recording()

hmp2_feature_plot <- function(expr, module){
  data <- expr %>%
    rownames_to_column(var = "Sample") %>%
    left_join(hmp2.meta,
              by = "Sample")
  
  lm <- aov_ez(id = "Sample", 
               dv = module, 
               between = "diagnosis2",
               data = data)
  
  print(lm)
  
  emm <- emmeans(lm, specs = pairwise ~ diagnosis2) 
  
  cld <- multcomp::cld(emm, Letters = letters)
  cld <- cld %>%
    mutate(.group = str_remove_all(.group, " "))
  
  module2 <- ensym(module)
  
  data %>% 
    left_join(cld, by = "diagnosis2") %>%
    mutate(diagnosis2 = factor(diagnosis2,
                               levels = c("nonIBD", "IBD"))) %>%
    ggplot(aes(x = diagnosis2, y = !!module2)) + 
    geom_quasirandom(aes(fill = diagnosis2, shape = diagnosis2),
                     size = 4,
                     show.legend = F,
                     alpha = 0.7) +
    stat_summary(color = "black",
                 fun = "mean",
                 geom = "crossbar",
                 show.legend = F) +
    stat_summary(color = "black",
                 fun.data = mean_se,
                 geom = "errorbar",
                 show.legend = F,
                 width = 0.6,
                 linewidth = 1.2) +
    geom_text(aes(color = diagnosis2, label = .group,
                  y = max(!!module2 + 0.3 * sd(!!module2, na.rm = T))),
              size = 8,
              show.legend = F) + 
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_paletteer_d(palette = "ggthemes::excel_Slice") + #change colors and shapes however you'd like 
    scale_color_paletteer_d(palette = "ggthemes::excel_Slice") +
    labs(y = "Average z-scored bias-corrected abundaces") +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_blank())
}
gg_record(device = "png",
          units = "in", dpi = 600,
          height = 6, width = 6)

hmp2_feature_plot(hmp2_pd_module, "up")

gg_stop_recording()

theme <- theme(plot.title = element_text(hjust = 0.5, size = 20))

annotate_figure(
  ggarrange(uf_module_plot(uf_pd_module, "up") +
              ggtitle("Enriched in Wallen PD") + theme,
            uf_module_plot(uf_pd_module, "down") +
              ggtitle("Depleted in Wallen PD") + theme,
            uf_module_plot(uf_ibd_module, "up") +
              ggtitle("Enriched in HMP2 IBD") + theme,
            uf_module_plot(uf_ibd_module, "down") +
              ggtitle("Depleted in HMP2 IBD") + theme,
            nrow = 1, ncol = 4),
  top = text_grob("Abundance of disease-associated pathways in the UFPF dataset",
                  size = 24, face = "bold")
) -> p1

annotate_figure(
  ggarrange(wallen_feature_plot(wallen_pd_module, "up") +
              ggtitle("Enriched in Wallen PD") + theme,
            wallen_feature_plot(wallen_pd_module, "down") +
              ggtitle("Depleted in Wallen PD") + theme,
            wallen_feature_plot(wallen_ibd_module, "up") +
              ggtitle("Enriched in HMP2 IBD") + theme,
            wallen_feature_plot(wallen_ibd_module, "down") +
              ggtitle("Depleted in HMP2 IBD") + theme,
            nrow = 1, ncol = 4),
  top = text_grob("Abundance of disease-associated pathways in the Wallen dataset",
                  size = 24, face = 'bold')
) -> p2

annotate_figure(
  ggarrange(hmp2_feature_plot(hmp2_pd_module, "up") +
              ggtitle("Enriched in Wallen PD") + theme,
            hmp2_feature_plot(hmp2_pd_module, "down") +
              ggtitle("Depleted in Wallen PD") + theme,
            hmp2_feature_plot(hmp2_ibd_module, "up") +
              ggtitle("Enriched in HMP2 IBD") + theme,
            hmp2_feature_plot(hmp2_ibd_module, "down") +
              ggtitle("Depleted in HMP2 IBD") + theme,
            nrow = 1, ncol = 4),
  top = text_grob("Abundance of disease-associated pathways in the HMP2 dataset",
                  size = 24, face = "bold")
) -> p3

gg_record(device = "png",
          units = "in", dpi = 600,
          height = 15, width = 20)

ggarrange(p1, p2, p3,
          nrow = 3, ncol = 1)
ggsave("figures/modules pathways w second batch.jpg",
       units = "in", dpi = 600,
       height = 15, width = 20)

gg_stop_recording()


# UFPF data only 
gg_record(device = "png",
          units = "in", dpi = 600,
          height = 8, width = 18)

p1
ggsave("figures/UFPF_only_species_modules_IBDcombined_pathways.jpg",
       units = "in", dpi = 600,
       height = 8, width = 18)

gg_stop_recording()