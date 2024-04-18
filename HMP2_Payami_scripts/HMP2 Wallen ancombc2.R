library(ANCOMBC)
# ran in hipergator 

# Wallen PD 
# species level 
phyloseq_object_s <- readRDS("HMP2_Payami/Phyloseq Objects/Wallen PD species phyloseq object.rds")

ancom_PD <- ancombc2(phyloseq_object_s, 
                     fix_formula = "Case_status + collection_method + total_sequences",
                     tax_level = "Species",
                     p_adj_method = "BH",
                     group = "Case_status",        
                     lib_cut=0,
                     struc_zero=FALSE,
                     neg_lb=FALSE,
                     alpha=0.05,
                     global = FALSE, 
                     pairwise = FALSE)

saveRDS(ancom_PD, "HMP2_Payami/Wallen PD species ancombc2 output.rds")


# genus level 
phyloseq_object_g <- readRDS("HMP2_Payami/Phyloseq Objects/Wallen PD genus phyloseq object.rds")

ancom_PD <- ancombc2(phyloseq_object_g, 
                     fix_formula = "Case_status + collection_method + total_sequences",
                     tax_level = "Genus",
                     p_adj_method = "BH",
                     group = "Case_status",        
                     lib_cut=0,
                     struc_zero=FALSE,
                     neg_lb=FALSE,
                     alpha=0.05,
                     global = FALSE, 
                     pairwise = FALSE)

saveRDS(ancom_PD, "HMP2_Payami/Wallen PD genus ancombc2 output.rds")


# MetaCyc Pathways 
path_phyloseq_object <- readRDS("HMP2_Payami/Phyloseq Objects/Wallen PD pathways phyloseq object.rds")

Payami_pathways <- ancombc2(path_phyloseq_object, 
                            fix_formula = "Case_status + total_sequences + collection_method",
                            tax_level = "Genus",
                            p_adj_method = "BH",
                            group = "Case_status",        
                            lib_cut=0,
                            struc_zero=FALSE,
                            neg_lb=FALSE,
                            alpha=0.05,
                            global = FALSE, 
                            pairwise = FALSE)

saveRDS(Payami_pathways, "HMP2_Payami/ANCOMBC2/Wallen PD ancombc2 metacyc pathways.rds")



# -------------------------------------------------------------------------
# HMP2 IBD  
# species level 
IBD_object <- readRDS("HMP2_Payami/Phyloseq Objects/HMP2 IBD species phyloseq object.rds")

IBD_object@sam_data$diagnosis2 <- as.factor(IBD_object@sam_data$diagnosis2)
IBD_object@sam_data$diagnosis2 <- relevel(IBD_object@sam_data$diagnosis2, ref = "nonIBD")

ancom1 <- ancombc2(IBD_object, 
                   fix_formula = "diagnosis2 + reads_filtered",
                   tax_level = "Species",
                   p_adj_method = "BH",
                   group = "diagnosis2",        
                   lib_cut=0,
                   struc_zero=FALSE,
                   neg_lb=FALSE,
                   alpha=0.05,
                   global = FALSE, 
                   pairwise = FALSE)

saveRDS(ancom1, "HMP2_Payami/ANCOMBC2/HMP2 IBD species ancombc2 output.rds")


# genus level 
IBD_object <- readRDS("HMP2_Payami/Phyloseq Objects/HMP2 IBD genus phyloseq object.rds")

IBD_object@sam_data$diagnosis2 <- as.factor(IBD_object@sam_data$diagnosis2)
IBD_object@sam_data$diagnosis2 <- relevel(IBD_object@sam_data$diagnosis2, ref = "nonIBD")

ancom2 <- ancombc2(IBD_object, 
                   fix_formula = "diagnosis2 + reads_filtered",
                   tax_level = "Genus",
                   p_adj_method = "BH",
                   group = "diagnosis2",        
                   lib_cut=0,
                   struc_zero=FALSE,
                   neg_lb=FALSE,
                   alpha=0.05,
                   global = FALSE, 
                   pairwise = FALSE)

saveRDS(ancom2, "HMP2_Payami/ANCOMBC2/HMP2 IBD genus ancombc2 output.rds")


# MetaCyc Pathways 
path_phyloseq_object <- readRDS("HMP2_Payami/Phyloseq Objects/HMP2 IBD pathways phyloseq object.rds")

path_phyloseq_object@sam_data$diagnosis2 <- as.factor(path_phyloseq_object@sam_data$diagnosis2)
path_phyloseq_object@sam_data$diagnosis2 <- relevel(path_phyloseq_object@sam_data$diagnosis2, ref = "nonIBD")

HMP2_pathways <- ancombc2(path_phyloseq_object, 
                          fix_formula = "diagnosis2 + reads_filtered",
                          tax_level = "Genus",
                          p_adj_method = "BH",
                          group = "diagnosis2",        
                          lib_cut=0,
                          struc_zero=FALSE,
                          neg_lb=FALSE,
                          alpha=0.05,
                          global = FALSE, 
                          pairwise = FALSE)

saveRDS(HMP2_pathways, "HMP2_Payami/ANCOMBC2/HMP2 IBD pathway ancombc2 output.rds")

