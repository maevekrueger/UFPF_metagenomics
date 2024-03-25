# creating tables of ANCOMBC2 output UFPF data 
library(reshape2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(knitr)
library(kableExtra)

# Species-Level
ancom <- readRDS("UFPF/ANCOMBC2/ancombc2 species.rds")

res_prim = ancom$res
res_pair = ancom$res_pair      # 233 species examined 

table <- res_pair[, c(1:7, 11:16)]

table <- table %>%
  rename(
    species = taxon,
    lfc_IBD = lfc_Diagnosis2IBD,
    lfc_PD = lfc_Diagnosis2PD,
    lfc_PDvsIBD = lfc_Diagnosis2PD_Diagnosis2IBD, 
    se_IBD = se_Diagnosis2IBD,
    se_PD = se_Diagnosis2PD,
    se_PDvsIBD = se_Diagnosis2PD_Diagnosis2IBD,
    pval_IBD = p_Diagnosis2IBD,
    pval_PD = p_Diagnosis2PD,
    pval_PDvsIBD = p_Diagnosis2PD_Diagnosis2IBD,
    adj_pval_IBD = q_Diagnosis2IBD,
    adj_pval_PD = q_Diagnosis2PD,
    adj_pval_PDvsIBD = q_Diagnosis2PD_Diagnosis2IBD
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF ANCOMBC2 Species Taxonomic Output") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

save_kable(table_formatted, "UFPF/ANCOMBC2/table all species ancombc2 output.html")

# Save as Excel
write.xlsx(as.data.frame(table), "UFPF/ANCOMBC2/table all species ancombc2 output.xlsx")


# -----------------------------------------------------------------------------------
# GENUS-LEVEL
ancom <- readRDS("UFPF/ANCOMBC2/ancombc2 genus.rds")

res_prim = ancom$res
res_pair = ancom$res_pair             # 140 genera examined

table <- res_pair[, c(1:7, 11:16)]

table <- table %>%
  rename(
    genus = taxon,
    lfc_IBD = lfc_Diagnosis2IBD,
    lfc_PD = lfc_Diagnosis2PD,
    lfc_PDvsIBD = lfc_Diagnosis2PD_Diagnosis2IBD, 
    se_IBD = se_Diagnosis2IBD,
    se_PD = se_Diagnosis2PD,
    se_PDvsIBD = se_Diagnosis2PD_Diagnosis2IBD,
    pval_IBD = p_Diagnosis2IBD,
    pval_PD = p_Diagnosis2PD,
    pval_PDvsIBD = p_Diagnosis2PD_Diagnosis2IBD,
    adj_pval_IBD = q_Diagnosis2IBD,
    adj_pval_PD = q_Diagnosis2PD,
    adj_pval_PDvsIBD = q_Diagnosis2PD_Diagnosis2IBD
  )

table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF ANCOMBC2 Genus Taxonomic Output") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

save_kable(table_formatted, "UFPF/ANCOMBC2/ANCOMBC2 tables/table genus taxonomic output.html")

# Save as Excel
write.xlsx(as.data.frame(table), "UFPF/ANCOMBC2/ANCOMBC2 tables/table genus taxonomic output.xlsx")



# -------------------------------------------------------------------------------
# Pathway output
ancom <- readRDS("UFPF/ANCOMBC2/ancombc2 pathways output.rds")

res_prim = ancom$res
res_pair = ancom$res_pair             # 350 pathways examined

sig_paths <- res_pair %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_"))))    # 132 sig pathways

# create a table with only significant pathways ***
table <- sig_paths[, c(1:7, 11:16)]

table <- table %>%
  rename(
    pathway = taxon,
    lfc_IBD = lfc_Diagnosis2IBD,
    lfc_PD = lfc_Diagnosis2PD,
    lfc_PDvsIBD = lfc_Diagnosis2PD_Diagnosis2IBD, 
    se_IBD = se_Diagnosis2IBD,
    se_PD = se_Diagnosis2PD,
    se_PDvsIBD = se_Diagnosis2PD_Diagnosis2IBD,
    pval_IBD = p_Diagnosis2IBD,
    pval_PD = p_Diagnosis2PD,
    pval_PDvsIBD = p_Diagnosis2PD_Diagnosis2IBD,
    adj_pval_IBD = q_Diagnosis2IBD,
    adj_pval_PD = q_Diagnosis2PD,
    adj_pval_PDvsIBD = q_Diagnosis2PD_Diagnosis2IBD
  )

table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF ANCOMBC2 Significant Pathways") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

save_kable(table_formatted, "UFPF/ANCOMBC2/ANCOMBC2 tables/table pathway output.html")

# Save as Excel
write.xlsx(as.data.frame(table), "UFPF/ANCOMBC2/ANCOMBC2 tables/table pathway output.xlsx")



# separate PD and IBD sig paths 
sig_paths_PD <- res_pair[res_pair$diff_Diagnosis2PD == TRUE, ]     # 1 sig 
sig_paths_IBD <- res_pair[res_pair$diff_Diagnosis2IBD == TRUE, ]   # 93 sig 

# create a table with only significant pathways ***
table1 <- sig_paths_PD[, c(1:7, 11:16)]

table1 <- table1 %>%
  rename(
    pathway = taxon,
    lfc_IBD = lfc_Diagnosis2IBD,
    lfc_PD = lfc_Diagnosis2PD,
    lfc_PDvsIBD = lfc_Diagnosis2PD_Diagnosis2IBD, 
    se_IBD = se_Diagnosis2IBD,
    se_PD = se_Diagnosis2PD,
    se_PDvsIBD = se_Diagnosis2PD_Diagnosis2IBD,
    pval_IBD = p_Diagnosis2IBD,
    pval_PD = p_Diagnosis2PD,
    pval_PDvsIBD = p_Diagnosis2PD_Diagnosis2IBD,
    adj_pval_IBD = q_Diagnosis2IBD,
    adj_pval_PD = q_Diagnosis2PD,
    adj_pval_PDvsIBD = q_Diagnosis2PD_Diagnosis2IBD
  )

table2 <- sig_paths_IBD[, c(1:7, 11:16)]

table2 <- table2 %>%
  rename(
    pathway = taxon,
    lfc_IBD = lfc_Diagnosis2IBD,
    lfc_PD = lfc_Diagnosis2PD,
    lfc_PDvsIBD = lfc_Diagnosis2PD_Diagnosis2IBD, 
    se_IBD = se_Diagnosis2IBD,
    se_PD = se_Diagnosis2PD,
    se_PDvsIBD = se_Diagnosis2PD_Diagnosis2IBD,
    pval_IBD = p_Diagnosis2IBD,
    pval_PD = p_Diagnosis2PD,
    pval_PDvsIBD = p_Diagnosis2PD_Diagnosis2IBD,
    adj_pval_IBD = q_Diagnosis2IBD,
    adj_pval_PD = q_Diagnosis2PD,
    adj_pval_PDvsIBD = q_Diagnosis2PD_Diagnosis2IBD
  )

table_formatted1 <- table1 %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF PD Significant Pathways") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted1)

table_formatted2 <- table2 %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "UFPF IBD Significant Pathways") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted2)

# Save as Excel
write.xlsx(as.data.frame(table1), "UFPF/ANCOMBC2/ANCOMBC2 tables/table PD sig pathways.xlsx")
write.xlsx(as.data.frame(table2), "UFPF/ANCOMBC2/ANCOMBC2 tables/table IBD sig pathways.xlsx")
