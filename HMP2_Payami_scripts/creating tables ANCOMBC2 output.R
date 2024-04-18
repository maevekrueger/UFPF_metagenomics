# creating tables of ANCOMBC2 output Wallen + HMP2 data 
library(reshape2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(knitr)
library(kableExtra)
library(tidyverse)

# Species-Level Wallen PD 
ancom_PD <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD species ancombc2 output.rds")

res_prim_PD = ancom_PD$res        # 195 species examined

sig_taxa_PD <- res_prim_PD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))        # 11 significant

# table of just significant output 
table <- sig_taxa_PD[, c(1, 3, 7, 15, 19)]

table <- table %>%
  rename(
    species = taxon,
    lfc_PD = lfc_Case_statusPD,
    se_PD = se_Case_statusPD,
    pval_PD = p_Case_statusPD,
    adj_pval_PD = q_Case_statusPD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "Wallen PD ANCOMBC2 Significant Species") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)


# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table Wallen PD significant species.xlsx")


# table of all output 
table <- res_prim_PD[, c(1, 3, 7, 15, 19)]

table <- table %>%
  rename(
    species = taxon,
    lfc_PD = lfc_Case_statusPD,
    se_PD = se_Case_statusPD,
    pval_PD = p_Case_statusPD,
    adj_pval_PD = q_Case_statusPD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "Wallen PD ANCOMBC2 Species Output") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table Wallen PD species ancombc2 output.xlsx")


# -------------------------
# GENUS-LEVEL Wallen PD 
ancom_PD <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD genus ancombc2 output.rds")

res_prim_PD = ancom_PD$res        # 89 genera examined

sig_taxa_PD <- res_prim_PD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))        # 5 significant

# just significant genera
table <- sig_taxa_PD[, c(1, 3, 7, 15, 19)]

table <- table %>%
  rename(
    genus = taxon,
    lfc_PD = lfc_Case_statusPD,
    se_PD = se_Case_statusPD,
    pval_PD = p_Case_statusPD,
    adj_pval_PD = q_Case_statusPD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "Wallen PD ANCOMBC2 Significant Genera") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table Wallen PD significant genus.xlsx")


# table of all genera 
table <- res_prim_PD[, c(1, 3, 7, 15, 19)]

table <- table %>%
  rename(
    genus = taxon,
    lfc_PD = lfc_Case_statusPD,
    se_PD = se_Case_statusPD,
    pval_PD = p_Case_statusPD,
    adj_pval_PD = q_Case_statusPD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "Wallen PD ANCOMBC2 Genus Output") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table Wallen PD genus ancombc2 output.xlsx")


# ---------------
# Wallen PD PATHWAYS 
ancom_PD <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD ancombc2 metacyc pathways.rds")

res_prim_PD = ancom_PD$res        # 356 examined

sig_path_PD <- res_prim_PD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))        # 198 significant

# only depicting significant pathways
table <- sig_path_PD[, c(1, 3, 7, 15, 19)]

table <- table %>%
  rename(
    pathway = taxon,
    lfc_PD = lfc_Case_statusPD,
    se_PD = se_Case_statusPD,
    pval_PD = p_Case_statusPD,
    adj_pval_PD = q_Case_statusPD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "Wallen PD ANCOMBC2 Significant Pathways") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table Wallen PD pathways ancombc2 output.xlsx")



# --------------------------------------------------------------------------------
# Species-Level HMP2 IBD 
ancom_IBD <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD species ancombc2 output.rds")
# fixed effects: age + sex + diagnosis 
res_prim_IBD = ancom_IBD$res

sig_taxa_IBD <- res_prim_IBD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))    # 64 significant species 

# table of only significant species 
table <- sig_taxa_IBD[, c(1, 3, 6, 12, 15)]

table <- table %>%
  rename(
    species = taxon,
    lfc_IBD = lfc_diagnosis2IBD,
    se_IBD = se_diagnosis2IBD,
    pval_IBD = p_diagnosis2IBD,
    adj_pval_IBD = q_diagnosis2IBD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "HMP2 IBD ANCOMBC2 Significant Species") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

#save_kable(table_formatted, "HMP2_Payami/ANCOMBC2/table HMP2 IBD species ancombc2 output.html")

# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table HMP2 IBD significant species.xlsx")


# table of all species 
table <- res_prim_IBD[, c(1, 3, 6, 12, 15)]

table <- table %>%
  rename(
    species = taxon,
    lfc_IBD = lfc_diagnosis2IBD,
    se_IBD = se_diagnosis2IBD,
    pval_IBD = p_diagnosis2IBD,
    adj_pval_IBD = q_diagnosis2IBD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "HMP2 IBD ANCOMBC2 Species Output") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table HMP2 IBD species ancombc2 output.xlsx")


# -------------------------
# GENUS-LEVEL HMP2 IBD 
ancom_IBD <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD genus ancombc2 output.rds")

res_prim_IBD = ancom_IBD$res        # 52 genera examined

sig_taxa_IBD <- res_prim_IBD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))        # 26 significant

# table of significant genera 
table <- sig_taxa_IBD[, c(1, 3, 6, 12, 15)]

table <- table %>%
  rename(
    genus = taxon,
    lfc_IBD = lfc_diagnosis2IBD,
    se_IBD = se_diagnosis2IBD,
    pval_IBD = p_diagnosis2IBD,
    adj_pval_IBD = q_diagnosis2IBD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "HMP2 IBD ANCOMBC2 Significant Genera") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table HMP2 IBD significant genus.xlsx")


# table of all genera 
table <- res_prim_IBD[, c(1, 3, 6, 12, 15)]

table <- table %>%
  rename(
    genus = taxon,
    lfc_IBD = lfc_diagnosis2IBD,
    se_IBD = se_diagnosis2IBD,
    pval_IBD = p_diagnosis2IBD,
    adj_pval_IBD = q_diagnosis2IBD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "HMP2 IBD ANCOMBC2 Genus Output") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table HMP2 IBD genus ancombc2 output.xlsx")


# ----------------------
# HMP2 IBD PATHWAYS
ancom_IBD <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD pathway ancombc2 output.rds")

res_prim_IBD = ancom_IBD$res        # 294 examined

sig_path_IBD <- res_prim_IBD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))        # 162 significant

# only depicting significant paths
table <- sig_path_IBD[, c(1, 3, 6, 12, 15)]

table <- table %>%
  rename(
    pathway = taxon,
    lfc_IBD = lfc_diagnosis2IBD,
    se_IBD = se_diagnosis2IBD,
    pval_IBD = p_diagnosis2IBD,
    adj_pval_IBD = q_diagnosis2IBD,
  )

# Create table using kable and kableExtra
table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "HMP2 IBD ANCOMBC2 Significant Pathways") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)

# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table HMP2 IBD pathways ancombc2 output.xlsx")



# ------------
# SHARED pathways 
depleted_paths <- readRDS("HMP2_Payami/ANCOMBC2/depleted pathways shared.rds")
enriched_paths <- readRDS("HMP2_Payami/ANCOMBC2/enriched pathways shared.rds")

combined_paths <- rbind(depleted_paths, enriched_paths)

table <- combined_paths[, c(1:5, 7, 8, 10, 11)]

table <- table %>%
  rename(
    pathway = taxon,
    lfc_IBD = lfc_diagnosis2IBD,
    se_IBD = se_diagnosis2IBD,
    pval_IBD = p_diagnosis2IBD,
    adj_pval_IBD = q_diagnosis2IBD,
    lfc_PD = lfc_Case_statusPD,
    se_PD = se_Case_statusPD,
    pval_PD = p_Case_statusPD,
    adj_pval_PD = q_Case_statusPD
  )

# rearrange 
table <- table[, c(1:2, 6, 3, 7, 4, 8, 5, 9:ncol(table))]

table_formatted <- table %>%
  kable(format = "html", table.attr = 'class="table table-striped table-hover"', caption = "Wallen PD + HMP2 IBD Shared Pathways") %>%
  kable_styling("striped", full_width = FALSE)

print(table_formatted)


# Save as Excel
write.xlsx(as.data.frame(table), "HMP2_Payami/ANCOMBC2/ANCOMBC2 tables/table shared pathways ancombc2 output.xlsx")
