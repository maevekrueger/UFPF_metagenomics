library(tidyverse)
library(broom)
library(openxlsx)

# exmaining the influence of potential confounding variables 

# HMP2 IBD data 
IBD_metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")
IBD_metadata <- column_to_rownames(IBD_metadata, var = "Sample")

# removing some columns I don't need
columns_to_remove <- c(14:19)
IBD_metadata <- IBD_metadata[, -columns_to_remove]

IBD_metadata <- IBD_metadata %>%
  rename(Age = consent_age, Diagnosis = diagnosis, Sex = sex, Reads = reads_filtered,
         Immunosuppressants = Immunosuppressants..e.g..oral.corticosteroids., 
         Diagnosis_Age = Age.at.diagnosis, Site = site_name, 
         Colonoscopy_past_2_weeks = X2..In.the.past.2.weeks..have.you.undergone.a.colonoscopy.or.other.procedure,
         Diarrhea_past_2_weeks = X4..In.the.past.2.weeks..have.you.had.diarrhea., 
         Bowel_surgery = X6..Have.you.ever.had.bowel.surgery.)

# label probiotic column N or Y for probiotics taken within 7 days of collection
IBD_metadata <- IBD_metadata %>%
  mutate(Probiotic = case_when(
    grepl("^No", Probiotic) ~ "N",
    Probiotic == "" ~ NA_character_,
    TRUE ~ "Y"
  ))

order <- c("nonIBD", "IBD", "Overall")
IBD_metadata$diagnosis2 <- factor(IBD_metadata$diagnosis2, levels = order)

IBD_metadata <- IBD_metadata %>%
  mutate(Age = as.numeric(Age)) 

IBD_metadata$Reads <- scale(IBD_metadata$Reads)
IBD_metadata$Age <- scale(IBD_metadata$Age)


# specify response variables 
ancom_abundances <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD species bias corrected abund.rds")

species <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD species ancombc2 output.rds")
res_prim_IBD = species$res        
sig_taxa_IBD <- res_prim_IBD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))        # 64 significant
species_names <- sig_taxa_IBD$taxon
filtered_ancom_abundances <- ancom_abundances[, colnames(ancom_abundances) %in% species_names]

# pull out predictor variables according to the demographics/metadata table stats 
selected_columns <- IBD_metadata %>%
  select(Sex, Age, diagnosis2, Reads, Site, Probiotic, Immunosuppressants, 
         Antibiotics)

# Reorder selected_columns to match the order of row names 
reordered <- selected_columns[match(rownames(filtered_ancom_abundances), rownames(selected_columns)), ]


combined <- cbind(filtered_ancom_abundances, reordered)

# Remove rows with missing metadata values
combined <- combined[complete.cases(combined[, 65:72]), ]

combined$Sex <- factor(combined$Sex)
combined$diagnosis2 <- factor(combined$diagnosis2)
combined$Site <- factor(combined$Site)

convert_values <- function(x) {
  x <- ifelse(x %in% c("No", "N", "never"), 0,
              ifelse(x %in% c("Yes", "Y", "yes"), 1, NA))
  return(x)
}
combined_HMP2 <- combined %>%
  mutate_at(vars(70:72), list(~convert_values(.)))


# Fit a multivariate regression model
model <- lm(cbind(Bifidobacterium_adolescentis, Bacteroides_eggerthii, Bacteroides_fragilis, 
                  Bacteroides_massiliensis, Bacteroides_nordii, Bacteroides_stercoris, Bacteroides_thetaiotaomicron, 
                  Bacteroidales_bacterium_ph8, Barnesiella_intestinihominis, Coprobacter_fastidiosus, 
                  Odoribacter_laneus, Odoribacter_unclassified, Parabacteroides_distasonis, 
                  Parabacteroides_goldsteinii, Parabacteroides_unclassified, Prevotella_copri, 
                  Alistipes_indistinctus, Alistipes_onderdonkii, Alistipes_putredinis, 
                  Alistipes_senegalensis, Alistipes_shahii, Clostridium_asparagiforme, 
                  Clostridium_bolteae, Clostridium_clostridioforme, Clostridium_hathewayi, 
                  Clostridium_symbiosum, Clostridiales_bacterium_1_7_47FAA, Flavonifractor_plautii, 
                  Eubacterium_eligens, Eubacterium_hallii, Eubacterium_ramulus, Eubacterium_rectale, 
                  Eubacterium_ventriosum, Anaerostipes_hadrus, Ruminococcus_gnavus, Ruminococcus_obeum, 
                  Butyrivibrio_crossotus, Coprococcus_catus, Coprococcus_comes, Coprococcus_sp_ART55_1,
                  Dorea_formicigenerans, Dorea_unclassified, Lachnospiraceae_bacterium_3_1_46FAA, 
                  Lachnospiraceae_bacterium_5_1_63FAA, Roseburia_hominis, Roseburia_intestinalis, 
                  Roseburia_inulinivorans, Peptostreptococcaceae_noname_unclassified, 
                  Faecalibacterium_prausnitzii, Ruminococcus_bromii, Ruminococcus_callidus, 
                  Subdoligranulum_unclassified, Coprobacillus_unclassified, Eubacterium_biforme, 
                  Veillonella_atypica, Veillonella_dispar, Veillonella_parvula, Veillonella_unclassified, 
                  Oxalobacter_formigenes, Sutterella_wadsworthensis, Bilophila_unclassified, 
                  Escherichia_coli, Haemophilus_parainfluenzae, Akkermansia_muciniphila) ~ ., data = combined_HMP2)

# Print the summary of the model
summary(model)

saveRDS(summary(model), "Regression scripts/HMP2_model_summary.rds")

# create a table of regression stats 
HMP2_tidy <- tidy(model)
write.xlsx(HMP2_tidy, "Regression scripts/HMP2 IBD regression summary fixed.xlsx")

# adjust the p values for multiple comparisons w/ Benjamini-Hochberg correction
p_values <- HMP2_tidy[, 6]
p_values <- unlist(p_values)
adjusted_p_values <- p.adjust(p_values, method = "BH")
HMP2_tidy$adjusted_p_values <- adjusted_p_values

# Save the updated dataframe to Excel
write.xlsx(HMP2_tidy, "Regression scripts/HMP2 IBD regression summary fixed.xlsx")


# Filtering to see which species retained association with IBD 
# Filter where "term" is diagnosis2IBD and "adjusted_p_values" is <= 0.05
filtered_table <- HMP2_tidy %>%
  filter(term == "diagnosis2IBD" & adjusted_p_values <= 0.05)

# response column (species names)
filtered_species <- filtered_table$response
print(filtered_species)   # 35 species




# ----------------------------------------------------------------------------
# without confounders 
# HMP2 IBD data 
IBD_metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")
IBD_metadata <- column_to_rownames(IBD_metadata, var = "Sample")

# removing some columns I don't need
columns_to_remove <- c(14:19)
IBD_metadata <- IBD_metadata[, -columns_to_remove]

IBD_metadata <- IBD_metadata %>%
  rename(Diagnosis = diagnosis, Reads = reads_filtered)

order <- c("nonIBD", "IBD", "Overall")
IBD_metadata$diagnosis2 <- factor(IBD_metadata$diagnosis2, levels = order)

IBD_metadata$Reads <- scale(IBD_metadata$Reads)

# specify response variables 
ancom_abundances <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD species bias corrected abund.rds")

species <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD species ancombc2 output.rds")
res_prim_IBD = species$res        
sig_taxa_IBD <- res_prim_IBD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))        # 64 significant
species_names <- sig_taxa_IBD$taxon
filtered_ancom_abundances <- ancom_abundances[, colnames(ancom_abundances) %in% species_names]

# pull out predictor variables according to the demographics/metadata table stats 
selected_columns <- IBD_metadata %>%
  select(diagnosis2, Reads)

# Reorder selected_columns to match the order of row names 
reordered <- selected_columns[match(rownames(filtered_ancom_abundances), rownames(selected_columns)), ]

combined <- cbind(filtered_ancom_abundances, reordered)

combined$diagnosis2 <- factor(combined$diagnosis2)

# Fit a multivariate regression model
model <- lm(cbind(Bifidobacterium_adolescentis, Bacteroides_eggerthii, Bacteroides_fragilis, 
                  Bacteroides_massiliensis, Bacteroides_nordii, Bacteroides_stercoris, Bacteroides_thetaiotaomicron, 
                  Bacteroidales_bacterium_ph8, Barnesiella_intestinihominis, Coprobacter_fastidiosus, 
                  Odoribacter_laneus, Odoribacter_unclassified, Parabacteroides_distasonis, 
                  Parabacteroides_goldsteinii, Parabacteroides_unclassified, Prevotella_copri, 
                  Alistipes_indistinctus, Alistipes_onderdonkii, Alistipes_putredinis, 
                  Alistipes_senegalensis, Alistipes_shahii, Clostridium_asparagiforme, 
                  Clostridium_bolteae, Clostridium_clostridioforme, Clostridium_hathewayi, 
                  Clostridium_symbiosum, Clostridiales_bacterium_1_7_47FAA, Flavonifractor_plautii, 
                  Eubacterium_eligens, Eubacterium_hallii, Eubacterium_ramulus, Eubacterium_rectale, 
                  Eubacterium_ventriosum, Anaerostipes_hadrus, Ruminococcus_gnavus, Ruminococcus_obeum, 
                  Butyrivibrio_crossotus, Coprococcus_catus, Coprococcus_comes, Coprococcus_sp_ART55_1,
                  Dorea_formicigenerans, Dorea_unclassified, Lachnospiraceae_bacterium_3_1_46FAA, 
                  Lachnospiraceae_bacterium_5_1_63FAA, Roseburia_hominis, Roseburia_intestinalis, 
                  Roseburia_inulinivorans, Peptostreptococcaceae_noname_unclassified, 
                  Faecalibacterium_prausnitzii, Ruminococcus_bromii, Ruminococcus_callidus, 
                  Subdoligranulum_unclassified, Coprobacillus_unclassified, Eubacterium_biforme, 
                  Veillonella_atypica, Veillonella_dispar, Veillonella_parvula, Veillonella_unclassified, 
                  Oxalobacter_formigenes, Sutterella_wadsworthensis, Bilophila_unclassified, 
                  Escherichia_coli, Haemophilus_parainfluenzae, Akkermansia_muciniphila) ~ ., data = combined)

# Print the summary of the model
summary(model)

saveRDS(summary(model), "Regression scripts/HMP2_model_summary_wo_confounders.rds")

# create a table of regression stats 
HMP2_tidy2 <- tidy(model)

# adjust the p values for multiple comparisons w/ Benjamini-Hochberg correction
p_values <- HMP2_tidy2[, 6]
p_values <- unlist(p_values)
adjusted_p_values <- p.adjust(p_values, method = "BH")
HMP2_tidy2$adjusted_p_values <- adjusted_p_values

# Save the updated dataframe to Excel
write.xlsx(HMP2_tidy2, "Regression scripts/HMP2 IBD regression summary wo confounders fixed.xlsx")


# Extract taxa names from HMP2_tidy
taxa_names <- HMP2_tidy$response[HMP2_tidy$term == "diagnosis2IBD"]

# Extract coefficients associated with the term "diagnosis2IBD" from HMP2_tidy
coef_with_confounding <- HMP2_tidy$estimate[HMP2_tidy$term == "diagnosis2IBD"]

# Extract coefficients associated with the term "diagnosis2IBD" from HMP2_tidy2
coef_without_confounding <- HMP2_tidy2$estimate[HMP2_tidy2$term == "diagnosis2IBD"]

# Calculate the difference in coefficients
difference <- coef_without_confounding - coef_with_confounding

# Calculate percent change
percent_change <- (difference / coef_without_confounding) * 100

# Create a new dataframe to store the results
results <- data.frame(
  Taxa = taxa_names,
  Coefficient_Without_Confounding = coef_without_confounding,
  Coefficient_With_Confounding = coef_with_confounding,
  Difference = difference,
  Percent_Change = percent_change
)

# Display the results
print(results)

write.xlsx(results, "Regression scripts/HMP2 IBD difference of confounders fixed.xlsx")

