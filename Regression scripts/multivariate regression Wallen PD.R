library(tidyverse)
library(broom)
library(openxlsx)

# exmaining the influence of potential confounding variables 

# Wallen PD data 
PD_metadata <- readRDS("HMP2_Payami/Wallen PD Metadata.rds")
PD_metadata <- PD_metadata %>%
  rename(Age = Age_at_collection, Diagnosis = Case_status)
PD_metadata <- PD_metadata %>%
  mutate(Sex = ifelse(Sex == "M", "Male", ifelse(Sex == "F", "Female", Sex)))

PD_metadata$total_sequences <- scale(PD_metadata$total_sequences)
PD_metadata$Age <- scale(PD_metadata$Age)

# for some reason "SC0637" "SC0558" are not present in the bias corrected abundance output
# so I am removing these from the metadata so it doesn't cause problems when we combine dataframes
PD_metadata <- PD_metadata[!rownames(PD_metadata) %in% c("SC0637", "SC0558"), ]

# specify response variables 
ancom_abundances <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD species bias corrected abund.rds")

species <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD species ancombc2 output.rds")
res_prim_PD = species$res        # 195 species examined
sig_taxa_PD <- res_prim_PD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))        # 11 significant
species_names <- sig_taxa_PD$taxon
filtered_ancom_abundances <- ancom_abundances[, colnames(ancom_abundances) %in% species_names]

# pull out predictor variables according to the demographics/metadata table stats 
selected_columns <- PD_metadata %>%
  select(Sex, Age, Diagnosis, total_sequences, Laxatives, Pain_med, Depression_anxiety_mood_med, 
         Birth_control_or_estrogen, Antihistamines, Probiotic, Sleep_aid)


# Reorder selected_columns to match the order of row names in 
reordered <- selected_columns[match(rownames(filtered_ancom_abundances), rownames(selected_columns)), ]

combined <- cbind(filtered_ancom_abundances, reordered)

# Remove rows with missing metadata values
combined <- combined[complete.cases(combined[, 16:22]), ]

combined$Sex <- factor(combined$Sex)
combined$Diagnosis <- factor(combined$Diagnosis)

convert_values <- function(x) {
  x <- ifelse(x %in% c("N", "N ", "never"), 0,
              ifelse(x %in% c("Y", "Y ", "yes"), 1, NA))
  return(x)
}
combined_Wallen <- combined %>%
  mutate_at(vars(16:22), list(~convert_values(.)))


# Fit a multivariate regression model
model <- lm(cbind(Actinomyces_oris, Bifidobacterium_dentium, Streptococcus_mutans,
                  Anaerostipes_hadrus, Blautia_wexlerae, Eisenbergiella_tayi,
                  Roseburia_intestinalis, Clostridium_leptum, Ruminococcaceae_bacterium_D5, 
                  Ruminococcus_lactaris, Escherichia_coli) ~ ., data = combined_Wallen)

# Print the summary of the model
summary(model)

saveRDS(summary(model), "Regression scripts/Wallen_model_summary.rds")

# create a table of regression stats 
Wallen_tidy <- tidy(model)

# adjust the p values for multiple comparisons w/ Benjamini-Hochberg correction
p_values <- Wallen_tidy[, 6]
p_values <- unlist(p_values)
adjusted_p_values <- p.adjust(p_values, method = "BH")
Wallen_tidy$adjusted_p_values <- adjusted_p_values

# Save
write.xlsx(Wallen_tidy, "Regression scripts/Wallen PD regression summary fixed.xlsx")


# Filtering to see which species retained association with PD 
# Filter where "term" is DiagnosisPD and "adjusted_p_values" is <= 0.05
filtered_table <- Wallen_tidy %>%
  filter(term == "DiagnosisPD" & adjusted_p_values <= 0.05)

# response column (species names)
filtered_species <- filtered_table$response
print(filtered_species)   # 9 species





# -------------------------------------------------------------------------
# without confounders 

# Wallen PD data 
PD_metadata <- readRDS("HMP2_Payami/Wallen PD Metadata.rds")
PD_metadata <- PD_metadata %>%
  rename(Diagnosis = Case_status)

PD_metadata$total_sequences <- scale(PD_metadata$total_sequences)

# for some reason "SC0637" "SC0558" are not present in the bias corrected abundance output
# so I am removing these from the metadata so it doesn't cause problems when we combine dataframes
PD_metadata <- PD_metadata[!rownames(PD_metadata) %in% c("SC0637", "SC0558"), ]

# specify response variables 
ancom_abundances <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD species bias corrected abund.rds")

species <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD species ancombc2 output.rds")
res_prim_PD = species$res        # 195 species examined
sig_taxa_PD <- res_prim_PD %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))        # 11 significant
species_names <- sig_taxa_PD$taxon
filtered_ancom_abundances <- ancom_abundances[, colnames(ancom_abundances) %in% species_names]

# pull out predictor variables 
selected_columns <- PD_metadata %>%
  select(Diagnosis, total_sequences)

# Reorder selected_columns to match the order of row names 
reordered <- selected_columns[match(rownames(filtered_ancom_abundances), rownames(selected_columns)), ]

combined <- cbind(filtered_ancom_abundances, reordered)

combined$Diagnosis <- factor(combined$Diagnosis)

# Fit a multivariate regression model
model <- lm(cbind(Actinomyces_oris, Bifidobacterium_dentium, Streptococcus_mutans,
                  Anaerostipes_hadrus, Blautia_wexlerae, Eisenbergiella_tayi,
                  Roseburia_intestinalis, Clostridium_leptum, Ruminococcaceae_bacterium_D5, 
                  Ruminococcus_lactaris, Escherichia_coli) ~ ., data = combined)

# Print the summary of the model
summary(model)

saveRDS(summary(model), "Regression scripts/Wallen_model_summary_wo_confounders.rds")

# create a table of regression stats 
Wallen_tidy2 <- tidy(model)

# adjust the p values for multiple comparisons w/ Benjamini-Hochberg correction
p_values <- Wallen_tidy2[, 6]
p_values <- unlist(p_values)
adjusted_p_values <- p.adjust(p_values, method = "BH")
Wallen_tidy2$adjusted_p_values <- adjusted_p_values

# Save
write.xlsx(Wallen_tidy2, "Regression scripts/Wallen PD regression summary wo confounders fixed.xlsx")



# Extract taxa names 
taxa_names <- Wallen_tidy$response[Wallen_tidy$term == "DiagnosisPD"]

# Extract coefficients associated with the term "DiagnosisPD" from Wallen_tidy
coef_with_confounding <- Wallen_tidy$estimate[Wallen_tidy$term == "DiagnosisPD"]

# Extract coefficients associated with the term "DiagnosisPD" from Wallen_tidy2
coef_without_confounding <- Wallen_tidy2$estimate[Wallen_tidy2$term == "DiagnosisPD"]

# Calculate the difference in coefficients
difference <- coef_without_confounding - coef_with_confounding

# Calculate percent change
percent_change <- (difference / coef_with_confounding) * 100

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

write.xlsx(results, "Regression scripts/Wallen PD difference of confounders fixed.xlsx")
