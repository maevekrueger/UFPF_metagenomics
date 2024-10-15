library(tidyverse)
library(broom)
library(openxlsx)

# exmaining the influence of potential confounding variables 

# UFPF data 
# --> significant variables: Age, indigestion meds, anti-TNF drugs, Anti-inflammatory meds, 
#                            depression/anxiety meds, iron supplements 

# specify response variables 
ancom_abundances <- readRDS("UFPF/ANCOMBC2/bias corrected abund genus.rds")
ancom_abundances <- ancom_abundances[, c("Klebsiella", "Faecalimonas")]

# read in metadata to pull out predictor variables 
Metadata <- readRDS("UFPF/Metadata.rds")
colnames(Metadata)[c(2,6,12:53)] <- gsub("\\.", "_", colnames(Metadata)[c(2,6,12:53)])

# Scaling 
Metadata$Reads <- scale(Metadata$Reads)
Metadata$Age <- scale(Metadata$Age)

# pull out predictor variables according to the demographics/metadata table stats 
selected_columns <- Metadata %>%
  select(Sex, Age, Diagnosis2, Reads, Indigestion_meds, Anti_TNF, Anti_inflammatories__non_NSAID_, Depression_anxiety_meds, Iron_specific_supplement)
selected_columns <- selected_columns %>%
  rename(Anti_inflammatories = Anti_inflammatories__non_NSAID_)
selected_columns <- selected_columns %>%
  rename(Diagnosis = Diagnosis2)

combined <- cbind(ancom_abundances, selected_columns)

combined$Sex <- factor(combined$Sex)
combined$Diagnosis <- factor(combined$Diagnosis)

convert_values <- function(x) {
  x <- ifelse(x %in% c("N", "N ", "never"), 0,
              ifelse(x %in% c("Y", "Y ", "yes"), 1, NA))
  return(x)
}
combined <- combined %>%
  mutate_at(vars(7:11), list(~convert_values(.)))


# Fit a multivariate regression model
model <- lm(cbind(Klebsiella, Faecalimonas) ~ ., data = combined)

# Print the summary of the model
summary(model)

saveRDS(summary(model), "Regression scripts/UFPF_model_summary.rds")

# create a table of regression stats 
UFPF_tidy <- tidy(model)

# adjust the p values for multiple comparisons w/ Benjamini-Hochberg correction
p_values <- UFPF_tidy[, 6]
p_values <- unlist(p_values)
adjusted_p_values <- p.adjust(p_values, method = "BH")
UFPF_tidy$adjusted_p_values <- adjusted_p_values

# Save
write.xlsx(UFPF_tidy, "Regression scripts/UFPF regression summary.xlsx")



# -----------------------------------------------------------------------------
# without confounders 
ancom_abundances <- readRDS("UFPF/ANCOMBC2/bias corrected abund genus.rds")
ancom_abundances <- ancom_abundances[, c("Klebsiella", "Faecalimonas")]

# read in metadata to pull out predictor variables 
Metadata <- readRDS("UFPF/Metadata.rds")
colnames(Metadata)[c(2,6,12:53)] <- gsub("\\.", "_", colnames(Metadata)[c(2,6,12:53)])

# Scaling 
Metadata$Reads <- scale(Metadata$Reads)

# pull out predictor variables according to the demographics/metadata table stats 
selected_columns <- Metadata %>%
  select(Diagnosis2, Reads)
selected_columns <- selected_columns %>%
  rename(Diagnosis = Diagnosis2)

combined <- cbind(ancom_abundances, selected_columns)

combined$Diagnosis <- factor(combined$Diagnosis)

# Fit a multivariate regression model
model <- lm(cbind(Klebsiella, Faecalimonas) ~ ., data = combined)

# Print the summary of the model
summary(model)

saveRDS(summary(model), "Regression scripts/UFPF_model_summary_wo_confounders.rds")

# create a table of regression stats 
UFPF_tidy2 <- tidy(model)

# adjust the p values for multiple comparisons w/ Benjamini-Hochberg correction
p_values <- UFPF_tidy2[, 6]
p_values <- unlist(p_values)
adjusted_p_values <- p.adjust(p_values, method = "BH")
UFPF_tidy2$adjusted_p_values <- adjusted_p_values

# Save
write.xlsx(UFPF_tidy2, "Regression scripts/UFPF regression summary without confounders.xlsx")



# Extract taxa names from UFPF_tidy
taxa_names <- UFPF_tidy$response[UFPF_tidy$term == "DiagnosisIBD"]

# Extract coefficients associated with the term "DiagnosisPD" from UFPF_tidy
coef_with_confounding <- UFPF_tidy$estimate[UFPF_tidy$term == "DiagnosisIBD"]

# Extract coefficients associated with the term "DiagnosisPD" from UFPF_tidy2
coef_without_confounding <- UFPF_tidy2$estimate[UFPF_tidy2$term == "DiagnosisIBD"]

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

print(results)

# Save
write.xlsx(results, "Regression scripts/UFPF IBD difference of confounders.xlsx")
