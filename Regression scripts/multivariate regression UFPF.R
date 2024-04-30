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
