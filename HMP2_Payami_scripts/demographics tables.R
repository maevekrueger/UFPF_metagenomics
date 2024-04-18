# using the arsenal package instead 
library(arsenal)
library(tidyverse)
 
# HMP2 IBD 
# demographics table 
IBD_metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")
IBD_metadata <- column_to_rownames(IBD_metadata, var = "Sample")

# removing some columns I don't need
columns_to_remove <- c(14:19)
IBD_metadata <- IBD_metadata[, -columns_to_remove]

IBD_metadata <- IBD_metadata %>%
  rename(Age = consent_age, Diagnosis = diagnosis, Sex = sex, 
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

# reorder diagnosis groups
IBD_metadata$diagnosis2 <- factor(IBD_metadata$diagnosis2, levels = order)

IBD_metadata$Diagnosis <- factor(IBD_metadata$Diagnosis, levels = c("CD", "UC", "nonIBD"))
levels(IBD_metadata$Diagnosis)


# Filter out samples with missing values in specific variables
# recombine later in excel
IBD_metadata_filtered <- IBD_metadata %>%
  filter(!is.na(Probiotic))  

arsenal_table <- tableby(diagnosis2 ~ Age + Sex + Diagnosis + Site + Probiotic + Immunosuppressants + 
                           Antibiotics + Chemotherapy + Colonoscopy_past_2_weeks + 
                           Diarrhea_past_2_weeks + Bowel_surgery,
                         data = IBD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'HMP2 Demographics', pfootnote = TRUE)

table <- summary(tableby(diagnosis2 ~ Age + Sex + Diagnosis + Site + Probiotic + Immunosuppressants + 
                           Antibiotics + Chemotherapy + Colonoscopy_past_2_weeks + 
                           Diarrhea_past_2_weeks + Bowel_surgery, data=IBD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/HMP2 IBD metadata table 1.csv')


# bowel_surgery
IBD_metadata_filtered <- IBD_metadata %>%
  filter(Bowel_surgery != "") 

arsenal_table <- tableby(diagnosis2 ~ Age + Sex + Diagnosis + Bowel_surgery,
                         data = IBD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'HMP2 Demographics', pfootnote = TRUE)

table <- summary(tableby(diagnosis2 ~ Age + Sex + Diagnosis + Bowel_surgery, data=IBD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/HMP2 IBD metadata table 2.csv')


# colonscopy 
IBD_metadata_filtered <- IBD_metadata %>%
  filter(Colonoscopy_past_2_weeks != "") 

arsenal_table <- tableby(diagnosis2 ~ Age + Sex + Diagnosis + Colonoscopy_past_2_weeks,
                         data = IBD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'HMP2 Demographics', pfootnote = TRUE)

table <- summary(tableby(diagnosis2 ~ Age + Sex + Diagnosis + Colonoscopy_past_2_weeks, data=IBD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/HMP2 IBD metadata table 3.csv')


# diarrhea
IBD_metadata_filtered <- IBD_metadata %>%
  filter(Diarrhea_past_2_weeks != "") 

arsenal_table <- tableby(diagnosis2 ~ Age + Sex + Diagnosis + Diarrhea_past_2_weeks,
                         data = IBD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'HMP2 Demographics', pfootnote = TRUE)

table <- summary(tableby(diagnosis2 ~ Age + Sex + Diagnosis + Diarrhea_past_2_weeks, data=IBD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/HMP2 IBD metadata table 4.csv')

# -------------------------------------------------------------------------------
# Wallen PD 
# demographics table 
PD_metadata <- readRDS("HMP2_Payami/Wallen PD Metadata.rds")

PD_metadata <- PD_metadata %>%
  rename(Age = Age_at_collection, Diagnosis = Case_status)

PD_metadata <- PD_metadata %>%
  mutate(Sex = ifelse(Sex == "M", "Male", ifelse(Sex == "F", "Female", Sex)))

order <- c("Control", "PD")
PD_metadata$Diagnosis <- factor(PD_metadata$Diagnosis, levels = order)

arsenal_table <- tableby(Diagnosis ~ Age + Sex + Constipation + Diarrhea +
                           IBS + IBD + SIBO + Celiac_disease + Antibiotics_current +
                           Antibiotics_past_3_months + Laxatives + Indigestion_drugs +
                           Anti_inflammatory_drugs + Radiation_Chemo + Blood_thinners +
                           Cholesterol_med + Blood_pressure_med + Thyroid_med +
                           Asthma_or_COPD_med + Diabetes_med + Pain_med + 
                           Depression_anxiety_mood_med + Birth_control_or_estrogen +
                           Antihistamines + Probiotic + Co_Q_10 + Sleep_aid, data= PD_metadata) 

summary(arsenal_table, text=TRUE, title='Wallen PD Demographics', pfootnote = TRUE)     # can remove pfootnote if you don't want to include what statistical test was used

# lists the statistical test performed for each variable 
tests(arsenal_table)

# exporting table 
# to .CSV
table <- summary(tableby(Diagnosis ~ Age + Sex + Constipation + Diarrhea +
                           IBS + IBD + SIBO + Celiac_disease + Antibiotics_current +
                           Antibiotics_past_3_months + Laxatives + Indigestion_drugs +
                           Anti_inflammatory_drugs + Radiation_Chemo + Blood_thinners +
                           Cholesterol_med + Blood_pressure_med + Thyroid_med +
                           Asthma_or_COPD_med + Diabetes_med + Pain_med + 
                           Depression_anxiety_mood_med + Birth_control_or_estrogen +
                           Antihistamines + Probiotic + Co_Q_10 + Sleep_aid, data=PD_metadata), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(Constipation))  

arsenal_table <- tableby(Diagnosis ~ Constipation,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ Constipation, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 2.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(Diarrhea))  

arsenal_table <- tableby(Diagnosis ~ Diarrhea,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ Diarrhea, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 3.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(IBS))  

arsenal_table <- tableby(Diagnosis ~ IBS,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ IBS, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 4.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(IBD))  

arsenal_table <- tableby(Diagnosis ~ IBD,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ IBD, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 5.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(SIBO))  

arsenal_table <- tableby(Diagnosis ~ SIBO,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ SIBO, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 6.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(Celiac_disease))  

arsenal_table <- tableby(Diagnosis ~ Celiac_disease,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ Celiac_disease, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 7.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(Antibiotics_current))  

arsenal_table <- tableby(Diagnosis ~ Antibiotics_current,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ Antibiotics_current, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 7.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(Antibiotics_past_3_months))  

arsenal_table <- tableby(Diagnosis ~ Antibiotics_past_3_months,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ Antibiotics_past_3_months, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 7.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(Laxatives))  

arsenal_table <- tableby(Diagnosis ~ Laxatives,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ Laxatives, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 7.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(Indigestion_drugs))  

arsenal_table <- tableby(Diagnosis ~ Indigestion_drugs,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ Indigestion_drugs, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 7.csv')

# remove missing 
PD_metadata_filtered <- PD_metadata %>%
  filter(!is.na(Birth_control_or_estrogen))  

arsenal_table <- tableby(Diagnosis ~ Birth_control_or_estrogen,
                         data = PD_metadata_filtered)
summary(arsenal_table, text = TRUE, title = 'Wallen PD Demographics', pfootnote = TRUE)

table <- summary(tableby(Diagnosis ~ Birth_control_or_estrogen, data=PD_metadata_filtered), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/Old/Wallen PD metadata table 7.csv')

# ---- continued for rest of variables to confirm stats 