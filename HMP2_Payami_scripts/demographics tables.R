# using the arsenal package instead 
library(arsenal)
library(tidyverse)

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

# combining IBD
order <- c("nonIBD", "IBD", "Overall")
IBD_metadata$diagnosis2 <- factor(IBD_metadata$diagnosis2, levels = order)

IBD_metadata <- IBD_metadata %>%
  mutate(Age = as.numeric(Age)) 

# reorder diagnosis groups
IBD_metadata$diagnosis2 <- factor(IBD_metadata$diagnosis2, levels = order)

IBD_metadata$Diagnosis <- factor(IBD_metadata$Diagnosis, levels = c("CD", "UC", "nonIBD"))
levels(IBD_metadata$Diagnosis)


arsenal_table <- tableby(diagnosis2 ~ Age + Sex + Diagnosis + Site + Immunosuppressants + 
                           Antibiotics + Chemotherapy + Colonoscopy_past_2_weeks + 
                           Diarrhea_past_2_weeks + Bowel_surgery,
                         data= IBD_metadata) 

summary(arsenal_table, text=TRUE, title='HMP2 Demographics', pfootnote = TRUE)     # can remove pfootnote if you don't want to include what statistical test was used

# lists the statistical test performed for each variable 
tests(arsenal_table)


# exporting table 
# to .CSV
table <- summary(tableby(diagnosis2 ~ Age + Sex + Diagnosis + Site + Immunosuppressants + 
                           Antibiotics + Chemotherapy + Colonoscopy_past_2_weeks + 
                           Diarrhea_past_2_weeks + Bowel_surgery, data=IBD_metadata), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/HMP2 IBD metadata table.csv')

## write to a Word document
write2word(arsenal_table, "HMP2_Payami/Figures/HMP2 IBD metadata table.doc", title="Demographics")


# --------------------------------------------------------------------------------------
# Payami PD Data 

# demographics table 
PD_metadata <- readRDS("HMP2_Payami/Payami PD Metadata.rds")

PD_metadata <- PD_metadata %>%
  rename(Age = Age_at_collection, Diagnosis = Case_status)

PD_metadata <- PD_metadata %>%
  mutate(Sex = ifelse(Sex == "M", "Male", ifelse(Sex == "F", "Female", Sex)))

# combining IBD
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

write.csv(table, 'HMP2_Payami/Figures/Wallen PD metadata table.csv')

## write to a Word document
write2word(arsenal_table, "HMP2_Payami/Figures/Wallen PD metadata table.doc", title="Demographics")

