# using the arsenal package instead 
library(arsenal)
library(knitr)
library(dplyr)
library(tibble)

# demographics table 
IBD_metadata <- readRDS("HMP2_Payami/IBD Metadata Age Filtered.rds")

IBD_metadata <- column_to_rownames(IBD_metadata, var = "Sample")

# removing some columns I don't need
columns_to_remove <- c(14:19)
IBD_metadata <- IBD_metadata[, -columns_to_remove]

IBD_metadata <- IBD_metadata %>%
  rename(Age = consent_age, Diagnosis = diagnosis, Sex = sex)

# combining IBD
order <- c("nonIBD", "IBD", "Overall")
IBD_metadata$diagnosis2 <- factor(IBD_metadata$diagnosis2, levels = order)

IBD_metadata <- IBD_metadata %>%
  mutate(Age = as.numeric(Age)) 

# reorder diagnosis groups
IBD_metadata$diagnosis2 <- factor(IBD_metadata$diagnosis2, levels = order)

arsenal_table <- tableby(diagnosis2 ~ Age + Sex, data= IBD_metadata) 
summary(arsenal_table, text=TRUE, title='HMP2 Demographics', pfootnote = TRUE)     # can remove pfootnote if you don't want to include what statistical test was used

# lists the statistical test performed for each variable 
tests(arsenal_table)


# exporting table 
# to .CSV
table <- summary(tableby(diagnosis2 ~ Age + Sex, data=IBD_metadata), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/arsenal HMP2 IBD combined table.csv')
# to an HTML 
write2html(arsenal_table, "HMP2_Payami/Figures/arsenal HMP2 IBD combined table.html")
## write to a Word document
write2word(arsenal_table, "HMP2_Payami/Figures/arsenal HMP2 IBD combined table.doc", title="Demographics")


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

arsenal_table <- tableby(Diagnosis ~ Age + Sex, data= PD_metadata) 
summary(arsenal_table, text=TRUE, title='UAB PD Demographics', pfootnote = TRUE)     # can remove pfootnote if you don't want to include what statistical test was used

# lists the statistical test performed for each variable 
tests(arsenal_table)

# exporting table 
# to .CSV
table <- summary(tableby(Diagnosis ~ Age + Sex, data=PD_metadata), text = NULL)
table <- as.data.frame(table)

write.csv(table, 'HMP2_Payami/Figures/arsenal Payami PD table.csv')
# to an HTML 
write2html(arsenal_table, "HMP2_Payami/Figures/arsenal Payami PD table.html")
## write to a Word document
write2word(arsenal_table, "HMP2_Payami/Figures/arsenal Payami PD table.doc", title="Demographics")



# adding in metadata (semi-matching UFPF format)
# Change the label names 
labels(PD_metadata)  <- c(Age = 'Age, yrs', Celiac_disease = 'Celiac Disease', 
                       Antibiotics_current = 'Current Antibiotics', Indigestion_drugs = 'Indigestion meds', 
                       Anti_inflammatory_drugs = 'Anti-inflammatories', Radition_chemo = 'Radiation or Chemo',
                       Blood_thinners = 'Blood thinners', Cholesterol_med = 'Cholesterol meds', 
                       Blood_pressure_med = 'Blood pressure meds', Thyroid_med = 'Thyroid meds', 
                       Diabetes_med = 'Diabetes meds', Depression_anxiety_mood_med = 'Depression/Anxiety meds', 
                       Sleep_aid = 'Sleep aids')

# Making demo table with the arsenal package 
arsenal_table <- tableby(Diagnosis ~ Age + Sex + Race + Constipation + Diarrhea + IBS + IBD +
                           SIBO + Celiac_disease + Laxatives + Indigestion_drugs + 
                           Anti_inflammatory_drugs + Blood_thinners + 
                           Cholesterol_med + Blood_pressure_med + Thyroid_med + 
                           Diabetes_med + Depression_anxiety_mood_med + Probiotic + 
                           Sleep_aid, data=PD_metadata)
summary(arsenal_table, text=TRUE, title='Demographics and Metadata', pfootnote = TRUE)     # can remove pfootnote if you don't want to include what statistical test was used


table <- summary(tableby(Diagnosis ~ Age + Sex + Race + Constipation + Diarrhea + IBS + IBD +
                           SIBO + Celiac_disease + Laxatives + Indigestion_drugs + 
                           Anti_inflammatory_drugs + Blood_thinners + 
                           Cholesterol_med + Blood_pressure_med + Thyroid_med + 
                           Diabetes_med + Depression_anxiety_mood_med + Probiotic + 
                           Sleep_aid, data=PD_metadata), text = NULL)
table <- as.data.frame(table)

table2 <- table
colnames(table2) <- c("x", "control", "pd", "total", "p")

table2 %>%
  filter(x != "N") %>%
  mutate(x = str_remove_all(x, "Y")) -> table2

# removing N-Miss rows
table2 <- table2 %>%
  slice(-c(8, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63))

colnames(table2) <- colnames(table)

write.csv(table2, 'HMP2_Payami/Figures/arsenal Payami demo metadata table.csv')

# to an HTML 
write2html(table2, "HMP2_Payami/Figures/arsenal Payami PD metadata table.html")

