library(tidyverse)
library(arsenal)

# version 2 with edited NSAID and anti-inflammatory categories 
Metadata <- readRDS("UFPF/Metadata.rds")

colnames(Metadata)[c(2, 12:53)] <- gsub("\\.", "_", colnames(Metadata)[c(2, 12:53)])

# combining IBD
order2 <- c("Control", "IBD", "PD", "Overall")
Metadata$Diagnosis2 <- factor(Metadata$Diagnosis2, levels = order2)

Metadata <- Metadata %>%
  mutate(Age = as.numeric(Age)) 

# reorder the Smoking variables 
Metadata$Smoking <- factor(Metadata$Smoking, levels = c("never", "less than 10 years", "between 10-20 years", "20+ years"))
levels(Metadata$Smoking)

# ---------------------------
# Change the label names 
labels(Metadata)  <- c(Age = 'Age, yrs', Race_Ethnicity = 'Race', Carbidopa_Levo = 'Carbidopa/Levodopa', 
                       Other_PD_meds = 'Other PD meds', Indigestion_meds = 'Indigestion meds', 
                       Anti_inflammatories__non_NSAID_ = 'Anti-inflammatories (other)', 
                       Blood_thinners = 'Blood thinners', Cholesterol_meds = 'Cholesterol meds', 
                       Blood_pressure_meds = 'Blood pressure meds', Thyroid_meds = 'Thyroid meds', 
                       Diabetes_meds = 'Diabetes meds', Depression_anxiety_meds = 'Depression/Anxiety meds', 
                       Iron_specific_supplement = 'Iron supplement', Sleep_aids = 'Sleep aids', 
                       NSAIDs = 'NSAIDs', Total_Caffeine_Intake = 'Caffeine Intake, mg/lifetime')

# Making demo table with the arsenal package 
arsenal_table <- tableby(Diagnosis2 ~ Age + Sex + Race_Ethnicity + Smoking + 
                           Carbidopa_Levo + Other_PD_meds + Laxatives + Indigestion_meds + 
                           Anti_TNF + NSAIDs + Anti_inflammatories__non_NSAID_ + Blood_thinners + 
                           Cholesterol_meds + Blood_pressure_meds + Thyroid_meds + 
                           Diabetes_meds + Depression_anxiety_meds + Iron_specific_supplement + 
                           Estrogen + Antihistamines + Sleep_aids + Total_Caffeine_Intake, data=Metadata)
summary(arsenal_table, text=TRUE, title='Demographics and Metadata', pfootnote = TRUE)     # can remove pfootnote if you don't want to include what statistical test was used

# lists the statistical test performed for each variable 
tests(arsenal_table)


# exporting table 
# to .CSV
table <- summary(tableby(Diagnosis2 ~ Age + Sex + Race_Ethnicity + Smoking + 
                           Carbidopa_Levo + Other_PD_meds + Laxatives + Indigestion_meds + 
                           Anti_TNF + NSAIDs + Anti_inflammatories__non_NSAID_ + Blood_thinners + 
                           Cholesterol_meds + Blood_pressure_meds + Thyroid_meds + 
                           Diabetes_meds + Depression_anxiety_meds + Iron_specific_supplement + 
                           Estrogen + Antihistamines + Sleep_aids + 
                           NSAIDs + Total_Caffeine_Intake, data=Metadata), text = NULL)
table
table <- as.data.frame(table)

table2 <- table
colnames(table2) <- c("x", "Control", "IBD", "PD", "Total", "p-value")

table2 %>%
  filter(x != "N") %>%
  mutate(x = str_remove_all(x, "Y")) -> table2

table2 <- table2[-c(27, 53), ]   
table2$x <- gsub("yes", "", table2$x)

table2 <- table2[, c(-5)]

colnames(table2) <- colnames(table)
table2

knitr::kable(table2, format = "html")

write.csv(table2, 'UFPF/Figures/arsenal demo table w TNF.csv')

## write to a Word document
write2word(arsenal_table, "UFPF/Figures/arsenal demo table w TNF.doc", title="Demographics and Metadata")

