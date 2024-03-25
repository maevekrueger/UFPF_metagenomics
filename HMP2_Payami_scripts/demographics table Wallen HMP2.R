# using the arsenal package to create simple demo tables for Wallen + HMP2 datasets 
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
PD_metadata <- readRDS("HMP2_Payami/Wallen PD Metadata.rds")

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

write.csv(table, 'HMP2_Payami/Figures/arsenal Wallen PD table.csv')
# to an HTML 
write2html(arsenal_table, "HMP2_Payami/Figures/arsenal Wallen PD table.html")
## write to a Word document
write2word(arsenal_table, "HMP2_Payami/Figures/arsenal Wallen PD table.doc", title="Demographics")

