
# get bias corrected abundances from ancombc2 output to be used later in the module analysis 
#(count data put into ancombc2 included the unclassified estimation)

# Wallen PD data set 
# SPECIES LEVEL 
ancom3 <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD species ancombc2 output.rds")   
data <- ancom3$samp_frac
features <- ancom3$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

PD_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(PD_abundances, "HMP2_Payami/ANCOMBC2/Wallen PD species bias corrected abund.rds")


# GENUS LEVEL 
ancom_g <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD genus ancombc2 output.rds")

data <- ancom_g$samp_frac
features <- ancom_g$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

PD_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(PD_abundances, "HMP2_Payami/ANCOMBC2/Wallen PD genus bias corrected abund.rds")


# pathways
ancom_path <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD ancombc2 metacyc pathways.rds")

data <- ancom_path$samp_frac
features <- ancom_path$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

PD_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(PD_abundances, "HMP2_Payami/ANCOMCB2/Wallen PD paths bias corrected abund ancombc2.rds")


# KO groups 
# KO group data was not used for this analysis but it is available to use for future analyses
ancom_ko <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD ancombc2 KO groups output.rds")

data <- ancom_ko$samp_frac
features <- ancom_ko$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

PD_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(PD_abundances, "HMP2_Payami/ANCOMBC2/Wallen PD KO bias corrected abund ancombc2.rds")


# -----------------------------------------------------------------------
# HMP2 IBD data set 
# SPECIES LEVEL 
ancom <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD species ancombc2 output.rds") 

data <- ancom$samp_frac
features <- ancom$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

IBD_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(IBD_abundances, "HMP2_Payami/ANCOMBC2/HMP2 IBD species bias corrected abund.rds")


# GENUS LEVEL 
ancom_g <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD genus ancombc2 output.rds") 

data <- ancom_g$samp_frac
features <- ancom_g$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

IBD_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(IBD_abundances, "HMP2_Payami/ANCOMBC2/HMP2 IBD genus bias corrected abund.rds")


# KO groups 
ancom_ko <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD KO ancombc2 output.rds") 

data <- ancom_ko$samp_frac
features <- ancom_ko$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

IBD_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(IBD_abundances, "HMP2_Payami/ANCOMBC2/HMP2 IBD KO bias corrected abund ancombc2.rds")


# PATHWAYS 
ancom_p <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD pathway ancombc2 output.rds") 

data <- ancom_p$samp_frac
features <- ancom_p$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

IBD_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(IBD_abundances, "HMP2_Payami/ANCOMBC2/HMP2 IBD pathway bias corrected abund.rds")
