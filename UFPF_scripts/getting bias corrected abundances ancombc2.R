# UFPF 
# get bias corrected abundances from ancombc2 output 
# WITH unclassified estimation 
# SPECIES
ancom2 <- readRDS("UFPF/ANCOMBC2/ancombc2 species.rds") 

data <- ancom2$samp_frac
features <- ancom2$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

UFPF_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(UFPF_abundances, "UFPF/ANCOMBC2/bias corrected abund species.rds")


# GENUS-level
# get bias corrected abundances from ancombc2 output from running at the genus-level
ancom2 <- readRDS("UFPF/ANCOMBC2/ancombc2 genus.rds")  

data <- ancom2$samp_frac
features <- ancom2$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

UFPF_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(UFPF_abundances, "UFPF/ANCOMBC2/bias corrected abund genus.rds")



# PATHWAYS
# get bias corrected abundances from ancombc2 output from running at the genus-level
ancom2 <- readRDS("UFPF/ANCOMBC2/ancombc2 pathways output.rds")  

data <- ancom2$samp_frac
features <- ancom2$feature_table

# Add pesudo-count to avoid taking the log of 0
log_obs_abn = log(features + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - data)

UFPF_abundances <- as.data.frame(t(log_corr_abn))

saveRDS(UFPF_abundances, "UFPF/ANCOMBC2/bias corrected abund pathways.rds")


