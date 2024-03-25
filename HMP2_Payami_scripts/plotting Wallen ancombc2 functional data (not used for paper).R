library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)

# plotting Wallen PD functional ancombc2 data 
# did not end up using the KO group data for this paper but it is available for future analyses

PD_KO <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD ancombc2 KO groups output.rds")
PD_path <- readRDS("HMP2_Payami/ANCOMBC2/Wallen PD ancombc2 metacyc pathways.rds")

prim_KO <- PD_KO$res        # 3575 KO groups 
prim_path <- PD_path$res    # 356 pathways 

sig_taxa_KO <- prim_KO %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))      # 1247 PD significant KO groups

sig_taxa_path <- prim_path %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_Case_statusPD"))))      # 198 PD significant pathways


# Pathway Data 
sig_taxa_path <- sig_taxa_path[, c(1, 3, 7, 19, 23)]

sig_taxa_path <- sig_taxa_path[order(sig_taxa_path$q_Case_statusPD), ]

sig_taxa_path <- sig_taxa_path[-7, ]

sig_taxa_long <- sig_taxa_path %>%
  pivot_longer(
    cols = starts_with("lfc_"),  
    names_to = "Comparison",    
    values_to = "LFC"           
  )

# 
sig_taxa_long <- sig_taxa_long %>%
  select(-2:-3)

# create separate data frame with significance info 
diff_columns <- sig_taxa_path %>%
  select(starts_with("taxon"), starts_with("diff_Case_statusPD"))

# Convert to long format
diff_columns <- diff_columns %>%
  pivot_longer(
    cols = starts_with("diff_"),  
    names_to = "Comparison",    
    values_to = "Significance"  
  )

# create separate data frame with q values  
q_columns <- sig_taxa_path %>%
  select(starts_with("taxon"), starts_with("q_"))

# Convert to long format
q_columns <- q_columns %>%
  pivot_longer(
    cols = starts_with("q_"),  
    names_to = "Comparison",    
    values_to = "Adj P Value"  
  )

# create separate data frame with se values  
se_columns <- sig_taxa_path %>%
  select(starts_with("taxon"), starts_with("se_"))

# Convert to long format
se_columns <- se_columns %>%
  pivot_longer(
    cols = starts_with("se_"), 
    names_to = "Comparison",   
    values_to = "Standard_error"  
  )

# adding TRUE/FALSE info 
sig_taxa_long$'Adj P Value' <- q_columns$'Adj P Value'
sig_taxa_long$Significance <- diff_columns$Significance
sig_taxa_long$Standard_error <- se_columns$Standard_error

# plotting ALL SIGNIFICANT PATHWAYS ( way too many for a clean visual )
ggplot(sig_taxa_long, aes(x = LFC, y = taxon)) +
  geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue3") +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "none", 
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  ) +
  labs(
    title = "Pairwise Differential MetaCyc Pathway Analysis",
    x = "Log Fold Change",
    y = "MetaCyc Pathways",
    fill = "Comparison"
  ) 


# getting top (most significant paths) for each group 
# Rename groups in the "Comparison" column
sig_taxa_long <- sig_taxa_long %>%
  mutate(Comparison = case_when(
    Comparison == "lfc_Case_statusPD" ~ "PD vs Control",
    TRUE ~ Comparison  # Keep other values as they are
  ))

# Create separate data frames for each group
PD_10 <- head(sig_taxa_long, 10)
PD_20 <- head(sig_taxa_long, 20)

# add back matching KOs and other group's LFCs 
PD_10 <- PD_10 %>%
  left_join(sig_taxa_long, by = "taxon")
PD_10 <- PD_10 %>%
  select(-2:-7)
PD_10 <- PD_10 %>%
  rename_at(vars(2:7), ~ gsub(".y", "", .))

PD_20 <- PD_20 %>%
  left_join(sig_taxa_long, by = "taxon")
PD_20 <- PD_20 %>%
  select(-2:-7)
PD_20 <- PD_20 %>%
  rename_at(vars(2:7), ~ gsub(".y", "", .))

# plotting 
# Create a list of data frames
df_list <- list(PD_10, PD_20)
titles <- c("Top 10", "Top 20")


plot_list <- list()  # empty list to store the plots

for (i in 1:length(df_list)) {
  plot <- ggplot(df_list[[i]], aes(x = LFC, y = fct_reorder(taxon, as.numeric(`Adj P Value`)))) +
    geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue3") +
    geom_errorbarh(
      aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
      position = position_dodge(0.9),
      height = 0.25,  
      size = 0.5,   
      color = "gray73"  
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 15, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
    ) +
    labs(
      title = paste("PD Associated MetaCyc Pathways -", titles[i]),
      x = "Log Fold Change",
      y = "MetaCyc Pathways",
      fill = "Comparison"
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))  
  
  plot_list[[i]] <- plot  # Store the plot in the list
}

plot1 <- plot_list[[1]]
plot2 <- plot_list[[2]]
plot1 
plot2

# Save plot1 as an image file
ggsave("HMP2_Payami/Figures/Wallen PD top 10 sig metacyc pathways se bars.png", plot1, dpi = 800, units = "in", height = 10, width = 12)

# Save plot2 as an image file
ggsave("HMP2_Payami/Figures/Wallen PD top 20 sig metacyc pathways se bars.png", plot2, dpi = 800, units = "in", height = 10, width = 12)




# ------------------------------------------------------------------------------
# KO groups  
sig_taxa_KO <- sig_taxa_KO[, c(1, 3, 7, 19, 23)]

sig_taxa_KO <- sig_taxa_KO[order(sig_taxa_KO$q_Case_statusPD), ]

# remove meaningless rows 
sig_taxa_KO <- sig_taxa_KO[-3, ]

sig_taxa_long <- sig_taxa_KO %>%
  pivot_longer(
    cols = starts_with("lfc_"),  
    names_to = "Comparison",    
    values_to = "LFC"           
  )

# 
sig_taxa_long <- sig_taxa_long %>%
  select(-2:-3)

# create separate data frame with significance info 
diff_columns <- sig_taxa_KO %>%
  select(starts_with("taxon"), starts_with("diff_Case_statusPD"))

# Convert to long format
diff_columns <- diff_columns %>%
  pivot_longer(
    cols = starts_with("diff_"),  
    names_to = "Comparison",    
    values_to = "Significance"  
  )

# create separate data frame with q values  
q_columns <- sig_taxa_KO %>%
  select(starts_with("taxon"), starts_with("q_"))

# Convert to long format
q_columns <- q_columns %>%
  pivot_longer(
    cols = starts_with("q_"),  
    names_to = "Comparison",    
    values_to = "Adj P Value"  
  )

# create separate data frame with se values  
se_columns <- sig_taxa_KO %>%
  select(starts_with("taxon"), starts_with("se_"))

# Convert to long format
se_columns <- se_columns %>%
  pivot_longer(
    cols = starts_with("se_"),  
    names_to = "Comparison",    
    values_to = "Standard_error"  
  )

# adding TRUE/FALSE info 
sig_taxa_long$'Adj P Value' <- q_columns$'Adj P Value'
sig_taxa_long$Significance <- diff_columns$Significance
sig_taxa_long$Standard_error <- se_columns$Standard_error

# plotting ALL SIGNIFICANT KO GROUPS ( way too many for a clean visual )
ggplot(sig_taxa_long, aes(x = LFC, y = taxon)) +
  geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue3") +
  geom_errorbarh(
    aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
    position = position_dodge(0.9),
    height = 0.25,  
    size = 1.2,    
    color = "gray"  
  ) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "none", 
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
  ) +
  labs(
    title = "Pairwise Differential KO Group Analysis",
    x = "Log Fold Change",
    y = "KO Groups",
    fill = "Comparison"
  ) 


# getting top (most significant KOs) for each group 
# Rename groups in the "Comparison" column
sig_taxa_long <- sig_taxa_long %>%
  mutate(Comparison = case_when(
    Comparison == "lfc_Case_statusPD" ~ "PD vs Control",
    TRUE ~ Comparison  # Keep other values as they are
  ))

# Create separate data frames for each group
PD_10 <- head(sig_taxa_long, 10)
PD_20 <- head(sig_taxa_long, 20)

# add back matching KOs and other group's LFCs 
PD_10 <- PD_10 %>%
  left_join(sig_taxa_long, by = "taxon")
PD_10 <- PD_10 %>%
  select(-2:-7)
PD_10 <- PD_10 %>%
  rename_at(vars(2:7), ~ gsub(".y", "", .))

PD_20 <- PD_20 %>%
  left_join(sig_taxa_long, by = "taxon")
PD_20 <- PD_20 %>%
  select(-2:-7)
PD_20 <- PD_20 %>%
  rename_at(vars(2:7), ~ gsub(".y", "", .))

# plotting 
# Create a list of data frames
df_list <- list(PD_10, PD_20)
titles <- c("Top 10", "Top 20")

plot_list <- list()  # empty list to store the plots

for (i in 1:length(df_list)) {
  plot <- ggplot(df_list[[i]], aes(x = LFC, y = fct_reorder(taxon, as.numeric(`Adj P Value`)))) +
    geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue3") +
    geom_errorbarh(
      aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
      position = position_dodge(0.9),
      height = 0.25, 
      size = 0.5,    
      color = "gray73"  
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      strip.text = element_text(color = "white", face = "bold", size = rel(1.5))
    ) +
    labs(
      title = paste("PD Associated KO Groups -", titles[i]),
      x = "Log Fold Change",
      y = "KO Groups",
      fill = "Comparison"
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))  
  
  plot_list[[i]] <- plot  # Store the plot in the list
}

plot1 <- plot_list[[1]]
plot2 <- plot_list[[2]]
plot1 
plot2

# Save plot1 as an image file
ggsave("HMP2_Payami/Figures/Wallen PD top 10 sig KO groups se bars.png", plot1, dpi = 800, units = "in", height = 10, width = 12)

# Save plot2 as an image file
ggsave("HMP2_Payami/Figures/Wallen PD top 20 sig KO groups se bars.png", plot2, dpi = 800, units = "in", height = 10, width = 12)

