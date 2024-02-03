library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)

ancom_ko <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD combined KO ancombc2 output.rds")

res_prim_ko = ancom_ko$res    # 2767 KOs examined 

sig_ko <- res_prim_ko %>%     # 1703 significant, IBD associated 
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))

# KO groups first 
sig_ko <- sig_ko[, c(1, 3, 7, 19, 23)]

sig_ko <- sig_ko[order(sig_ko$q_diagnosis2IBD), ]

sig_ko_long <- sig_ko %>%
  pivot_longer(
    cols = starts_with("lfc_"),  
    names_to = "Comparison",    
    values_to = "LFC"           
  )

sig_ko_long <- sig_ko_long[, -c(2:4)]

# create separate data frame with significance info 
diff_columns <- sig_ko %>%
  select(starts_with("taxon"), starts_with("diff_diagnosis2IBD"))

# Convert to long format
diff_columns <- diff_columns %>%
  pivot_longer(
    cols = starts_with("diff_"),  # Columns starting with "diff_"
    names_to = "Comparison",    # New column name for comparisons
    values_to = "Significance"  # New column name for significance values
  )

# create separate data frame with q values  
q_columns <- sig_ko %>%
  select(starts_with("taxon"), starts_with("q_"))

# Convert to long format
q_columns <- q_columns %>%
  pivot_longer(
    cols = starts_with("q_"),  # Columns starting with "q_"
    names_to = "Comparison",    # New column name for comparisons
    values_to = "Adj P Value"  # New column name for significance values
  )

# create separate data frame with se values  
se_columns <- sig_ko %>%
  select(starts_with("taxon"), starts_with("se_"))

# Convert to long format
se_columns <- se_columns %>%
  pivot_longer(
    cols = starts_with("se_"),  # Columns starting with "q_"
    names_to = "Comparison",    # New column name for comparisons
    values_to = "Standard_error"  # New column name for significance values
  )

# adding TRUE/FALSE info 
sig_ko_long$'Adj P Value' <- q_columns$'Adj P Value'
sig_ko_long$Significance <- diff_columns$Significance
sig_ko_long$Standard_error <- se_columns$Standard_error

# plotting ALL SIGNIFICANT KO GROUPS ( way too many for a clean visual )
ggplot(sig_ko_long, aes(x = LFC, y = taxon)) +
  geom_bar(stat = "identity", position = "dodge", fill = "green4") +
  geom_errorbarh(
    aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
    position = position_dodge(0.9),
    height = 0.25,  # Adjust the height of the error bars
    size = 1.2,    # Adjust the size of the error bars
    color = "gray"  # Set the color of the error bars
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
    title = "HMP2 IBD TOp KO Groups",
    x = "Log Fold Change",
    y = "KO Groups",
    fill = "Comparison"
  ) 


# getting top (most significant KOs) for each group 
# Rename groups in the "Comparison" column
sig_ko_long <- sig_ko_long %>%
  mutate(Comparison = case_when(
    Comparison == "lfc_diagnosis2IBD" ~ "IBD vs nonIBD",
    TRUE ~ Comparison  # Keep other values as they are
  ))

# Create separate data frames for each group
IBD_10 <- head(sig_ko_long, 10)
IBD_20 <- head(sig_ko_long, 20)


# plotting 
# Create a list of data frames
df_list <- list(IBD_10, IBD_20)

# Create a list of corresponding titles
titles <- c("Top 10", "Top 20")


plot_list <- list()  # Initialize an empty list to store the plots

for (i in 1:length(df_list)) {
  plot <- ggplot(df_list[[i]], aes(x = LFC, y = fct_reorder(taxon, as.numeric(`Adj P Value`)))) +
    geom_bar(stat = "identity", position = "dodge", fill = "darkseagreen4") +
    geom_errorbarh(
      aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
      position = position_dodge(0.9),
      height = 0.25,  # Adjust the height of the error bars
      size = 0.5,    # Adjust the size of the error bars
      color = "gray73"  # Set the color of the error bars
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
      title = paste("HMP2 IBD Associated KO Groups -", titles[i]),
      x = "Log Fold Change",
      y = "KO Groups",
      fill = "Comparison"
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))  # Adjust the width as needed
  
  plot_list[[i]] <- plot  # Store the plot in the list
}

plot1 <- plot_list[[1]]
plot2 <- plot_list[[2]]
plot1 
plot2

# Save plot1 as an image file
ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 IBD 10 sig KO.png", plot1, dpi = 800, units = "in", height = 10, width = 12)

# Save plot2 as an image file
ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 IBD top 20 sig KO.png", plot2, dpi = 800, units = "in", height = 10, width = 12)




# ------------------------------------------------------------------------------------------
# plotting pathway data IBD vs nonIBD 
ancom_path <- readRDS("HMP2_Payami/ANCOMBC2/HMP2 IBD combined pathway ancombc2 output.rds")

res_prim_path = ancom_path$res      # 294 pathways examined 

sig_path <- res_prim_path %>%
  rowwise() %>%
  filter(any(c_across(starts_with("diff_diagnosis2IBD"))))     # 158 pathways significant 


sig_path <- sig_path[, c(1, 3, 7, 19, 23)]

sig_path <- sig_path[order(sig_path$q_diagnosis2IBD), ]

sig_path_long <- sig_path %>%
  pivot_longer(
    cols = starts_with("lfc_"),  
    names_to = "Comparison",    
    values_to = "LFC"           
  )


sig_path_long <- sig_path_long[, -c(2:4)]

# create separate data frame with significance info 
diff_columns <- sig_path %>%
  select(starts_with("taxon"), starts_with("diff_diagnosis2IBD"))

# Convert to long format
diff_columns <- diff_columns %>%
  pivot_longer(
    cols = starts_with("diff_"),  # Columns starting with "diff_"
    names_to = "Comparison",    # New column name for comparisons
    values_to = "Significance"  # New column name for significance values
  )

# create separate data frame with q values  
q_columns <- sig_path %>%
  select(starts_with("taxon"), starts_with("q_"))

# Convert to long format
q_columns <- q_columns %>%
  pivot_longer(
    cols = starts_with("q_"),  # Columns starting with "q_"
    names_to = "Comparison",    # New column name for comparisons
    values_to = "Adj P Value"  # New column name for significance values
  )

# create separate data frame with se values  
se_columns <- sig_path %>%
  select(starts_with("taxon"), starts_with("se_"))

# Convert to long format
se_columns <- se_columns %>%
  pivot_longer(
    cols = starts_with("se_"),  # Columns starting with "q_"
    names_to = "Comparison",    # New column name for comparisons
    values_to = "Standard_error"  # New column name for significance values
  )

# adding TRUE/FALSE info 
sig_path_long$'Adj P Value' <- q_columns$'Adj P Value'
sig_path_long$Significance <- diff_columns$Significance
sig_path_long$Standard_error <- se_columns$Standard_error

# plotting ALL SIGNIFICANT PATHWAYS ( way too many for a clean visual )
ggplot(sig_path_long, aes(x = LFC, y = taxon)) +
  geom_bar(stat = "identity", position = "dodge", fill = "green4") +
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


# getting top (most significant Paths) for each group 
# Rename groups in the "Comparison" column
sig_path_long <- sig_path_long %>%
  mutate(Comparison = case_when(
    Comparison == "lfc_diagnosis2IBD" ~ "IBD vs nonIBD",
    TRUE ~ Comparison  # Keep other values as they are
  ))

# Create separate data frames for each group
IBD_10 <- head(sig_path_long, 10)
IBD_20 <- head(sig_path_long, 20)

# plotting 
# Create a list of data frames
df_list <- list(IBD_10, IBD_20)

# Create a list of corresponding titles
titles <- c("Top 10", "Top 20")


plot_list <- list()  # Initialize an empty list to store the plots

for (i in 1:length(df_list)) {
  plot <- ggplot(df_list[[i]], aes(x = LFC, y = fct_reorder(taxon, as.numeric(`Adj P Value`)))) +
    geom_bar(stat = "identity", position = "dodge", fill = "darkseagreen4") +
    geom_errorbarh(
      aes(xmin = LFC - Standard_error, xmax = LFC + Standard_error),
      position = position_dodge(0.9),
      height = 0.25,  # Adjust the height of the error bars
      size = 0.5,    # Adjust the size of the error bars
      color = "gray73"  # Set the color of the error bars
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
      title = paste("HMP2 IBD Associated MetaCyc Pathways -", titles[i]),
      x = "Log Fold Change",
      y = "MetaCyc Pathways",
      fill = "Comparison"
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))  # Adjust the width as needed
  
  plot_list[[i]] <- plot  # Store the plot in the list
}

plot1 <- plot_list[[1]]
plot2 <- plot_list[[2]]
plot1 
plot2

# Save plot1 as an image file
ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 IBD top 10 sig metacyc pathways.png", plot1, dpi = 800, units = "in", height = 10, width = 12)

# Save plot2 as an image file
ggsave("HMP2_Payami/Figures/Combined HMP2 IBD/HMP2 IBD top 20 sig metacyc pathways.png", plot2, dpi = 800, units = "in", height = 10, width = 12)
