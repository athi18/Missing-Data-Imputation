# Install required libraries if not already installed
required_packages <- c("UpSetR", "dplyr", "reshape2", "broom", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load necessary libraries
library(UpSetR)
library(dplyr)
library(reshape2)
library(broom)
library(tidyr)

# Define the function to classify genes with relaxed thresholds
classify_genes <- function(df, log2fc_threshold = 0.2, pvalue_threshold = 0.1) {
  df$id <- 1:nrow(df)
  df_long <- melt(df, id.vars = 'id', variable.name = 'Gene', value.name = 'Expression')
  
  # Convert Expression column to numeric, coercing invalid entries to NA
  df_long$Expression <- suppressWarnings(as.numeric(as.character(df_long$Expression)))
  
  # Remove rows with NA, NaN, or Inf values
  df_long <- df_long %>% filter(!is.na(Expression) & !is.infinite(Expression) & !is.nan(Expression))
  
  # Add a 'Group' column to simulate different imputation methods (replace with actual groups if available)
  set.seed(123) # For reproducibility
  df_long$Group <- sample(c('Method1', 'Method2', 'Method3'), nrow(df_long), replace = TRUE)
  
  # Perform ANOVA for each gene
  anova_results <- df_long %>%
    group_by(Gene) %>%
    filter(!is.na(Expression)) %>%
    do(tidy(aov(Expression ~ Group, data = .)))
  
  # Extract p-values and calculate adjusted p-values
  p_values <- anova_results %>% filter(term == 'Group') %>% pull(p.value)
  adjusted_p_values <- p.adjust(p_values, method = 'fdr')
  
  # Calculate log2 fold changes
  log2fc <- df_long %>%
    group_by(Gene, Group) %>%
    summarise(mean_expression = mean(Expression, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Group, values_from = mean_expression) %>%
    mutate(log2FC = log2(ifelse(Method1 != 0, Method2 / Method1, NA)))  # Handle division by zero
  
  # Filter out rows with NA/NaN/Inf values in the log2FC column
  log2fc <- log2fc %>% filter(!is.na(log2FC) & !is.infinite(log2FC) & !is.nan(log2FC))
  
  # Merge results to handle different row numbers
  anova_results_filtered <- anova_results %>% filter(term == 'Group') %>% select(Gene, p.value)
  log2fc_filtered <- log2fc %>% select(Gene, log2FC)
  
  # Merge by Gene to ensure consistency in row numbers
  merged_results <- merge(anova_results_filtered, log2fc_filtered, by = "Gene", all = TRUE)
  merged_results$AdjustedPValue <- adjusted_p_values[1:nrow(merged_results)]
  
  # Classify genes based on relaxed thresholds
  classified_genes <- merged_results %>%
    mutate(Classification = case_when(
      log2FC >= log2fc_threshold & p.value < pvalue_threshold ~ 'Upregulated',
      log2FC <= -log2fc_threshold & p.value < pvalue_threshold ~ 'Downregulated',
      TRUE ~ 'Not Significant'
    ))
  
  return(classified_genes)
}

# List of datasets and method names
datasets <- list(
  "C:/Users/Sammy/Documents/Missing data/Imputations/amelia_imp_5_.csv" = "amelia_imp_5",
  "C:/Users/Sammy/Documents/Missing data/Imputations/5 percent/ClustImpute_5_2024-06-23.csv" = "ClustImpute_5",
  "C:/Users/Sammy/Documents/Missing data/Imputations/KNN_imputation_5_2024-06-23.csv" = "KNN_imputation_5",
  "C:/Users/Sammy/Documents/Missing data/Imputations/5 percent/mice_imp_5_.csv" = "mice_imp_5",
  "C:/Users/Sammy/Documents/Missing data/Imputations/5 percent/missForest_imputation_5_2024-06-23.csv" = "missForest_imputation_5",
  "C:/Users/Sammy/Documents/Missing data/Imputations/PCA_imputation_2024-06-24.csv" = "PCA_imputation_5",
  "C:/Users/Sammy/Documents/Missing data/Imputations/RegLLs_imputation_5_2024-06-23.csv" = "RegLLs_imputation_5",
  "C:/Users/Sammy/Documents/Missing data/Imputations/5 percent/SVD_imputation_5_2024-06-23.csv" = "SVD_imputation_5",
  "C:/Users/Sammy/Documents/Missing data/Imputations/XGB_imp_5.csv" = "XGB_imp_5"
)

# Collect classified genes from all methods
all_classified_genes <- list()

for (file_path in names(datasets)) {
  method_name <- datasets[[file_path]]
  df <- read.csv(file_path)
  classified_genes <- classify_genes(df)
  
  # Store classified genes
  all_classified_genes[[method_name]] <- classified_genes$Gene
}

# Prepare data for the UpSet plot
upset_data <- fromList(all_classified_genes)

# Generate the UpSet plot
upset(
  upset_data,
  sets = names(all_classified_genes),
  keep.order = TRUE,
  order.by = "freq",
  main.bar.color = "blue",
  sets.bar.color = "red",
  text.scale = c(2, 1.5, 1.5, 1.5, 1.5, 1.5),
  number.angles = 30,
  point.size = 3.5,
  line.size = 2.5
)
