# Load necessary libraries
required_packages <- c("ggplot2", "dplyr", "ggrepel", "reshape2", "broom", "tidyr", "UpSetR", "RColorBrewer")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(dplyr)
library(ggrepel)
library(reshape2)
library(broom)
library(tidyr)
library(UpSetR)
library(RColorBrewer)

# Define the function to generate the global color scheme for genes
generate_gene_colors <- function(genes) {
  unique_genes <- unique(genes)
  gene_colors <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique_genes)), unique_genes)
  return(gene_colors)
}

# List of datasets and method names
datasets <- list(
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/amelia_imp_30_.csv" = "amelia_imp_30",
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/ClustImpute_30_2024-06-24.csv" = "ClustImpute_30",
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/mice_imp_30_2024-06-24.csv" = "mice_imp_30",
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/missForest_imputation_5_2024-06-23.csv" = "missForest_imputation_30",
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/PCA_imputation_30_2024-06-24.csv" = "PCA_imputation_30",
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/RegLLs_imputation_30_2024-06-24.csv" = "RegLLs_imputation_30",
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/SVD_imp_30_2024-06-24.csv" = "SVD_imputation_30",
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/XGB_imp_30.csv" = "XGB_imp_30"
)

# Collect all unique genes from all datasets
all_genes <- unique(unlist(lapply(names(datasets), function(file_path) {
  df <- read.csv(file_path)
  return(names(df))
})))

# Generate a consistent color scheme for all genes
gene_colors <- generate_gene_colors(all_genes)

# Function to generate volcano plot for a given dataset and collect significant genes
generate_volcano_plot_and_collect_genes <- function(df, method_name, gene_colors, log2fc_threshold = 0.05, pvalue_threshold = 0.75) {
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
    mutate(log2FC = log2(ifelse(!is.na(Method1) & Method1 != 0, Method2 / Method1, NA)))  # Handle division by zero and NA
  
  # Filter out rows with NA/NaN/Inf values in the log2FC column
  log2fc <- log2fc %>% filter(!is.na(log2FC) & !is.infinite(log2FC) & !is.nan(log2FC))
  
  # Merge results by Gene to handle different row numbers
  merged_results <- anova_results %>%
    filter(term == 'Group') %>%
    select(Gene, p.value) %>%
    left_join(log2fc %>% select(Gene, log2FC), by = "Gene")
  
  merged_results$AdjustedPValue <- adjusted_p_values[1:nrow(merged_results)]
  
  results <- merged_results
  
  results$Expression <- 'Not Significant'
  results$Expression[results$log2FC >= log2fc_threshold & results$AdjustedPValue <= pvalue_threshold] <- 'Upregulated'
  
  # Create the volcano plot
  p <- ggplot(results, aes(log2FC, -log10(AdjustedPValue), color = Gene)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = gene_colors) +
    xlab("log2FC") +
    ylab("-log10(Adjusted P-Value)") +
    geom_vline(xintercept = log2fc_threshold, linetype="dotted", color="black") +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype="dotted", color="black") +
    ggtitle(paste("Volcano Plot -", method_name)) +
    theme_minimal() +
    theme(legend.position = "none", # Hide legend
          plot.title = element_text(hjust = 0.5, face="bold", size=14),
          axis.title = element_text(face="bold", size=12),
          axis.text = element_text(size=10))
  
  # Identify significant genes for labeling
  significant_genes <- results %>% filter(Expression == 'Upregulated')
  
  # Add labels to the significant genes
  p <- p +
    geom_label_repel(data = significant_genes,
                     mapping = aes(log2FC, -log10(AdjustedPValue), label = Gene),
                     size = 3, box.padding = 0.3, point.padding = 0.3, segment.color = 'grey50',
                     max.overlaps = 20)  # Increased max.overlaps
  
  # Ensure the plot is printed
  print(p)
  
  # Save the plot as a PNG file
  ggsave(paste0("Volcano_Plot_", method_name, ".png"), p, width=10, height=7, dpi=300)
  
  # Collect significant genes
  significant_genes <- results$Gene[results$AdjustedPValue <= pvalue_threshold]
  
  return(significant_genes)
}

# Collect classified genes from all methods
all_classified_genes <- list()

for (file_path in names(datasets)) {
  method_name <- datasets[[file_path]]
  df <- read.csv(file_path)
  classified_genes <- generate_volcano_plot_and_collect_genes(df, method_name, gene_colors)
  
  # Store classified genes
  all_classified_genes[[method_name]] <- classified_genes
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
