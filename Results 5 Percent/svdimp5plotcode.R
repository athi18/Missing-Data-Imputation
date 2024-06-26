# Load necessary libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}
if (!requireNamespace("broom", quietly = TRUE)) {
  install.packages("broom")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

library(ggplot2)
library(dplyr)
library(ggrepel)
library(reshape2)
library(broom)
library(tidyr)

# Load the dataset
df <- read.csv("C:/Users/Sammy/Documents/Missing data/Imputations/5 percent/SVD_imputation_5_2024-06-23.csv")

# Generate random colors for genes
set.seed(123) # For reproducibility
gene_colors <- setNames(sample(colors(), nrow(df), replace = TRUE), df$Gene)

# Function to generate volcano plot for a given dataset and collect significant genes
generate_volcano_plot_and_collect_genes <- function(df, method_name, gene_colors) {
  df$id <- 1:nrow(df)
  df_long <- reshape2::melt(df, id.vars = 'id', variable.name = 'Gene', value.name = 'Expression')
  
  # Convert Expression column to numeric, coercing invalid entries to NA
  df_long$Expression <- suppressWarnings(as.numeric(as.character(df_long$Expression)))
  
  # Remove rows with NA, NaN, or Inf values
  df_long <- df_long %>% filter(!is.na(Expression) & !is.infinite(Expression) & !is.nan(Expression))
  
  # Debug print to check data
  print("Data after melting and filtering:")
  print(head(df_long))
  print(summary(df_long$Expression))
  
  # Add a 'Group' column to simulate different imputation methods (replace with actual groups if available)
  set.seed(123) # For reproducibility
  df_long$Group <- sample(c('Method1', 'Method2', 'Method3'), nrow(df_long), replace = TRUE)
  
  # Perform ANOVA for each gene
  anova_results <- df_long %>%
    group_by(Gene) %>%
    filter(!is.na(Expression)) %>%
    do(tidy(aov(Expression ~ Group, data = .)))
  
  # Debug print to check ANOVA results
  print("ANOVA results:")
  print(head(anova_results))
  
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
  
  # Debug print to check log2FC results
  print("log2FC results:")
  print(head(log2fc))
  
  # Combine results into a dataframe
  results <- data.frame(
    Gene = anova_results$Gene[anova_results$term == 'Group'],
    log2FC = log2fc$log2FC,
    PValue = p_values,
    AdjustedPValue = adjusted_p_values
  )
  
  results$Expression <- 'Not Significant'
  results$Expression[results$log2FC >= 0.15 & results$AdjustedPValue <= 0.1] <- 'Upregulated'  # Adjust threshold
  
  # Debug print to check results
  print("Results with classifications:")
  print(head(results))
  
  # Create the volcano plot
  p <- ggplot(results, aes(log2FC, -log10(AdjustedPValue), color = Gene)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = gene_colors) +
    xlab("log2FC") +
    ylab("-log10(Adjusted P-Value)") +
    geom_vline(xintercept = 0.15, linetype="dotted", color="black") +
    geom_hline(yintercept = -log10(0.1), linetype="dotted", color="black") +  # Adjust threshold
    ggtitle(paste("Volcano Plot -", method_name)) +
    theme_minimal() +
    theme(legend.position = "none", # Hide legend
          plot.title = element_text(hjust = 0.5, face="bold", size=14),
          axis.title = element_text(face="bold", size=12),
          axis.text = element_text(size=10))
  
  # Identify top significant genes for labeling
  top_genes <- results %>% arrange(AdjustedPValue) %>% head(10)
  
  # Add labels to the top genes with increased max.overlaps
  p <- p +
    geom_label_repel(data = top_genes,
                     mapping = aes(log2FC, -log10(AdjustedPValue), label = Gene),
                     size = 3, box.padding = 0.3, point.padding = 0.3, segment.color = 'grey50',
                     max.overlaps = 20)  # Increased max.overlaps
  
  # Ensure the plot is printed
  print(p)
  
  # Save the plot as a PNG file
  ggsave(paste0("Volcano_Plot_", method_name, ".png"), p, width=10, height=7, dpi=300)
  
  # Collect significant genes
  significant_genes <- results$Gene[results$AdjustedPValue <= 0.1]
  
  return(significant_genes)
}

# Generate volcano plot for the dataset and collect significant genes
method_name <- "SVD_imputation_5"
significant_genes <- generate_volcano_plot_and_collect_genes(df, method_name, gene_colors)
