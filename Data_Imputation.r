# data sctructure - Gene Expression Matrix
# Column (SRR7538444) represent sample
# Row (ENSMUSG00002075870) represent the genes
# DATA IMPUTATION HANDBOOK 
# Contributors - Athira Anandan, Usama Khan, Isobel Kate Donneley, Kaustubh Murali Prakash
# Group B
# For more info - https://github.com/Sirczechs-Kaustubh/Missing-Data-Imputation

######################## Data Visualization ####################

# Install the required Packages
install.packages("ggplot2")
install.packages("reshape2")
install.packages("pheatmap")
install.packages("naniar")
install.packages("missRanger")
install.packages("mice")
install.packages("visdat")

library(ggplot2)
library(reshape2)

#Load the data 
d1 <- read.csv("C:/Users/ROHIT/Downloads/countdata5.csv")
d2 <- read.csv("C:/Users/ROHIT/Downloads/countdata10.csv")
d3 <- read.csv("C:/Users/ROHIT/Downloads/countdata30.csv")

# Checking Data Structure and summary

# For 5% Dataset
str(d1)
head(d1)
summary(d1)

# For 10% Dataset
str(d2)
head(d2)
summary(d2)

# For 30% Dataset
str(d3)
head(d3)
summary(d3)

# Checking missing data patterns in each dataset
missing_values<- data.frame(
  Dataset = c("d1", "d2", "d3"),
  Missing_values <- c(sum(is.na(d1)), sum(is.na(d2)), sum(is.na(d3)))
)
# Print result
print(missing_values)

# Visualizing the missing data
# Create Label Vectors
column_labels <- paste0("Identifiers", 1:39)
row_labels <- paste0("Genome", 1:nrow(d1))

column_labels2 <- paste0("Identifiers2", 1:39)
row_labels2 <- paste0("Genome2", 1:nrow(d2))

column_labels3 <- paste0("Identifiers3", 1:39)
row_labels3 <- paste0("Genome3", 1:nrow(d3))


# Convert Data to Long Format
d1_long <- melt(as.matrix(d1))
names(d1_long) <- c("Genome", "Identifiers", "Value")

d2_long <- melt(as.matrix(d2))
names(d2_long) <- c("Genome2", "Identifiers2", "Value2")

d3_long <- melt(as.matrix(d3))
names(d3_long) <- c("Genome3", "Identifiers3", "Value")
#transforming to long format helps with easy grouping and plotting of graphs


# Add labels
d1_long$Genome <- row_labels[d1_long$Genome]
d1_long$Identifiers <- column_labels[d1_long$Identifiers]

d2_long$Genome2 <- row_labels[d2_long$Genome2]
d2_long$Identifiers2 <- column_labels[d2_long$Identifiers2]

d3_long$Genome3 <- row_labels[d3_long$Genome3]
d3_long$Identifiers3 <- column_labels[d3_long$Identifiers3]


# Plotting using ggplot2

#For 5% Dataset
ggplot(d1_long, aes(x = Genome, y = Value, color = Identifiers)) + geom_point() +labs(x = "Genome", y = "Values", color = "Identifiers") +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_x_discrete(labels=row_labels)

#For 10% Dataset
ggplot(d2_long, aes(x = Genome2, y = Value, color = Identifiers2)) + geom_point() +labs(x = "Genome", y = "Values", color = "Identifiers") +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_x_discrete(labels=row_labels)
#For 30% Dataset
ggplot(d3_long, aes(x = Genome3, y = Value, color = Identifiers3)) + geom_point() +labs(x = "Genome", y = "Values", color = "Identifiers") +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_x_discrete(labels=row_labels)

#Conclusion : Here Genome is used to colour the points based on columns.
# Histogram with density line
#For 5% Dataset
ggplot(d1_long, aes(x = Genome)) + geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "yellow", color = "black", alpha = 0.5) +  geom_density(color = "red", linewidth = 1) + ggtitle("Histogram of Values with Density Line of 5% Dataset") + xlab("Genome") + ylab("Density")
#For 10% Dataset
ggplot(d2_long, aes(x = Genome2)) + geom_histogram(aes(y = after_stat(density2)), bins = 100, fill = "yellow", color = "black", alpha = 0.5) +  geom_density(color = "red", linewidth = 1) + ggtitle("Histogram of Values with Density Line of 10% Dataset") + xlab("Genome") + ylab("Density")
#For 30% Dataset
ggplot(d3_long, aes(x = Genome3)) + geom_histogram(aes(y = after_stat(density3)), bins = 100, fill = "yellow", color = "black", alpha = 0.5) +  geom_density(color = "red", linewidth = 1) + ggtitle("Histogram of Values with Density Line of 30% Dataset") + xlab("Genome") + ylab("Density")

# Heatmap to check Correlation between missing data
library(pheatmap)

# Function to convert all columns to numeric, setting non-convertible values to NA
convert_to_numeric <- function(df) {
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  return(df)
}

# Apply the conversion function to each dataset
d1_num <- convert_to_numeric(d1)
d2_num <- convert_to_numeric(d2)
d3_num <- convert_to_numeric(d3)

# Calculate the correlation matrices for each dataset
d1_cor <- cor(d1_num, use = "pairwise.complete.obs")
d2_cor <- cor(d2_num, use = "pairwise.complete.obs")
d3_cor <- cor(d3_num, use = "pairwise.complete.obs")

# Replace NA/NaN/Inf values with zeros
d1_cor[is.na(d1_cor)] <- 0
d2_cor[is.na(d2_cor)] <- 0
d3_cor[is.na(d3_cor)] <- 0

# Plot heatmaps for each correlation matrix with additional customization
pheatmap(d1_cor, main = "Correlation Heatmap for Dataset d1",
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean")

pheatmap(d2_cor, main = "Correlation Heatmap for Dataset d2",
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean")

pheatmap(d3_cor, main = "Correlation Heatmap for Dataset d3",
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean")

# Conclusion:
# The heatmaps show that the correlation between variables varies across the three datasets.
# In general, the variables in d1 have the highest correlation, followed by d2 and then d3.
# This suggests that the missing data in each dataset may be affecting the relationships between the variables.
# Checking Data structure and summary for long format
# For 5% Dataset
str(d1_long)
summary(d1_long)
# For 10% Dataset
str(d2_long)
summary(d2_long)
# For 30% Dataset
str(d3_long)
summary(d3_long)           
# Package for analyzing missing data
library(visdat)
# Visual summary of missing data
vis_miss(d1, warn_large_data = FALSE)
vis_miss(d2, warn_large_data = FALSE)
vis_miss(d3, warn_large_data = FALSE)
# Visualizing the Missing Pattern
library(naniar)
# Number of missing values in each variable
#For 5% Dataset
gg_miss_var(d1)
#For 10% Dataset
gg_miss_var(d2)
#For 30% Dataset
gg_miss_var(d3)
# Visualizing using Upset plot
#For 30% Dataset
gg_miss_upset(d1)
#For 30% Dataset
gg_miss_upset(d2)
#For 30% Dataset
gg_miss_upset(d3)
#Conclusion:
#Upset plots provide overview of intersections among the sets.
#High intesection size shows overlap of sets while moderate and low intersections show specific relationships and unique elements.
# Package for exploring missing data pattern
library(mice)
# Exploratory Data Analysis
md.pattern(d1)
md.pattern(d2)
md.pattern(d3)
# Packgae for evaluating the Missing Pattern
library(missRanger)
# Little's MCAR (Missing Completely At Random) test
#For 5% Dataset
mcar_test(d1)
#For 30% Dataset
mcar_test(d2)
#For 30% Dataset
mcar_test(d3)
# Conclusion:
# If the p-value is less than 0.05, we reject the null hypothesis and conclude that the data is not missing completely at random.
# If the p-value is greater than or equal to 0.05, we fail to reject the null hypothesis and cannot conclude that the data is not missing completely at random.

# Based on the results of the MCAR tests, we can conclude that:

# - The data in d1 is not missing completely at random (p-value < 0.05).
# - The data in d2 is missing completely at random (p-value > 0.05).
# - The data in d3 is not missing completely at random (p-value < 0.05).

######################## Data Imputation #######################
# Loading Dataset
data1 <- read.table("~/Missing Data/countdata5.txt", header = T, row.names = 1)
data2 <- read.table("~/Missing Data/countdata5.txt", header = T, row.names = 1)
data3 <- read.table("~/Missing Data/countdata5.txt", header = T, row.names = 1)


# # Extract complete cases
# complete_data <- data1[complete.cases(data), ]
# sampled_data <- sample_n(complete_data,500)
# write.csv(sampled_data, "complete_simulated.csv", row.names = TRUE)
# 
# 
# 
# # Assuming that all the missing data comes from the same dataset but with different missing percentage
# data1 <- ampute(sampled_data,prop=0.05,mech = "MCAR")
# data2 <- ampute(sampled_data,prop=0.1,mech = "MCAR")
# data3 <- ampute(sampled_data,prop=0.35,mech = "MCAR")
# 
# write.csv(df1$amp, "data_5pc.csv", row.names = TRUE)
# write.csv(df2$amp, "data_10pc.csv", row.names = TRUE)
# write.csv(df3$amp, "data_30pc.csv", row.names = TRUE)
# 
# data1_amp <- data1$amp
# data2_amp <- data2$amp
# data3_amp <- data3$amp

############## Various Imputation Techniques ################
# 1) Mice
library(mice)
imputed_data <- mice(data1, method = "norm")
imputed_dataset_1 <- complete(imputed_data, 1)
file_name <- paste("mice_imp_5_sim",".csv", sep = "")
write.csv(imputed_dataset_1, file_name)


imputed_data <- mice(data2, method = "norm")
imputed_dataset_1 <- complete(imputed_data, 1)
file_name <- paste("mice_imp_15_sim",".csv", sep = "")
write.csv(imputed_dataset_1, file_name)


imputed_data <- mice(data3, method = "norm")
imputed_dataset_1 <- complete(imputed_data, 1)
file_name <- paste("mice_imp_30_sim",".csv", sep = "")
write.csv(imputed_dataset_1, file_name)

# 2) PCA Impute
library(missMDA) # Loading Library
pca_imp <- imputePCA(data1, ncp = 1)
data_imputed <- pca_imp$completeObs
file_name <- paste("pca_imp_5_sim",".csv", sep = "")
write.csv(data_imputed, file_name)

pca_imp <- imputePCA(data2, ncp = 1)
data_imputed <- pca_imp$completeObs
file_name <- paste("pca_imp_10_sim",".csv", sep = "")
write.csv(data_imputed, file_name)

pca_imp <- imputePCA(data3, ncp = 1)
data_imputed <- pca_imp$completeObs
file_name <- paste("pca_imp_30_sim",".csv", sep = "")
write.csv(data_imputed, file_name)

# 3) Amelia
install.packages("Amelia")
library(Amelia)
file_name <- paste("amelia_imp_",".csv", sep = "5")
imputed_data <- amelia(data1,m=5,ncpus=1,frontend=FALSE,p2s=1)
write.csv(imputed_data[["imputations"]][["imp1"]],file=file_name)

file_name <- paste("amelia_imp_",".csv", sep = "10")
imputed_data <- amelia(data2,m=5,ncpus=1,frontend=FALSE,p2s=1)
write.csv(imputed_data[["imputations"]][["imp1"]],file=file_name)

file_name <- paste("amelia_imp_",".csv", sep = "30")
imputed_data <- amelia(data3,m=5,ncpus=1,frontend=FALSE,p2s=1)
write.csv(imputed_data[["imputations"]][["imp1"]],file=file_name)

# 4) missForest
library(missForest)
file_name <- paste("missforest_imp_5_",".csv", sep = "")
imputed_data <- missForest(data1)
write.csv(imputed_data$ximp, file=file_name)

file_name <- paste("missforest_imp_10_",".csv", sep = "")
imputed_data <- missForest(data2)
write.csv(imputed_data$ximp, file=file_name)

file_name <- paste("missforest_imp_30_",".csv", sep = "")
imputed_data <- missForest(data3)
write.csv(imputed_data$ximp, file=file_name)

# 5) KNN Impute
library(performanceEstimation)
file_name <- paste("KNN_imp_",".csv", sep = "5")
imputed_data <- knnImp(data1,k=5)
write.csv(imputed_data,file=file_name)

file_name <- paste("KNN_imp_",".csv", sep = "10")
imputed_data <- knnImp(data2,k=5)
write.csv(imputed_data,file=file_name)

file_name <- paste("KNN_imp_",".csv", sep = "30")
imputed_data <- knnImp(data3,k=5)
write.csv(imputed_data,file=file_name)

# 6) Regression Imp
BiocManager::install("pcaMethods")
library(pcaMethods)
file_name <- paste("Reg_imp_",".csv", sep = "5")
imputed_data <- llsImpute(data1, 
                          k = 10, correlation = "pearson", 
                          allVariables = TRUE)
write.csv(imputed_data@completeObs,file_name)


file_name <- paste("Reg_imp_",".csv", sep = "10")
imputed_data <- llsImpute(data2, 
                          k = 10, correlation = "pearson", 
                          allVariables = TRUE)
write.csv(imputed_data@completeObs,file_name)


file_name <- paste("Reg_imp_",".csv", sep = "30")
imputed_data <- llsImpute(data3, 
                          k = 10, correlation = "pearson", 
                          allVariables = TRUE)
write.csv(imputed_data@completeObs,file_name)

# 7) SVD Impute
library(softImpute)
file_name <- paste("SVD_imp_",".csv", sep = "5")
s_matrix <- as.matrix(data1)
soft_data = softImpute(s_matrix, rank=35, lambda = 30, type = "svd", maxit = 30)
# compute the factorization
S_sft <- soft_data$u %*% diag(soft_data$d) %*% t(soft_data$v)
# replace missing values by computed values
S_sft[which(!is.na(s_matrix))] <- s_matrix[which(!is.na(s_matrix))]
#converting the matrix back to a dataframe
s_dataframe <- as.data.frame(S_sft)
rownames(s_dataframe) <- rownames(data1)
colnames(s_dataframe) <- colnames(data1)
write.csv(s_dataframe,file_name)


file_name <- paste("SVD_imp_",".csv", sep = "10")
s_matrix <- as.matrix(data2)
soft_data = softImpute(s_matrix, rank=35, lambda = 30, type = "svd", maxit = 30)
# compute the factorization
S_sft <- soft_data$u %*% diag(soft_data$d) %*% t(soft_data$v)
# replace missing values by computed values
S_sft[which(!is.na(s_matrix))] <- s_matrix[which(!is.na(s_matrix))]
#converting the matrix back to a dataframe
s_dataframe <- as.data.frame(S_sft)
rownames(s_dataframe) <- rownames(data2)
colnames(s_dataframe) <- colnames(data2)
write.csv(s_dataframe,file_name)


file_name <- paste("SVD_imp_",".csv", sep = "30")
s_matrix <- as.matrix(data3)
soft_data = softImpute(s_matrix, rank=35, lambda = 30, type = "svd", maxit = 30)
# compute the factorization
S_sft <- soft_data$u %*% diag(soft_data$d) %*% t(soft_data$v)
# replace missing values by computed values
S_sft[which(!is.na(s_matrix))] <- s_matrix[which(!is.na(s_matrix))]
#converting the matrix back to a dataframe
s_dataframe <- as.data.frame(S_sft)
rownames(s_dataframe) <- rownames(data3)
colnames(s_dataframe) <- colnames(data3)
write.csv(s_dataframe,file_name)


# 8) K Means Cluster Impute

file_name <- paste("Cluster_imp_",".csv", sep = "5")
# Use to find optimal cluster
#results <- MCS(data.matrix(data), nc = 38, method1 = "kmeans", 
#method2 = "ward.D2", index = "rand")

# Impute missing values using ClustImpute with cluster information
data_imputed <- ClustImpute(data1,nr_cluster = 7, nr_iter = 10)
write.csv(data_imputed$complete_data,file_name)


file_name <- paste("Cluster_imp_",".csv", sep = "10")
# Impute missing values using ClustImpute with cluster information
data_imputed <- ClustImpute(data2,nr_cluster = 7, nr_iter = 10)
write.csv(data_imputed$complete_data,file_name)



file_name <- paste("Cluster_imp_",".csv", sep = "30")
#results <- MCS(data.matrix(data), nc = 38, method1 = "kmeans", 
#method2 = "ward.D2", index = "rand")

# Impute missing values using ClustImpute with cluster information
data_imputed <- ClustImpute(data3,nr_cluster = 7, nr_iter = 10)
write.csv(data_imputed$complete_data,file_name)


# 9) XGB Impute
library(mixgb)
file_name <- paste("XGB_imp_",".csv", sep = "5")
params <- list(subsample = 0.7)
imputed_data <- mixgb(data=data1, m=1,
                      xgb.params = params, 
                      nrounds = 100,
                      early_stopping_rounds = 50)
write.csv(imputed_data[[1]],file_name)


file_name <- paste("XGB_imp_",".csv", sep = "10")
params <- list(subsample = 0.7)
imputed_data <- mixgb(data=data2, m=1,
                      xgb.params = params, 
                      nrounds = 100,
                      early_stopping_rounds = 50)
write.csv(imputed_data[[1]],file_name)


file_name <- paste("XGB_imp_",".csv", sep = "30")
params <- list(subsample = 0.7)
imputed_data <- mixgb(data=data3, m=1,
                      xgb.params = params, 
                      nrounds = 100,
                      early_stopping_rounds = 50)
write.csv(imputed_data[[1]],file_name)


############# Metrics Calculations ######################


# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)

# Define file paths
original_data_path <- '~/Missing Data/Simulated Data/complete_simulated.csv'
imputed_data_dir <- '~/Missing Data/Simulated Data/Imputed Data/'

# Read the original dataset
original <- read.csv(original_data_path, row.names = 1, header = TRUE)

# Get the list of imputed data files
imputed_files <- list.files(path = imputed_data_dir, full.names = TRUE)

# Read datasets with varying levels of induced missingness
miss_5 <- read.csv("~/Missing Data/Simulated Data/data_5pc.csv", row.names = 1, header = TRUE)
miss_10 <- read.csv("~/Missing Data/Simulated Data/data_10pc.csv", row.names = 1, header = TRUE)
miss_30 <- read.csv("~/Missing Data/Simulated Data/data_30pc.csv", row.names = 1, header = TRUE)

# Find indices of missing values
miss_index_5 <- which(is.na(miss_5), arr.ind = TRUE)
miss_index_10 <- which(is.na(miss_10), arr.ind = TRUE)
miss_index_30 <- which(is.na(miss_30), arr.ind = TRUE)

# Function to evaluate imputation
evaluate_imputation <- function(original, imputed, missing_index) {
  mae <- mean(abs(original[missing_index] - imputed[missing_index]))
  mse <- mean((original[missing_index] - imputed[missing_index])^2)
  rmse <- sqrt(mse)
  data_range <- max(original[missing_index], na.rm = TRUE) - min(original[missing_index], na.rm = TRUE)
  nrmse <- rmse / data_range
  return(data.frame(MAE = mae, MSE = mse, RMSE = rmse, NRMSE = nrmse))
}

# Initialize lists to store results
results_5 <- tibble()
results_10 <- tibble()
results_30 <- tibble()

# Process each imputed file
for (i in imputed_files) {
  imputed <- read.csv(i, row.names = 1, header = TRUE)
  
  filename <- basename(i)
  method <- str_extract(filename, "^[^_]+")
  missingness <- as.numeric(str_extract(filename, "\\d+"))
  
  if (missingness == 5) {
    metrics <- evaluate_imputation(original, imputed, miss_index_5)
  } else if (missingness == 10) {
    metrics <- evaluate_imputation(original, imputed, miss_index_10)
  } else if (missingness == 30) {
    metrics <- evaluate_imputation(original, imputed, miss_index_30)
  }
  
  metrics$Method <- method
  metrics$Missingness <- missingness
  
  if (missingness == 5) {
    results_5 <- bind_rows(results_5, metrics)
  } else if (missingness == 10) {
    results_10 <- bind_rows(results_10, metrics)
  } else if (missingness == 30) {
    results_30 <- bind_rows(results_30, metrics)
  }
}

# Write results to CSV files
write_csv(results_5, "imputation_evaluation_results_5_percent.csv")
write_csv(results_10, "imputation_evaluation_results_10_percent.csv")
write_csv(results_30, "imputation_evaluation_results_30_percent.csv")

library(readr)
library(dplyr)
library(tidyr)
results_5 <- read_csv("imputation_evaluation_results_5_percent.csv")
results_10 <- read_csv("imputation_evaluation_results_10_percent.csv")
results_30 <- read_csv("imputation_evaluation_results_30_percent.csv")


results_long <- pivot_longer(results_30, cols = c(MAE, MSE, RMSE, NRMSE), names_to = "Metric", values_to = "Value")

# Plotting
library(ggplot2)
library(RColorBrewer)

# Ensure your 'results_long' dataset is loaded and ready to plot
# Your existing plot code with added color scale for 'Set1' palette
plot <- ggplot(results_long, aes(x = Method, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Metric, scales = "free") +
  labs(title = "Performance of Imputation Methods at 30% Missingness",
       x = "Method",
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")  # Applying Set1 palette for the fill color of bars

# Print the plot to view it
print(plot)

ggsave("30_Percent_Metrics.png", units="in", width=8, height=6, dpi=500)

############################# Result Plots ##########################
############## Plots for 5 % ###############

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

# Please change the File Paths
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
############## Plots for 10 % ###############

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
  "C:/Users/Sammy/Documents/Missing data/Imputations/10 percent/amelia_imp_10_.csv" = "amelia_imp_10",
  "C:/Users/Sammy/Documents/Missing data/Imputations/10 percent/ClustImpute_10_2024-06-24.csv" = "ClustImpute_10",
  "C:/Users/Sammy/Documents/Missing data/Imputations/10 percent/KNN_imp_10_.csv" = "KNN_imp_10",
  "C:/Users/Sammy/Documents/Missing data/Imputations/10 percent/mice_imp_10_2024-06-24.csv" = "mice_imp_10",
  "C:/Users/Sammy/Documents/Missing data/Imputations/10 percent/missforest_imp_10_.csv" = "missforest_imp_10",
  "C:/Users/Sammy/Documents/Missing data/Imputations/10 percent/PCA_imputation_10_2024-06-24.csv" = "PCA_imputation_10",
  "C:/Users/Sammy/Documents/Missing data/Imputations/10 percent/RegLLs_imputation_10_2024-06-24.csv" = "RegLLs_imputation_10",
  "C:/Users/Sammy/Documents/Missing data/Imputations/10 percent/SVD_imp_10_2024-06-24.csv" = "SVD_imputation_10",
  "C:/Users/Sammy/Documents/Missing data/Imputations/10 percent/XGB_imp_10.csv" = "XGB_imp_10"
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
######################### Plots for 30 % ###############################
# Load necessary libraries
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
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/XGB_imp_30.csv" = "XGB_imp_30",
  "C:/Users/Sammy/Documents/Missing data/Imputations/30 percent/Knn_Impute_30.csv" = "Knn_Impute_30"
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

