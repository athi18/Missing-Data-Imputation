# Gene Expression Imputation and DGE Analysis
This project evaluates different imputation methods for handling missing data in gene expression matrices and explores their impact on differential gene expression (DGE) analysis.

# Data Description

This project uses a Gene Expression Matrix with 5%, 10% and 30% missingness.

# Imputation Methods

We compared nine different imputation methods to recover missing values in the gene expression data:

- [x] XGB Impute
- [x] K Means Cluster Impute
- [x] SVD Impute
- [x] Regression Impute llsImpute
- [x] KNN Impute
- [x] MissForest
- [x] PCA Impute
- [x] Amelia 
- [x] Mice


We evaluated the performance of each imputation method using the following metrics:

Mean Squared Error (MSE)
Root Mean Squared Error (RMSE)
Mean Absolute Error (MAE)
The evaluation script calculates these metrics for each method across three datasets with varying missing data percentages (5%, 10%, 30%).

# Results

The evaluation results suggest that XGB, Regression, KNN, PCA and MissForest worked exceptionally well.
