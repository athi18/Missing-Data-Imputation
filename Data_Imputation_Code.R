# data sctructure - Gene Expression Matrix
# Column (SRR7538444) represent sample
# Row (ENSMUSG00002075870) represent the genes

install.packages("naniar")
library(naniar)
df <- read.csv("results_DE-5%.csv")
df1 <- read.table("~/Desktop/Coursework/Missing Data/countdata5.txt", header = T, row.names = 1)
df1
str(df1)
res = mcar_test(df1)
res
res$p.value
vis_miss(df1, warn_large_data = FALSE)
gg_miss_var(df1)

# Imputation 1 using missForest package
install.packages("missForest")
library(missForest)
df1_imp <- missForest(df1)

# Mice
install.packages("mice")
library(mice)
imputed_data <-  mice(df1, method="cart")
full_data <- complete(imputed_data) 
summary(full_data)

# amelia
install.packages("Amelia")
library(Amelia)
library(dplyr)
# lets try out on the sample data first
s = sample_n(df1,200)
set.seed(123)  # For reproducibility
preprocessed_data <- amelia(s,m=1,ncpus=1,frontend=FALSE,p2s=1)

# KNN Impute
install.packages('performanceEstimation')
library(performanceEstimation)
knn_s1 = knnImp(s,k=5)
sum(is.na(knn_s))


# Regression Impute - llsImpute

BiocManager::install("pcaMethods")
library(pcaMethods)
imputed_data_reg <- llsImpute(s, 
                          k = 10, correlation = "pearson", 
                          allVariables = TRUE)
imputed_data_reg@completeObs
#Imputation using pearson correlation

# Imputation using XGBoost or MixGB
install.packages("mixgb")
library(mixgb)
params <- list(subsample = 0.7)
mixgb_data <- mixgb(data=s, m=1,
                    xgb.params = params, 
                    nrounds = 50,
                    early_stopping_rounds = 10)



# Matrix Factorization or SVD
install.packages("softImpute")
library(softImpute)
s_matrix <- as.matrix(s)
soft_data = softImpute(s_matrix, rank=35, lambda = 30, type = "svd", maxit = 30)
# compute the factorization
S_sft <- soft_data$u %*% diag(soft_data$d) %*% t(soft_data$v)
# replace missing values by computed values
S_sft[which(!is.na(s_matrix))] <- s_matrix[which(!is.na(s_matrix))]
#converting the matrix back to a dataframe
s_dataframe <- as.data.frame(S_sft)

# Preserve row names
rownames(s_dataframe) <- rownames(s)
colnames(s_dataframe) <- colnames(s)



# PCA for Data Imputation
install.packages("missMDA")
library(missMDA)
ncp_pca <- estim_ncpPCA(s,method.cv="loo")$ncp
pca_imp <- imputePCA(s, ncp = ncp_pca)
data_pca <- pca_imp$completeObs

 