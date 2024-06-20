########### Utility Function ##########
# Function to install and load packages
install_and_load <- function() {
  packages <- c("missForest", "mice", "Amelia", "dplyr", 
                "performanceEstimation", "pcaMethods", 
                "mixgb", "softImpute", "missMDA", "MCSim", "ClustImpute")
  
  # Check for missing packages
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  
  # Install missing packages
  if (length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)
  }
  
  # Load all packages
  sapply(packages, require, character.only = TRUE)
}

################# Imputation Function ###################
multi_impute <- function(data,file_ref){
  # 1) Miss Forest Imputation
  file_name <- paste("missforest_imp_",".csv", sep = file_ref)
  imputed_data <- missForest(data,maxiter = 10,ntree = 250)
  write.csv(imputed_data$ximp, file=file_name)
  
  # 2) Mice - Cart Imputation
  file_name <- paste("mice_imp_",".csv", sep = file_ref)
  imputed_data <-  mice(data, method="cart")
  completed_dataset <- complete(imputed_data)
  write.csv(completed_dataset,file=file_name)
  
  # 3) Amelia
  file_name <- paste("amelia_imp_",".csv", sep = file_ref)
  imputed_data <- amelia(data,m=5,ncpus=1,frontend=FALSE,p2s=1)
  write.csv(imputed_data[["imputations"]][["imp1"]],file=file_name)
  
  # 4) KNN Impute
  file_name <- paste("KNN_imp_",".csv", sep = file_ref)
  imputed_data <- knnImp(data,k=5)
  write.csv(imputed_data,file=file_name)
  
  # 5) Regression Impute - llsImpute
  file_name <- paste("Reg_imp_",".csv", sep = file_ref)
  imputed_data <- llsImpute(data, 
                                k = 10, correlation = "pearson", 
                                allVariables = TRUE)
  write.csv(imputed_data@completeObs,file_name)
  
  # 6) XGB Impute
  file_name <- paste("XGB_imp_",".csv", sep = file_ref)
  params <- list(subsample = 0.7)
  imputed_data <- mixgb(data=data, m=1,
                      xgb.params = params, 
                      nrounds = 50,
                      early_stopping_rounds = 10)
  write.csv(imputed_data[[1]],file_name)
  
  # 7) SVD Impute
  file_name <- paste("SVD_imp_",".csv", sep = file_ref)
  s_matrix <- as.matrix(data)
  soft_data = softImpute(s_matrix, rank=35, lambda = 30, type = "svd", maxit = 30)
      # compute the factorization
  S_sft <- soft_data$u %*% diag(soft_data$d) %*% t(soft_data$v)
      # replace missing values by computed values
  S_sft[which(!is.na(s_matrix))] <- s_matrix[which(!is.na(s_matrix))]
      #converting the matrix back to a dataframe
  s_dataframe <- as.data.frame(S_sft)
  
      # Preserve row names
  rownames(s_dataframe) <- rownames(data)
  colnames(s_dataframe) <- colnames(data)
  write.csv(s_dataframe,file_name)
  
  # 8) PCA Impute
  file_name <- paste("PCA_imp_",".csv", sep = file_ref)
  #ncp_pca <- estim_ncpPCA(s,method.cv="loo")$ncp
  pca_imp <- imputePCA(data, ncp = 1)
  data_imputed <- pca_imp$completeObs
  write.csv(data_imputed,file_name)
  
  # 9) K Means Cluster Impute
  file_name <- paste("Cluster_imp_",".csv", sep = file_ref)
  results <- MCS(data.matrix(data), nc = 38, method1 = "kmeans", 
                 method2 = "ward.D2", index = "rand")
  
        # Impute missing values using ClustImpute with cluster information
  data_imputed <- ClustImpute(data,nr_cluster = 7, nr_iter = 10)
  write.csv(data_imputed$complete_data,file_name)
}

################# Main #####################

data1 <- read.table("~/Desktop/Coursework/Missing Data/countdata5.txt", header = T, row.names = 1)
data2 <- read.table("~/Desktop/Coursework/Missing Data/countdata10.txt", header = T, row.names = 1)
data3 <- read.table("~/Desktop/Coursework/Missing Data/countdata30.txt", header = T, row.names = 1)

# Loads and Installs Necessary packages
install_and_load()
# Imputes 5% Missing Dataset
multi_impute(data=data1,s="5")
# Imputes 10% Missing Dataset
multi_impute(data=data2,s="10")
# Imputes 30% Missing Dataset
multi_impute(data=data3,s="30")

