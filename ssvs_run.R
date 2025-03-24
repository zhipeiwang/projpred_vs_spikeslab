
library(parallel)
library(SSVS)


#### Part 1: Set up conditions -----
# Define the conditions grid
N <- c(50, 75)
corr <- c(0.2, 0.8)
TE <- c("clustered", "mixed")
r <- 1:100
# r <- sample(1:500, size = 100, replace = FALSE)  # Randomly sample 100 values from 1:500


conditions <- expand.grid(N = N, corr = corr, TE = TE, r = r)

#### Part 2: Set up workflow function -----
# Workflow function for SSVS
workflow_ssvs <- function(pos, conditions) {
  seed <- 123  # Set seed for reproducibility
  # Start timing
  start.time <- Sys.time()

  source("workflow_functions.R")
  
  # Extract condition-specific information
  sample_size <- conditions$N[pos]
  replication <- conditions$r[pos]
  corr <- conditions$corr[pos]
  TE <- conditions$TE[pos]
  
  # Load and prepare the dataset
  data <- load_and_prepare_data(
    sample_size = sample_size,
    replication = replication,
    corr = corr,
    TE = TE
  )
  
  train_data <- data$train
  test_data <- data$test
  
  # Extract original training means and standard deviations BEFORE modifying train_data
  train_means <- colMeans(train_data[, -1])
  train_sds <- apply(train_data[, -1], 2, sd)
  # ssvs() function scales the data internally, so we don't need to scale the training data
  # Standardize test data using its own test statistics
  test_data[, -1] <- scale(test_data[, -1]) %>% as.data.frame() # Standardize all columns except for 'y'
  # Standardize test data using training statistics
  test_data_scaledOnTrain <- data$test # Copy the original test data
  test_data_scaledOnTrain[, -1] <- scale(test_data_scaledOnTrain[, -1], center = train_means, scale = train_sds) %>% as.data.frame() # Standardize all columns except for 'y'
  
  # Define a unique name for saving the results
  name <- paste0("N", sample_size, "_corr", corr, "_TE", TE, "_r", replication)
  
  # Run SSVS
  set.seed(seed)  # Ensure reproducibility
  ssvs_out <- ssvs(data = train_data, x = paste0("X", 1:75), y = "y", progress = TRUE)
  
  # Get summary and convert to dataframe
  ssvs_summary <- summary(ssvs_out, interval = 0.95)
  ssvs_df_out <- as.data.frame(ssvs_summary)
  
  # Add the column `selected_pred_ssvs`
  ssvs_df_out$selected_pred_ssvs <- ifelse(ssvs_df_out$MIP > 0.5, 1, 0)
  
  # Save the dataframe to an RData file
  save(ssvs_df_out, file = paste0("output_study2/ssvs/ssvs_df_out/ssvs_", name, ".RData"))
  
  # Compute predictions
  ssvs_pred_results <- list()
  # Extract necessary matrices
  beta_draws <- ssvs_out$beta[, ssvs_df_out$selected_pred_ssvs == 1]
  intercept_draws <- ssvs_out$int
  sigma_draws <- 1 / sqrt(ssvs_out$taue)
  X_test <- as.matrix(test_data[, ssvs_df_out$Variable[ssvs_df_out$selected_pred_ssvs == 1]])
  X_test_scaledOnTrain <- as.matrix(test_data_scaledOnTrain[, ssvs_df_out$Variable[ssvs_df_out$selected_pred_ssvs == 1]])
  
  # Point predictions using posterior means
  point_pred_avgBeta <- mean(intercept_draws) + 
    X_test %*% ssvs_df_out$`Avg Beta`[ssvs_df_out$selected_pred_ssvs == 1]
  
  point_pred_avgNonZeroBeta <- mean(intercept_draws) + 
    X_test %*% ssvs_df_out$`Avg Nonzero Beta`[ssvs_df_out$selected_pred_ssvs == 1]
  
  # Point predictions using training-scaled test data
  point_pred_testScaledOnTrain_avgBeta <- mean(intercept_draws) + 
    X_test_scaledOnTrain %*% 
    ssvs_df_out$`Avg Beta`[ssvs_df_out$selected_pred_ssvs == 1]
  
  point_pred_testScaledOnTrain_avgNonZeroBeta <- mean(intercept_draws) + 
    X_test_scaledOnTrain %*% 
    ssvs_df_out$`Avg Nonzero Beta`[ssvs_df_out$selected_pred_ssvs == 1]
  
  # Linear predictor with uncertainty
  linpred <- colMeans(intercept_draws + beta_draws %*% t(X_test))
  
  linpred_testScaledOnTrain <- colMeans(intercept_draws + beta_draws %*% t(X_test_scaledOnTrain))
  
  # ---- The predictive distribution for a new observation ----
  # Number of posterior draws
  num_draws <- nrow(ssvs_out$beta)  
  
  # Compute predictive distribution using matrix multiplication (broadcasting the intercept)
  set.seed(seed)
  pred <- intercept_draws + beta_draws %*% t(X_test) + matrix(rnorm(num_draws * nrow(X_test), mean = 0, sd = sigma_draws), nrow = num_draws)
  set.seed(seed)
  pred_testScaledOnTrain <- intercept_draws + beta_draws %*% t(X_test_scaledOnTrain) + matrix(rnorm(num_draws * nrow(X_test_scaledOnTrain), mean = 0, sd = sigma_draws), nrow = num_draws)
  
  # Save the predictions to an RData file
  ssvs_pred_results <- list(
    y = test_data$y,  # True values
    point_pred_avgBeta = as.vector(point_pred_avgBeta),
    point_pred_avgNonZeroBeta = as.vector(point_pred_avgNonZeroBeta),
    point_pred_testScaledOnTrain_avgBeta = as.vector(point_pred_testScaledOnTrain_avgBeta),
    point_pred_testScaledOnTrain_avgNonZeroBeta = as.vector(point_pred_testScaledOnTrain_avgNonZeroBeta),
    linpred = linpred,
    linpred_testScaledOnTrain = linpred_testScaledOnTrain,
    pred = pred, # Full posterior predictive distribution (num_draws x nrow(X_test))
    pred_testScaledOnTrain = pred_testScaledOnTrain # Full posterior predictive distribution (num_draws x nrow(X_test_scaledOnTrain))
  )
  save(ssvs_pred_results, file = paste0("output_study2/ssvs/ssvs_pred_results/ssvs_pred_results_", name, ".RData"))
  
  # End timing
  end.time <- Sys.time()
  time.taken <- round(difftime(end.time, start.time, units = "secs"), 2)
  
  # Save timing information
  save(time.taken, file = paste0("output_study2/ssvs/time_taken/time_", name, ".RData"))
}


##### Part 3: Run in parallel -----
# Set up the cluster
nworkers <- 4 
cl <- makeCluster(nworkers, type = "PSOCK")

# Load libraries and workflow function on each worker
clusterEvalQ(cl, {
  library(tidyverse)
  library(SSVS)
  library(data.table)
  source("workflow_functions.R")
})

# Export required objects to the workers
clusterExport(cl, varlist = c("conditions", "workflow_ssvs", "seed"), envir = environment())

# Run the workflow in parallel
system.time(out <- parallel::clusterApplyLB(cl, 1:nrow(conditions), workflow_ssvs, conditions = conditions))

# Stop the cluster
stopCluster(cl)

# time elapsed: 292.06
