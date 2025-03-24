
library(parallel)
library(bayesplot)
library(posterior)
library(tidyverse)
library(brms)
library(projpred)
library(data.table)
library(cmdstanr)


seed <- 123
set.seed(seed)

#### Part 1: Set up conditions -----
# Define the conditions grid
N <- c(50, 75)
corr <- c(0.2, 0.8)
TE <- c("clustered", "mixed")
r <- 1:100
# r <- sample(1:500, size = 100, replace = FALSE)  # Randomly sample 100 values from 1:500


conditions <- expand.grid(N = N, corr = corr, TE = TE, r = r)

#### Part 2: Set up workflow function -----

workflow_fun <- function(pos, conditions) {
  
  # Start timing
  start.time <- Sys.time()
  
  source("workflow_functions.R")
  
  # Extract condition and replication information
  sample_size <- conditions$N[pos]
  replication <- conditions$r[pos]
  corr <- conditions$corr[pos]
  TE <- conditions$TE[pos]
  # K <- conditions$K[pos]
  
  # Load and prepare training data
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
  # Standardize training data
  train_data[, -1] <- scale(train_data[, -1]) %>% as.data.frame()  # Standardize all columns except for 'y'
  
  # Standardize test data using training statistics
  test_data_scaledOnTrain <- test_data # Copy the original test data
  test_data_scaledOnTrain[, -1] <- scale(test_data[, -1], center = train_means, scale = train_sds) %>% as.data.frame() # Standardize all columns except for 'y'
  # Standardize test data using its own test statistics
  test_data[, -1] <- scale(test_data[, -1]) %>% as.data.frame() # Standardize all columns except for 'y'
  
  # Create name for saving output
  name <- paste0("N", sample_size, "_corr", corr, "_TE", TE, "_r", replication)
  
  # Fit the reference model
  refm_fit <- fit_reference_model(train_data)
  
  # Save convergence info reference model and save the reference model summary
  summ_ref <- summary(refm_fit)
  summ_ref$fixed$div <- rstan::get_num_divergent(refm_fit$fit)
  save(summ_ref, file = paste0("output_study2/ppvs/summ_ref/summ_ref_", name, ".RData"))
  
  # Run projection predictive variable selection
  ppvs_out <- run_projpred(refm_fit, K = 5, nterms_max = 16)
  
  # Save `cvvs` if any heuristic failed
  if (any(is.na(ppvs_out$suggested_sizes))) {
    cat("Some suggest_size heuristics failed for", name, ".\n")
    cvvs_out <- ppvs_out$cvvs
    save(cvvs_out, file = paste0("output_study2/ppvs/cvvs4suggest_size_failed/cvvs_", name, ".RData"))
  }
  
  # Remove `cvvs` from `ppvs_out` to save memory
  ppvs_out$cvvs <- NULL
  
  # Add selection based on projpred to reference model summary
  df_out <- summ_ref$fixed
  for (heuristic in names(ppvs_out$suggested_sizes)) {
    suggested_size <- ppvs_out$suggested_sizes[[heuristic]]
    if (is.na(suggested_size)) {
      df_out[[paste0("selected_pred_ppvs_", heuristic)]] <- NA
    } else {
      selected_pred_ppvs <- head(ppvs_out$ranking[["fulldata"]], suggested_size)
      df_out[[paste0("selected_pred_ppvs_", heuristic)]] <- 0
      df_out[[paste0("selected_pred_ppvs_", heuristic)]][rownames(df_out) %in% selected_pred_ppvs] <- 1
      df_out[[paste0("selected_pred_ppvs_", heuristic)]][rownames(df_out) == "Intercept"] <- 1
    }
  }
  
  save(ppvs_out, file = paste0("output_study2/ppvs/ppvs_out/ppvs_out_", name, ".RData"))
  save(df_out, file = paste0("output_study2/ppvs/df_out/df_out_", name, ".RData"))

  
  # Compute predictions for each heuristic
  
  # # Create an empty list to store predictions from all heuristics
  # prj_pred_results_list <- list()
  # 
  # # Store true test y-values (for easy PMSE calculation later)
  # prj_pred_results_list[["y"]] <- test_data$y 
  # 
  # for (heuristic in names(ppvs_out$suggested_sizes)) {
  #   suggested_size <- ppvs_out$suggested_sizes[[heuristic]]
  #   
  #   if (is.na(suggested_size)) {
  #     prj_pred_results_list[[heuristic]] <- NULL
  #     next
  #   }
  #   
  #   # Extract selected predictors
  #   predictors_final <- head(ppvs_out$ranking[["fulldata"]], suggested_size)
  #   
  #   # Get the reference model object
  #   refm_obj <- get_refmodel(refm_fit)
  #   
  #   # Generate projections
  #   prj <- project(
  #     refm_obj,
  #     predictor_terms = predictors_final,
  #     verbose = FALSE
  #   )
  #   prj_mat <- as.matrix(prj)
  #   prj_drws <- as_draws_matrix(prj_mat)
  #   prj_smmry <- summarize_draws(prj_drws, "mean")
  #   
  #   # Predict on test data
  #   
  #   ## Point prediction
  #   prj_point_pred <- prj_smmry$mean[1] + as.matrix(test_data[, predictors_final]) %*% prj_smmry$mean[-c(1, nrow(prj_smmry))]
  #   prj_point_pred_testScaledOnTrain <- prj_smmry$mean[1] + as.matrix(test_data_scaledOnTrain[, predictors_final]) %*% prj_smmry$mean[-c(1, nrow(prj_smmry))]
  #   
  #   ## Linear predictor with uncertainty
  #   prj_linpred <- proj_linpred(
  #     prj,
  #     newdata = test_data[, -1], # drop the dependent 'y' column
  #     integrated = TRUE
  #   )
  #   prj_linpred_testScaledOnTrain <- proj_linpred(
  #     prj,
  #     newdata = test_data_scaledOnTrain[, -1], # drop the dependent 'y' column
  #     integrated = TRUE
  #   )
  #   
  #   ## Predictive distribution for a new observation
  #   prj_pred <- proj_predict(prj, newdata = test_data[, -1], .seed = 123)
  #   prj_pred_testScaledOnTrain <- proj_predict(prj, newdata = test_data_scaledOnTrain[, -1], .seed = 123)
  #   
  #   # Store results for this heuristic inside prj_pred_results_list
  #   prj_pred_results_list[[heuristic]] <- list(
  #     point_pred = as.vector(prj_point_pred),
  #     point_pred_testScaledOnTrain = as.vector(prj_point_pred_testScaledOnTrain),
  #     linpred = as.vector(prj_linpred$pred),
  #     linpred_testScaledOnTrain = as.vector(prj_linpred_testScaledOnTrain$pred),
  #     prj_pred = prj_pred,  # 400 × 1000 posterior predictive samples
  #     prj_pred_testScaledOnTrain = prj_pred_testScaledOnTrain  # 400 × 1000 predictions scaled on training stats
  #   )
  # }
  # 
  
  # Create an empty list to store predictions from all heuristics
  prj_pred_results_list <- list()
  
  # Store true test y-values (for easy PMSE calculation later)
  prj_pred_results_list[["y"]] <- test_data$y 
  
  # Step 1: Group heuristics by suggested size
  # Convert suggested sizes to a named vector
  suggested_sizes_vec <- unlist(ppvs_out$suggested_sizes)
  
  # Remove NA values to avoid issues
  suggested_sizes_vec <- suggested_sizes_vec[!is.na(suggested_sizes_vec)]
  
  # Create a proper grouping: heuristics mapped to suggested sizes
  heuristic_groups <- split(names(suggested_sizes_vec), suggested_sizes_vec)
  
  # Step 2: Loop over unique suggested sizes
  prj_cache <- list()  # Cache projections to avoid recomputation
  
  for (suggested_size in unique(ppvs_out$suggested_sizes)) {
    
    if (is.na(suggested_size)) {
      # If a heuristic fails, store NULL for all heuristics in that group
      for (heuristic in heuristic_groups[[as.character(suggested_size)]]) {
        prj_pred_results_list[[heuristic]] <- NULL
      }
      next
    }
    
    # Extract selected predictors (same for all heuristics in this group)
    predictors_final <- head(ppvs_out$ranking[["fulldata"]], suggested_size)
    
    # Check if projection has already been computed
    predictors_key <- paste(sort(predictors_final), collapse = "_")
    
    if (!predictors_key %in% names(prj_cache)) {
      # Compute projection only once per unique predictor set
      prj_cache[[predictors_key]] <- project(
        refm_obj,
        predictor_terms = predictors_final,
        verbose = FALSE
      )
    }
    
    # Reuse the stored projection
    prj <- prj_cache[[predictors_key]]
    
    # Convert projection to draws
    prj_mat <- as.matrix(prj)
    prj_drws <- as_draws_matrix(prj_mat)
    prj_smmry <- summarize_draws(prj_drws, "mean")
    
    # Predictions
    # point prediction
    prj_point_pred <- prj_smmry$mean[1] + as.matrix(test_data[, predictors_final]) %*% prj_smmry$mean[-c(1, nrow(prj_smmry))]
    prj_point_pred_testScaledOnTrain <- prj_smmry$mean[1] + as.matrix(test_data_scaledOnTrain[, predictors_final]) %*% prj_smmry$mean[-c(1, nrow(prj_smmry))]
    
    # Linear predictor with uncertainty
    prj_linpred <- proj_linpred(prj, newdata = test_data[, -1], integrated = TRUE)
    prj_linpred_testScaledOnTrain <- proj_linpred(prj, newdata = test_data_scaledOnTrain[, -1], integrated = TRUE)
    
    # The predictive distribution for a new observation (similar to posterior_predict)
    prj_pred <- proj_predict(prj, newdata = test_data[, -1], .seed = 123)
    prj_pred_testScaledOnTrain <- proj_predict(prj, newdata = test_data_scaledOnTrain[, -1], .seed = 123)
    
    # Step 3: Assign results to all heuristics in this group
    for (heuristic in heuristic_groups[[as.character(suggested_size)]]) {
      prj_pred_results_list[[heuristic]] <- list(
        point_pred = as.vector(prj_point_pred),
        point_pred_testScaledOnTrain = as.vector(prj_point_pred_testScaledOnTrain),
        linpred = as.vector(prj_linpred$pred),
        linpred_testScaledOnTrain = as.vector(prj_linpred_testScaledOnTrain$pred),
        prj_pred = prj_pred,  # 400 × 1000 posterior predictive samples
        prj_pred_testScaledOnTrain = prj_pred_testScaledOnTrain
      )
    }
  }
  
  
  # Save predictions for all heuristics
  save(prj_pred_results_list, file = paste0("output_study2/ppvs/prj_pred/prj_pred_", name, ".RData"))
  
  # End timing
  end.time <- Sys.time()
  time.taken <- round(difftime(end.time, start.time, units = "secs"), 2)
  
  # Save timing information
  save(time.taken, file = paste0("output_study2/ppvs/time_taken/time_", name, ".RData"))
}

workflow_fun(12, conditions)
mean((prj_pred_results_list$y - prj_pred_results_list$default$point_pred)^2)

##### Part 3: Run in parallel -----


# Set the number of workers (number of parallel processes)
nworkers <- 2 # number of cores to use, max = 20
cl <- makeCluster(nworkers, type = "PSOCK")  # Create a PSOCK cluster

# Ensure all required libraries and the script are loaded on each worker
clusterEvalQ(cl, {
  options(future.globals.maxSize = 2 * 1024^3)  # Increase memory limit
  library(tidyverse)
  library(brms)
  library(projpred)
  library(data.table)
  # library(cmdstanr)
  library(posterior)
  source("workflow_functions.R")  # Load the script containing functions
})

# Export required objects to the workers
clusterExport(cl, varlist = c("conditions", "workflow_fun", "seed"),
              envir = environment())

# Run the workflow in parallel
#test_conditions <- conditions[1:4, ]

# system.time(out <- parLapplyLB(cl, 1:nrow(conditions), workflow_fun, conditions = conditions))
system.time({
  out <- parLapplyLB(cl, 1:nrow(conditions), function(i) {
    workflow_fun(i, conditions)  # Only pass row index, not the entire data frame
  })
})
# Stop the cluster after execution
stopCluster(cl)





# # Get the list of saved df_out files
# saved_files <- list.files("output/df_out/", full.names = FALSE)
# 
# # Extract the replication numbers from filenames
# saved_reps <- gsub(".*_r(\\d+)\\.RData$", "\\1", saved_files)
# saved_reps <- as.numeric(saved_reps)
# 
# 
# # Find the replications that were not completed
# setdiff(conditions$r, saved_reps)

completed_files <- list.files("output_study2/ppvs/time_taken/", pattern = "^time_N", full.names = FALSE)
library(stringr)

# Extract information from filenames
completed_conditions <- data.frame(
  N = as.numeric(str_extract(completed_files, "(?<=time_N)\\d+")),
  corr = as.numeric(str_extract(completed_files, "(?<=_corr)0\\.\\d+")),
  TE = str_extract(completed_files, "(?<=_TE)[a-z]+"),
  r = as.numeric(str_extract(completed_files, "(?<=_r)\\d+"))
)

# View completed conditions
head(completed_conditions)
nrow(completed_conditions)

# Find rows in `conditions` that are NOT in `completed_conditions`
remaining_conditions <- dplyr::anti_join(conditions, completed_conditions, by = c("N", "corr", "TE", "r"))

# View how many are left
nrow(remaining_conditions)