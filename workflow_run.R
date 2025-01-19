
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
N <- c(100)
corr <- c(0.8)
loc_TE <- c("first10", "mixed")
# r <- 1:100
r <- sample(1:500, size = 100, replace = FALSE)  # Randomly sample 100 values from 1:500


conditions <- expand.grid(N = N, corr = corr, loc_TE = loc_TE, r = r)
conditions$K <- ifelse(conditions$N == 100, 5, 2)
conditions <- conditions %>%
  filter(!(corr == 0 & loc_TE == "mixed"))

#### Part 2: Set up workflow function -----
# workflow_fun <- function(pos, conditions){
#   
#   source("workflow_functions.R")
#   
#   # extract condition and replication information
#   sample_size <- conditions$N[pos]
#   replication <- conditions$r[pos]
#   corr <- conditions$corr[pos]
#   loc_TE <- conditions$loc_TE[pos]
#   K <- conditions$K[pos]
#   
#   # load and prepare data for the current condition and replication
#   data <- load_and_prepare_data(file_path = NULL, # File path is constructed dynamically
#                                 sample_size = sample_size,
#                                 replication = replication,
#                                 corr = corr,
#                                 loc_TE = loc_TE)
#   
#   
#   # create name for saving output
#   name <- paste0("N", sample_size, "_corr", corr, "_loc_TE", loc_TE, "_r", replication)
#   
#   # fit the reference model
#   refm_fit <- fit_reference_model(data)
#   
#   # save convergence info reference model and save the reference model summary
#   summ_ref <- summary(refm_fit)
#   summ_ref$fixed$div <- rstan::get_num_divergent(refm_fit$fit)
#   save(summ_ref, file = paste0("output/summ_ref/summ_ref_", name, ".RData"))
#   
#   # run projection predictive variable selection
#   ppvs_out <- run_projpred(refm_fit, K = K, nterms_max = 11)
#   
#   # Save `cvvs` if any heuristic failed
#   if (any(is.na(ppvs_out$suggested_sizes))) {
#     cat("Some suggest_size heuristics failed for", name, ".\n")
#     cvvs_out <- ppvs_out$cvvs
#     save(cvvs_out, file = paste0("output/cvvs4suggest_size_failed/cvvs_", name, ".RData"))
#   }
#   
#   # If suggest_size succeeded, remove the cvvs object from ppvs_out
#   ppvs_out$cvvs <- NULL
#   
#   # add selection based on projpred to reference model summary
#   # Proceed to create `df_out` even if some heuristics failed
#   df_out <- summ_ref$fixed
#   for (heuristic in names(ppvs_out$suggested_sizes)) {
#     suggested_size <- ppvs_out$suggested_sizes[[heuristic]]
#     if (is.na(suggested_size)) {
#       # If `suggest_size` failed, add a column with all `NA`
#       df_out[[paste0("selected_pred_ppvs_", heuristic)]] <- NA
#     } else {
#       # If `suggest_size` succeeded, populate the column
#       selected_pred_ppvs <- head(ppvs_out$ranking[["fulldata"]], suggested_size)
#       df_out[[paste0("selected_pred_ppvs_", heuristic)]] <- 0
#       df_out[[paste0("selected_pred_ppvs_", heuristic)]][rownames(df_out) %in% selected_pred_ppvs] <- 1
#     }
#     # Ensure the intercept is always selected
#     df_out[[paste0("selected_pred_ppvs_", heuristic)]][rownames(df_out) == "Intercept"] <- 1
#   }
#   
#   # Remove `cvvs` from `ppvs_out` to save memory
#   ppvs_out$cvvs <- NULL
#   
#   # save outputs
#   save(ppvs_out, file = paste0("output/ppvs_out/ppvs_out_", name, ".RData"))
#   save(df_out, file = paste0("output/df_out/df_out_", name, ".RData"))
# }





















workflow_fun <- function(pos, conditions) {
  
  # Start timing
  start.time <- Sys.time()
  
  source("workflow_functions.R")
  
  # Extract condition and replication information
  sample_size <- conditions$N[pos]
  replication <- conditions$r[pos]
  corr <- conditions$corr[pos]
  loc_TE <- conditions$loc_TE[pos]
  K <- conditions$K[pos]
  
  # Load and prepare training data
  data_train <- load_and_prepare_data(
    file_path = NULL,
    sample_size = sample_size,
    replication = replication,
    corr = corr,
    loc_TE = loc_TE
  )
  
  # Create name for saving output
  name <- paste0("N", sample_size, "_corr", corr, "_loc_TE", loc_TE, "_r", replication)
  
  # Fit the reference model
  refm_fit <- fit_reference_model(data_train)
  
  # Save convergence info reference model and save the reference model summary
  summ_ref <- summary(refm_fit)
  summ_ref$fixed$div <- rstan::get_num_divergent(refm_fit$fit)
  save(summ_ref, file = paste0("output/summ_ref/summ_ref_", name, ".RData"))
  
  # Run projection predictive variable selection
  ppvs_out <- run_projpred(refm_fit, K = K, nterms_max = 11)
  
  # Save `cvvs` if any heuristic failed
  if (any(is.na(ppvs_out$suggested_sizes))) {
    cat("Some suggest_size heuristics failed for", name, ".\n")
    cvvs_out <- ppvs_out$cvvs
    save(cvvs_out, file = paste0("output/cvvs4suggest_size_failed/cvvs_", name, ".RData"))
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
  
  save(ppvs_out, file = paste0("output/ppvs_out/ppvs_out_", name, ".RData"))
  save(df_out, file = paste0("output/df_out/df_out_", name, ".RData"))
  
  # Load test data for predictions
  data_test <- load_and_prepare_data(
    file_path = NULL,
    sample_size = 400,
    replication = replication,
    corr = corr,
    loc_TE = loc_TE
  )
  
  # Compute predictions for each heuristic
  prj_linpred_results <- data.frame(y = data_test$y)  # Start with the true y values
  refm_obj <- get_refmodel(refm_fit)
  
  for (heuristic in names(ppvs_out$suggested_sizes)) {
    suggested_size <- ppvs_out$suggested_sizes[[heuristic]]
    if (is.na(suggested_size)) {
      # Add an NA column for the heuristic
      prj_linpred_results[[paste0("prediction_", heuristic)]] <- NA
      next
    }
    
    # Extract selected predictors
    predictors_final <- head(ppvs_out$ranking[["fulldata"]], suggested_size)
    
    # Generate projections
    prj <- project(
      refm_obj,
      predictor_terms = predictors_final,
      verbose = FALSE
    )
    
    # Predict on test data
    prj_linpred <- proj_linpred(
      prj,
      newdata = data_test[, -1], # drop the dependent 'y' column
      integrated = TRUE
    )
    
    # Add predictions as a new column for the heuristic
    prj_linpred_results[[paste0("pred_", heuristic)]] <- as.vector(prj_linpred$pred)
  }
  
  # Save predictions for all heuristics
  save(prj_linpred_results, file = paste0("output/prj_linpred/prj_linpred_", name, ".RData"))
  
  # End timing
  end.time <- Sys.time()
  time.taken <- round(difftime(end.time, start.time, units = "secs"), 2)
  
  # Save timing information
  save(time.taken, file = paste0("output/time_taken/time_", name, ".RData"))
}

workflow_fun(2, conditions)
mean((prj_linpred_results$y - prj_linpred_results$pred_default)^2)

##### Part 3: Run in parallel -----


# Set the number of workers (number of parallel processes)
nworkers <- 4 # number of cores to use, max = 20
cl <- makeCluster(nworkers, type = "PSOCK")  # Create a PSOCK cluster

# Ensure all required libraries and the script are loaded on each worker
clusterEvalQ(cl, {
  library(tidyverse)
  library(brms)
  library(projpred)
  library(data.table)
  library(cmdstanr)
  source("workflow_functions.R")  # Load the script containing functions
})

# Export required objects to the workers
clusterExport(cl, varlist = c("conditions", "workflow_fun", "seed"),
              envir = environment())

# Run the workflow in parallel
#test_conditions <- conditions[1:4, ]

missing_conditions <- conditions[conditions$r %in% c(82, 87, 220, 412, 477, 483, 496, 291), ]
system.time(out <- clusterApplyLB(cl, 1:nrow(missing_conditions), workflow_fun, conditions = missing_conditions))

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

