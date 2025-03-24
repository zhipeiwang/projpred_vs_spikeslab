
# load, filter and standardize data
load_and_prepare_data <- function(sample_size, replication, corr, TE) {
  
  # Construct the file name based on the given conditions
  data_dir <- "simdata/"
  train_file <- paste0(data_dir, "train_n", sample_size, "_p75_corr", corr, "_TE", TE, ".RData")
  test_file  <- paste0(data_dir, "test_fortrain_n", sample_size, "_p75_corr", corr, "_TE", TE, ".RData")
  
  # Load the correct training dataset
  if (file.exists(train_file)) {
    load(train_file)  # This loads `df_train`
  } else {
    stop(paste("Training data file not found:", train_file))
  }
  
  # Load the correct test dataset
  if (file.exists(test_file)) {
    load(test_file)  # This loads `df_test`
  } else {
    stop(paste("Test data file not found:", test_file))
  }
  
  # Extract the correct replication
  if (replication > length(df_train) || replication > length(df_test)) {
    stop("Replication index exceeds available replications.")
  }
  
  # Select the appropriate replication
  train_data <- df_train[[replication]]
  test_data <- df_test[[replication]]
  
  return(list(train = train_data, test = test_data))
}



# fit the reference model and save for projpred
fit_reference_model <- function(data, seed = 123) {
  refm_fit <- brm(
    formula = y ~ .,
    data = data,
    family = gaussian(),
    prior = prior(horseshoe(), class = "b"),
    chains = 4,
    iter = 2000,
#    backend = "cmdstanr",
    seed = seed
  )
  return(refm_fit)
}

# run cross-validation for variable selection
run_projpred <- function(refm_fit, K, nterms_max = 16, seed = 123) {
  
  # get the reference model object
  refm_obj <- get_refmodel(refm_fit)
  
  cvvs <- cv_varsel(
    refm_obj,
    cv_method = "kfold",
    K = K,
    method = "forward",
    nclusters = 20,
    ndraws_pred = 400,
    nterms_max = nterms_max,
    verbose = FALSE,
    seed = seed
  )
  
  
  # Try to suggest size
  # Compute various suggest_size heuristics
  suggested_sizes <- list(
    default = tryCatch({
      suggest_size(cvvs, type = "upper")
    }, error = function(e) {
      cat("Error in suggest_size (default):", e$message, "\n")
      return(NA)
    }),
    lower = tryCatch({
      suggest_size(cvvs, type = "lower")
    }, error = function(e) {
      cat("Error in suggest_size (lower):", e$message, "\n")
      return(NA)
    }),
    thres_upper = tryCatch({
      suggest_size(cvvs, thres_elpd = -4, type = "upper")
    }, error = function(e) {
      cat("Error in suggest_size (thres_upper):", e$message, "\n")
      return(NA)
    }),
    thres_lower = tryCatch({
      suggest_size(cvvs, thres_elpd = -4, type = "lower")
    }, error = function(e) {
      cat("Error in suggest_size (thres_lower):", e$message, "\n")
      return(NA)
    })
  )
  
  
  ranking <- ranking(cvvs)
  summary_cvvs <- summary(cvvs)
  
  output_projpred <- list(
    suggested_sizes = suggested_sizes,
    ranking = ranking,
    summary_cvvs = summary_cvvs,
    cvvs = cvvs
  )
  
  return(output_projpred)
}

