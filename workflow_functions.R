
# load, filter and standardize data
# load_and_prepare_data <- function(file_path, sample_size, replication) {
#   data <- fread(file_path)[N == sample_size & r == replication]
#   data <- as.data.frame(data) %>% select(-V1, -N, -r)  # Adjust 'V1' as needed based on the data loading method
#   data[ , -1] <- scale(data[ , -1])  # Standardize all columns except dependent variable 'y'
#   return(data)
# }
load_and_prepare_data <- function(file_path, sample_size, replication, corr, loc_TE) {
  # Define the file name mapping
  file_info <- data.frame(
    corr = c(0, 0.4, 0.8, 0.4, 0.8),
    loc_TE = c("first10", "first10", "first10", "mixed", "mixed"),
    file_name = c(
      "sim_p50_t10",
      "sim_p50_t10_corr0.4_first10",
      "sim_p50_t10_corr0.8_first10",
      "sim_p50_t10_corr0.4_mixed",
      "sim_p50_t10_corr0.8_mixed"
    )
  )
  
  # Match the correct file based on corr and loc_TE
  file_row <- file_info[file_info$corr == corr & file_info$loc_TE == loc_TE, ]
  if (nrow(file_row) == 0) {
    stop("No matching file found for the specified corr and loc_TE.")
  }
  
  # Construct the full file path
  selected_file <- paste0("p50/", file_row$file_name, ".csv")
  if (!file.exists(selected_file)) {
    stop(paste("File not found:", selected_file))
  }
  
  # Load and filter the data
  data <- fread(selected_file)[N == sample_size & r == replication]
  
  # Convert to data frame and remove unnecessary columns
  data <- as.data.frame(data) %>% select(-V1, -N, -r)
  
  # Standardize all columns except the dependent variable 'y'
  data[ , -1] <- scale(data[ , -1])
  
  return(data)
}


# fit the reference model and save for projpred
fit_reference_model <- function(data, seed = 123) {
  refm_fit <- brm(
    formula = y ~ .,
    data = data,
    family = gaussian(),
    prior = prior(horseshoe(par_ratio = 0.2), class = "b"),
    chains = 4,
    iter = 2000,
    backend = "cmdstanr",
    seed = seed
  )
  return(refm_fit)
}

# run cross-validation for variable selection
run_projpred <- function(refm_fit, K, nterms_max = 11, seed = 123) {
  
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

