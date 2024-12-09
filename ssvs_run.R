
library(parallel)
library(SSVS)


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
# Workflow function for SSVS
workflow_ssvs <- function(pos, conditions) {

  source("workflow_functions.R")
  
  # Extract condition-specific information
  sample_size <- conditions$N[pos]
  replication <- conditions$r[pos]
  corr <- conditions$corr[pos]
  loc_TE <- conditions$loc_TE[pos]
  
  # Load and prepare the dataset
  data <- load_and_prepare_data(file_path = NULL, 
                                sample_size = sample_size,
                                replication = replication,
                                corr = corr,
                                loc_TE = loc_TE)
  
  # Define a unique name for saving the results
  name <- paste0("N", sample_size, "_corr", corr, "_loc_TE", loc_TE, "_r", replication)
  # Run SSVS
  set.seed(123)  # Ensure reproducibility
  ssvs_out <- ssvs(data = data, x = paste0("X", 1:50), y = "y", progress = TRUE)
  
  # Get summary and convert to dataframe
  ssvs_summary <- summary(ssvs_out, interval = 0.95)
  ssvs_df_out <- as.data.frame(ssvs_summary)
  
  # Add the column `selected_pred_ssvs`
  ssvs_df_out$selected_pred_ssvs <- ifelse(ssvs_df_out$MIP > 0.5, 1, 0)
  
  # Save the dataframe to an RData file
  save(ssvs_df_out, file = paste0("output/ssvs_df_out/ssvs_", name, ".RData"))
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
clusterExport(cl, varlist = c("conditions", "workflow_ssvs",  "load_and_prepare_data"), envir = environment())

# Run the workflow in parallel
system.time(out <- parallel::clusterApplyLB(cl, 1:nrow(conditions), workflow_ssvs, conditions = conditions))

# Stop the cluster
stopCluster(cl)

