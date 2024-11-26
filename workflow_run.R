
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
N <- c(100, 400)
corr <- c(0, 0.4, 0.8)
loc_TE <- c("first10", "mixed")
# r <- 1:100
r <- sample(1:500, size = 100, replace = FALSE)  # Randomly sample 100 values from 1:500


conditions <- expand.grid(N = N, corr = corr, loc_TE = loc_TE, r = r)
conditions$K <- ifelse(conditions$N == 100, 5, 2)
conditions <- conditions %>% 
  filter(!(corr == 0 & loc_TE == "mixed"))

#### Part 2: Set up workflow function -----
workflow_fun <- function(pos, conditions){
  
  source("workflow_functions.R")
  
  # extract condition and replication information
  sample_size <- conditions$N[pos]
  replication <- conditions$r[pos]
  corr <- conditions$corr[pos]
  loc_TE <- conditions$loc_TE[pos]
  K <- conditions$K[pos]
  
  # load and prepare data for the current condition and replication
  data <- load_and_prepare_data(file_path = NULL, # File path is constructed dynamically
                                sample_size = sample_size,
                                replication = replication,
                                corr = corr,
                                loc_TE = loc_TE)
  
  
  # create name for saving output
  name <- paste0("N", sample_size, "_corr", corr, "_loc_TE", loc_TE, "_r", replication)
  
  # fit the reference model
  refm_fit <- fit_reference_model(data)
  
  # save convergence info reference model and save the reference model summary
  summ_ref <- summary(refm_fit)
  summ_ref$fixed$div <- rstan::get_num_divergent(refm_fit$fit)
  save(summ_ref, file = paste0("output/summ_ref/summ_ref_", name, ".RData"))
  
  # run projection predictive variable selection
  ppvs_out <- run_projpred(refm_fit, K = K, nterms_max = 11)
  
  save(ppvs_out, file = paste0("output/ppvs_out/ppvs_out_", name, ".RData"))
  
  # add selection based on projpred to reference model summary
  selected_pred_ppvs <- head(ppvs_out$ranking[["fulldata"]], ppvs_out$suggested_size)
  df_out = summ_ref$fixed
  df_out$selected_pred_ppvs = 0
  df_out$selected_pred_ppvs[which(rownames(df_out) %in% selected_pred_ppvs)] = 1
  
  save(df_out, file = paste0("output/df_out/df_out_", name, ".RData"))
}

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
test_conditions <- conditions[1:4, ]

system.time(out <- clusterApplyLB(cl, 1:4, workflow_fun, conditions = test_conditions))

# Stop the cluster after execution
stopCluster(cl)





system.time(workflow_fun(3, conditions))

load("output/ppvs_out/ppvs_out_N100_corr0_loc_TEfirst10_r463.RData")














