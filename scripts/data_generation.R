
# load libraries
library(Matrix)
library(MASS)

#################################### Study 1 test data generation ####################################
##### data generation function #####--------------------------------------------------------------------------------
data_gen_study1 <- function(N, bc, pattern_TE) {
  p <- 50
  t <- 10
  mu <- rep(0, p)
  
  # Create block-diagonal correlation matrix
  Sigma0 <- bdiag(replicate(10, {
    block <- matrix(bc, nrow = 5, ncol = 5)
    diag(block) <- 1
    block
  }, simplify = FALSE))
  Sigma0 <- as.matrix(Sigma0)
  
  # Set beta coefficients
  if (pattern_TE == "mixed") {
    beta <- c(1.5, rep(0, 4), 0.9, rep(0, 4), 0.9, rep(0, 4), 0.3, rep(0, 4), 0.3, rep(0, 4), 
              1.5, rep(0, 4), 0.9, rep(0, 4), 0.9, rep(0, 4), 0.3, rep(0, 4), 0.3, rep(0, 4))
  } else if (pattern_TE == "clustered") {
    beta <- c(rep(c(1.5, 0.9, 0.9, 0.3, 0.3), 2), rep(0, p - t))
  } else {
    stop("Invalid pattern_TE value. Use 'mixed' or 'clustered'.")
  }
  
  # Generate one dataset
  X <- mvrnorm(n = N, mu = mu, Sigma = Sigma0)
  colnames(X) <- paste0("X", 1:p)
  y <- X %*% beta + rnorm(N)
  df <- data.frame(y, X)
  return(df)
}

##### data generation #####--------------------------------------------------------------------------------
set.seed(123)
bc_vals <- c(0, 0.4, 0.8)
patterns_all <- c("clustered", "mixed")
n_train <- 100
n_test <- 1000
r <- 100  # number of replications

for (bc in bc_vals) {
  patterns <- if (bc == 0) "clustered" else patterns_all
  
  for (pattern in patterns) {
    df_train <- vector("list", r)
    df_test <- vector("list", r)
    
    for (i in 1:r) {
      df_train[[i]] <- data_gen_study1(N = n_train, bc = bc, pattern_TE = pattern)
      df_test[[i]]  <- data_gen_study1(N = n_test, bc = bc, pattern_TE = pattern)
    }
    
    # Save in RData format
    
    train_filename <- paste0("data/p50/train_n", n_train, "_p50_corr", bc, "_TE", pattern, ".RData")
    test_filename  <- paste0("data/p50/test_fortrain_n", n_train, "_p50_corr", bc, "_TE", pattern, ".RData")
    
    save(df_train, file = train_filename)
    save(df_test, file = test_filename)
    
    cat("Saved:", train_filename, "and", test_filename, "\n")
  }
}


#################################### Study 2 data generation ####################################
##### data generation function #####--------------------------------------------------------------------------------
data_gen <- function(rho, n, pattern_TE) {
  
  # number of parameters in the model
  p <- 75
  
  # define xi given rho and pattern of true effect
  if (pattern_TE == "clustered"){
    if (rho == 0){
      xi <- 0.59
    } else if (rho == 0.2){
      xi <- 0.44
    } else if (rho == 0.8){
      xi <- 0.29
    } else if (rho == 0.5){
      xi <- 0.34 #
    } else {
      warning("Invalid rho provided.")
    }
  } else if (pattern_TE == "mixed") {
    xi <- 0.59
  } else {
    warning("Invalid pattern of true effects provided")
  }
  
  
  
  # build correlation matrix
  block_size <- 5
  num_matrices <- p / block_size
  listOfMatrices <- vector("list", num_matrices)
  for (i in 1:num_matrices) {
    listOfMatrices[[i]] <- matrix(
      rep(rho, block_size * block_size), nrow=block_size
    )
  }
  R <- matrix(bdiag(listOfMatrices), nrow=p)
  diag(R) <- 1
  
  # build associated covariate coefficients
  if (pattern_TE == "clustered"){
    beta1 <- rep(xi, block_size)
    beta2 <- rep(xi * 0.5, block_size)
    beta3 <- rep(xi * 0.25, block_size)
    beta4 <- rep(0, (num_matrices - 3) * block_size)
    beta <- c(beta1, beta2, beta3, beta4)
    
  } else if (pattern_TE == "mixed"){
    beta <- rep(0, p)
    true_indices <- seq(from = 1, by = 5, length.out = 15)
    beta[true_indices[1:5]] <- xi
    beta[true_indices[6:10]] <- xi * 0.5
    beta[true_indices[11:15]] <- xi * 0.25
    
  } else {
    warning("Invalid pattern of true effects provided.")
  }
  
  # # define zero mean vector
  sigma <- 1
  mu <- rep(0, p)
  # sds <- rep(sigma, p)
  # S <- .cor2cov(R, sds)
  
  # sample data
  # X <- mvrnorm(n=n, mu=mu, Sigma=S)
  X <- mvrnorm(n=n, mu=mu, Sigma=R) # sigma = 1 so here the cor and cov matrices are the same
  X <- as.data.frame(X)
  colnames(X) <- paste0("X", 1:ncol(X))
  y <- as.matrix(X) %*% beta + sigma * rnorm(n=n)
  data <- data.frame(y = y, X)
  
  return(data)
}


##### data generation #####--------------------------------------------------------------------------------

# set seed for reproducibility
set.seed(123)

# Define all conditions
rho_values <- c(0.2, 0.8)          # Correlation levels
n_train_values <- c(50, 75)        # Training sample sizes
pattern_TE_values <- c("clustered", "mixed")  # True effect patterns
n_test <- 1000                     # Fixed test set size
r <- 100                           # Number of replications per condition

# Loop through all conditions
for (rho in rho_values) {
  for (n_train in n_train_values) {
    for (pattern_TE in pattern_TE_values) {
      
      # Create lists to store 100 replications
      df_train <- vector("list", r)
      df_test <- vector("list", r)
      
      # Generate 100 replications
      for (i in 1:r) {
        df_train[[i]] <- data_gen(rho, n_train, pattern_TE)  # Generate training data
        df_test[[i]] <- data_gen(rho, n_test, pattern_TE)  # Generate test data
      }
      
      # Create unique filenames
      train_filename <- paste0("data/simdata/train_n", n_train, "_p75_corr", rho, "_TE", pattern_TE, ".RData")
      test_filename <- paste0("data/simdata/test_fortrain_n", n_train, "_p75_corr", rho, "_TE", pattern_TE, ".RData")
      
      # Save lists of 100 replications
      save(df_train, file = train_filename)
      save(df_test, file = test_filename)
      
      # Print progress update
      cat("Saved:", train_filename, "and", test_filename, "\n")
      
    }
  }
}
