# matrix rows are the draws, columns are the y_test

# check first the mean then the diff or other,
# then the quantiles
library(bayestestR)

# create a test marix
m_test = matrix(data = c(1,2,3,4,5,6,7,8), nrow = 4, ncol = 2)

# Example matrix
set.seed(123)
mat <- matrix(rnorm(1000), nrow = 100, ncol = 10) # 100 rows, 10 columns

# Function to compute mean and 95% confidence interval
summary_stats <- function(x) {
  mean_x <- mean(x)
  ci_x <- quantile(x, probs = c(0.025, 0.975))  # 95% interval
  return(c(mean = mean_x, lower = ci_x[1], upper = ci_x[2]))
}

# method number two
summary_stats2 <- function(x) {
  mean_x <- mean(x)
  ci_x <- bayestestR::ci(x, method = "HDI", ci = 0.95)
  return(c(mean = mean_x, lower = ci_x$CI_low, upper = ci_x$CI_high))
}

# Apply function to each column
result <- apply(ssvs_pred_results$pred, 2, summary_stats2)

# Convert to data frame for better readability
result_df <- as.data.frame(t(result))

# Print results
print(result_df)

# now we can check if a real value is inside the interval
y_real = rnorm(1000)

coverage = ((ssvs_pred_results$y > result_df$lower) & (ssvs_pred_results$y < result_df$upper))
mean(coverage)

