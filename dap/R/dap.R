#' DAP-PIR
#' @param X Genotype data
#' @param y Phenotype data
#' @param L Number of causal variants
#' @param prior_weights Prior weights
#' @param null_weight Null weight
#' @param residual_tau Residual variance
#' @param threshold Threshold for proposal density
#' @param phi2 Scaled prior effect size variance
#' @param phi2_vec Scaled prior effect size variance vector
#'
#' @import Rfast
#' @importFrom susieR susie
#' @importFrom dplyr %>% arrange
#' @return Fine mapping results
#'
#' @export
dap <- function(X, y, L,
                prior_weights = NULL,
                null_weight = NULL,
                residual_tau = NULL,
                threshold = 1e-6,
                phi2 = 0.36,
                phi2_vec = NULL) {

  # load data
  print("Loading data...")
  ## check input
  X <- scale(X, scale = FALSE)
  y <- scale(y, scale = FALSE)
  n <- length(y)
  p <- ncol(X)
  ## initialize
  if (is.null(null_weight)) {
    null_weight <- (1-1/p)^p
  }
  if (is.null(prior_weights)) {
    prior_weights <- rep(1/p, p)
  }
  if (is.null(phi2_vec)) {
    phi2_vec <- c(0.04, 0.16, 0.36, 0.64)
  }
  if (is.null(residual_tau)) {
    residual_tau <- 1/var(y)
  }

  # pseudo importance resampling
  print("Pseudo importance resampling...")
  rst <- susie(X, y, L, max_iter = 1000, coverage = 0.95, standardize = FALSE,
               null_weight = null_weight, prior_weights = prior_weights,
               estimate_residual_variance = FALSE, residual_variance = 1/residual_tau,
               estimate_prior_variance = FALSE, scaled_prior_variance = phi2)
  matrix <- t(rst$alpha)
  mcfg_mat <- pir(matrix, threshold = threshold)
  mcfg_mat <- mcfg_mat[, 1:p]
  mcfg_mat <- unique(mcfg_mat)
  m_size <- nrow(mcfg_mat)

  # deterministic approximation of posteriors
  print(paste0("Calculating posterior of ", m_size, " models..."))
  log10_posterior <- compute_log10_posterior(X, y, mcfg_mat, prior_weights, phi2_vec)

  # print results
  print("Summarizing fine mapping results...")
  max_log_posterior <- max(log10_posterior)
  log_nc <- max_log_posterior + log10(sum(10^(log10_posterior - max_log_posterior)))
  posterior_probs <- 10^(log10_posterior - log_nc)
  # Compute the PIP for each variant using matrix multiplication
  pip_vector <- posterior_probs %*% mcfg_mat
  pip_vector <- round(pip_vector, digits = 4)
  pip_vector <- as.vector(pip_vector)
  model_configurations <- apply(mcfg_mat, 1, function(row) {
    if (all(row == 0)) {
      return("NULL")
    } else {
      return(paste(which(row == 1), collapse = "+"))
    }
  })
  result_df <- data.frame(Model = model_configurations, Log10_Posterior = log10_posterior, stringsAsFactors = FALSE)
  result_df <- result_df %>% arrange(desc(Log10_Posterior))

  print("Done!")
  # Return a list containing the PIP vector and the dataframe
  return(list(pip = pip_vector, models = result_df, log10_nc = log_nc))
}
