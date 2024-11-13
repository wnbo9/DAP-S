#' DAP-PIR
#' @param X Genotype data
#' @param y Phenotype data
#' @param L Number of causal variants
#' @param prior_weights Prior weights
#' @param null_weight Null weight
#' @param residual_tau Residual variance
#' @param threshold Threshold for proposal density
#' @param phi2 Scaled prior effect size variance
#'
#' @import Rfast
#' @importFrom susieR susie
#' @return Fine mapping results
#'
#' @export
dap <- function(X, y, L,
                prior_weights = NULL,
                null_weight = NULL,
                residual_tau = NULL,
                threshold = NULL,
                phi2 = NULL) {

  # load data
  print('Loading input...')
  ## check input
  X <- scale(X, scale = FALSE)
  y <- scale(y, scale = FALSE)
  n <- length(y)
  p <- ncol(X)
  ## initialize
  if (is.null(null_weight)) {
    null_weight = (1-1/p)^p
  }
  if (is.null(prior_weights)) {
    prior_weights = rep(1/p, p)
  }
  if (is.null(phi2)) {
    phi2 = 0.6^2
  }
  if (is.null(residual_tau)) {
    residual_tau = 1/var(y)
  }
  if (is.null(threshold)) {
    threshold = 1e-6
  }





  # pseudo importance resampling
  print("Pseudo importance resampling...")
  rst <- susie(X, y, L, max_iter = 1000, coverage = 0.95, standardize = FALSE,
               null_weight = null_weight, prior_weights = prior_weights,
               estimate_residual_variance = FALSE, residual_variance = 1/residual_tau,
               estimate_prior_variance = FALSE, scaled_prior_variance = phi2)
  matrix <- rst$alpha

  print(dim(matrix))

  # deterministic approximation of posteriors
  print("Running DAP...")

  # results
}
