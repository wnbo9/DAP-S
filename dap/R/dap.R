#' DAP-S fine-mapping with individual-level data.
#' @param X Genotype matrix (n x p)
#' @param y Phenotype vector (n)
#' @param L Number of causal variants
#' @param susie_params Parameters for SuSiE: \code{scaled_prior_variance},
#' \code{residual_variance}, \code{null_weight}, \code{standardize},
#' \code{estimate_residual_variance}, \code{estimate_prior_variance}.
#' @param scaled_prior_variance Scaled prior variance. That is, the prior
#' variance of each non-zero element of b is \code{var(y)} *
#' \code{scaled_prior_variance}.
#' @param residual_variance Residual variance
#' @param prior_weights Prior weights of each SNP being causal
#' @param null_weight Null weight used in SuSiE
#' @param standardize Standardize each column of the genotype matrix
#' @param estimate_residual_variance Estimate residual variance in SuSiE
#' @param estimate_prior_variance Estimate prior variance in SuSiE
#' @param use_susie_variance_estimate Use SuSiE-estimated prior variance
#' @param grid Grid of prior variances, default is c(0.04, 0.16, 0.64)
#' @param exclusive Exclusive signal clusters
#' @param pir_threshold Threshold for PIR
#' @param r2_threshold Threshold for R2 when constructing credible sets
#' @param coverage Coverage for credible set
#' @import Rfast
#' @importFrom susieR susie
#' @importFrom dplyr %>% arrange desc
#' @return List of fine-mapping results.
#' @export
dap <- function(X, y, L = min(10, ncol(X)),
                susie_params = list(),
                prior_weights = NULL,
                standardize = TRUE,
                use_susie_variance_estimate = FALSE,
                grid = c(0.04, 0.16, 0.64),
                exclusive = TRUE,
                pir_threshold = 1e-6,
                r2_threshold = 0.25,
                coverage = NULL) {

  cat("Processing inputs...\n")
  X <- scale(X, scale = standardize)
  y <- scale(y, scale = FALSE)
  p <- ncol(X)
  if (is.null(colnames(X))) {
    snp_names <- paste0("SNP_", 1:p)
  } else {
    snp_names <- colnames(X)
  }

  cat("Initializing parameters...\n")
  if (is.null(prior_weights)) prior_weights <- rep(1 / p, p)
  if (is.null(coverage)) coverage <- 2

  default_susie_params <- list(
    scaled_prior_variance = 0.2,
    residual_variance = NULL,
    null_weight = prod(1 - prior_weights),
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE
  )
  susie_params <- modifyList(default_susie_params, susie_params)

  cat("Running SuSiE...\n")
  susie_fit <- susie(X, y, L, max_iter = 1000, coverage = 0.95,
                     scaled_prior_variance = susie_params$scaled_prior_variance,
                     residual_variance = susie_params$residual_variance,
                     prior_weights = prior_weights,
                     null_weight = susie_params$null_weight,
                     standardize = standardize,
                     estimate_residual_variance = susie_params$estimate_residual_variance,
                     estimate_prior_variance = susie_params$estimate_prior_variance)
  info <- param_setup(y, susie_fit, exclusive, use_susie_variance_estimate, grid)

  cat("Running DAP-S fine-mapping...\n")
  results <- dap_main(X, y, info$mat, pir_threshold,
                      prior_weights, info$phi2_mat,
                      r2_threshold, coverage, exclusive, snp_names)

  cat("Summarizing results...\n")
  output <- get_summarization(susie_params, info, results,
                              pir_threshold, r2_threshold,
                              prior_weights, snp_names)

  cat("Done!\n")
  return(output)
}