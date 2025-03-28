#' DAP-S fine-mapping using sufficient statistics data
#' @param XtX A p by p matrix \eqn{X'X} in which the columns of X are centered to have mean zero.
#' @param Xty A p-vector \eqn{X'y} in which y and the columns of X are centered to have mean zero.
#' @param yty A scalar \eqn{y'y} in which y is centered to have mean zero.
#' @param n The sample size.
#' @param L Number of causal variants
#' @param prior_weights Prior weights
#' @param null_weight Null weight
#' @param residual_tau Residual variance
#' @param threshold Threshold for proposal density
#' @param phi2 Scaled prior effect size variance
#' @param phi2_vec Scaled prior effect size variance vector
#' @param r2_threshold Genotype R2 threshold for LD
#' @param coverage Coverage of credible sets. When not set, it outputs signal clusters;
#' otherwise, it outputs credible sets at the specified coverage level
#'
#' @import Rfast Matrix
#' @importFrom susieR susie_suff_stat
#' @importFrom dplyr %>% arrange
#' @return Fine mapping results
#'
#' @export
dap_suff_stat <- function(XtX, Xty, yty, n, L = min(10, ncol(XtX)),
                          susie_params = list(),
                          prior_weights = NULL,
                          standardize = TRUE,
                          use_susie_variance_estimate = TRUE,
                          grid = c(0.04, 0.16, 0.64),
                          overlapping = TRUE,
                          pir_threshold = 1e-6,
                          r2_threshold = 0.25,
                          coverage = NULL,
                          twas_weight = FALSE) {

  cat("Processing inputs...\n")
  if (standardize) { # standardize XtX and Xty
    dXtX <- diag(XtX)
    csd <- sqrt(dXtX / (n - 1))
    csd[csd == 0] <- 1
    XtX <- t((1 / csd) * XtX) / csd
    Xty <- Xty / csd
  }

  p <- ncol(XtX)
  if (is.null(colnames(XtX))) {
    snp_names <- paste0("SNP_", 1:p)
  } else {
    snp_names <- colnames(XtX)
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
  print(isSymmetric(XtX))
  susie_fit <- susie_suff_stat(XtX, Xty, yty, n, L,
                               scaled_prior_variance = 0.2,
                               prior_weights = prior_weights,
                               null_weight = susie_params$null_weight,
                               standardize = standardize,
                               estimate_residual_variance = susie_params$estimate_residual_variance,
                               estimate_prior_variance = susie_params$estimate_prior_variance)
  info <- param_setup(yty/(n-1), susie_fit, overlapping,
                      use_susie_variance_estimate, grid)

  cat("Running DAP-S fine-mapping...\n")
  results <- dap_main(matrix = info$mat,
                      pir_threshold = pir_threshold,
                      prior_weights = prior_weights,
                      phi2_mat = info$phi2_mat,
                      r2_threshold = r2_threshold,
                      coverage = coverage,
                      overlapping = overlapping,
                      twas_weight = twas_weight,
                      snp_names = snp_names,
                      ss = 1,
                      XtX_input = XtX,
                      Xty_input = Xty,
                      yty_input = yty,
                      n_input = n)

  cat("Summarizing results...\n")
  output <- get_summarization(susie_params, info, results,
                              pir_threshold, r2_threshold,
                              prior_weights, snp_names,
                              coverage, twas_weight)

  cat("Done!\n")
  return(output)
}
