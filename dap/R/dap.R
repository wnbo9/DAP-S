#' DAP-S fine-mapping with individual-level data.
#' @param X Genotype matrix (n x p)
#' @param y Phenotype vector (n)
#' @param L Number of causal variants. Default is min(10, p).
#' @param susie_params Parameters for SuSiE: \code{scaled_prior_variance},
#' \code{residual_variance}, \code{null_weight}, \code{standardize},
#' \code{estimate_residual_variance}, \code{estimate_prior_variance}.
#' @param scaled_prior_variance Scaled prior variance. That is, the prior
#' variance of each non-zero element of b is \code{var(y)} *
#' \code{scaled_prior_variance}. Default is 0.2.
#' @param residual_variance Residual variance. Default is NULL.
#' @param null_weight Null weight used in SuSiE. Default is calculated as
#' \code{prod(1 - prior_weights)}.
#' @param prior_weights Vector of prior probabilities for each variant being
#'   causal. Default is uniform priors
#' @param standardize Standardize each column of the genotype matrix. Default is
#' TRUE.
#' @param use_susie_variance_estimate Use SuSiE-estimated prior variance to
#' guide the grid upper bound. Default is TRUE.
#' @param grid Grid of prior variances, default is c(0.04, 0.16, 0.64). The
#' upper bound is shifted to 1.28 if the maximum of SuSiE-estimated prior
#' variance is greater than 1.
#' @param overlapping Overlapping signal clusters. Default is TRUE.
#' @param pir_threshold Threshold for PIR. Default is 1e-6.
#' @param r2_threshold Threshold for genotype R2 when constructing signal
#' clusters. Default is 0.25.
#' @param coverage Coverage for credible set. If NULL, the algorithm will
#' construct signal clusters.
#' @param twas_weight Compute TWAS weights. Default is FALSE.
#' @import Rfast
#' @importFrom susieR susie
#' @importFrom dplyr %>% arrange desc
#' @return List of fine-mapping results.
#' @export
#'
#' @examples
#' # DAP-S fine-mapping example
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0, p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' res = dap(X, y, L = 10)
#'
#' # Marginal PIP
#' res$pip
#' # Signal clusters
#' res$sets
#' # Variants and PIP in signal clusters
#' res$variants
#'
#' # 95% credible set
#' res = dap(X, y, L = 10, coverage = 0.95)
#' res$sets
#'
#' # TWAS weights
#' res = dap(X, y, L = 10, twas_weight = TRUE)
#' res$twas_weights
#'
dap <- function(X, y, L = min(10, ncol(X)),
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
    null_weight = prod(1 - prior_weights)
  )
  susie_params <- modifyList(default_susie_params, susie_params)

  cat("Running SuSiE...\n")
  susie_fit <- susie(X, y, L,
                     scaled_prior_variance = susie_params$scaled_prior_variance,
                     residual_variance = susie_params$residual_variance,
                     prior_weights = prior_weights,
                     null_weight = susie_params$null_weight,
                     standardize = standardize)
  info <- param_setup(y, susie_fit, overlapping,
                      use_susie_variance_estimate, grid)

  cat("Running DAP-S fine-mapping...\n")
  results <- dap_main(X, y, info$mat, pir_threshold,
                      prior_weights, info$phi2_mat,
                      r2_threshold, coverage, overlapping,
                      twas_weight, snp_names)

  cat("Summarizing results...\n")
  output <- get_summarization(susie_params, info, results,
                              pir_threshold, r2_threshold,
                              prior_weights, snp_names,
                              coverage, twas_weight)

  cat("Done!\n")
  return(output)
}