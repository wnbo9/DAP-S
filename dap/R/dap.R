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
#' set.seed(1234)
#' n = 1000
#' p = 1000
#' beta = rep(0, p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' y = X %*% beta + rnorm(n)
#' X = scale(X,center = TRUE,scale = FALSE)
#' y = scale(y,center = TRUE,scale = FALSE)
#' y = as.vector(y)
#' 
#' # run DAP-S
#' res <- dap(X, y, L = 10, standardize = FALSE)
#' # Marginal PIP
#' res$pip
#' # Signal clusters
#' res$sets
#' # Variants and PIP in signal clusters
#' res$variants
#' 
#' # 95% credible set
#' res <- dap(X, y, L = 10, coverage = 0.95, standardize = FALSE)
#' res$sets
#' 
#' # DAP-S with sufficient statistics (XtX, Xty, yty, n), results identical to original individual-level data
#' suff_stat <- susieR::compute_suff_stat(X, y, standardize = FALSE)
#' rss <- with(suff_stat, dap_suff_stat(XtX, Xty, yty, n, L=10, standardize = FALSE))
#' plot(res$pip, rss$pip)
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#' 
#' # DAP-S with summary statistics
#' # (1) summary statistics (bhat, shat, R, var_y, n), results identical to original individual-level data
#' stat <- susieR::univariate_regression(X, y)
#' R <- cov2cor(suff_stat$XtX)
#' rss1 <- dap_rss(bhat = stat$betahat, shat = stat$sebetahat, R = R, n = n, var_y = var(y), estimate_residual_variance = TRUE)
#' plot(res$pip, rss1$pip)
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#' # (2) summary statistics (zhat, R, n), results identical to standardized individual-level data; similar to original individual-level data
#' zhat <- stat$betahat / stat$sebetahat
#' rss2 <- dap_rss(zhat, R, n = n, estimate_residual_variance = TRUE)
#' plot(res$pip, rss2$pip)
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#' res2 <- dap(X, y, L = 10, standardize = TRUE)
#' plot(res2$pip, rss2$pip)
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#' # (3) summary statistics (bhat, shat, R, n), results identical to standardized individual-level data; similar to original individual-level data
#' rss3 <- dap_rss(bhat = stat$betahat, shat = stat$sebetahat, R = R, n = n, estimate_residual_variance = TRUE)
#' plot(res2$pip, rss3$pip)
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' # TWAS weights
#' res <- dap(X, y, L = 10, twas_weight = TRUE)
#' res$twas_weights
#'
#' # Update DAP-S prior
#' vec = rep(1/p, p)
#' vec[5] = 10/p
#' vec = vec/sum(vec)
#' res2 <- dap_update(X=X, dap_result = res, prior_weights = vec)
#' res2$pip[1:5]

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
    null_weight = prod(1 - prior_weights),
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE
  )
  susie_params <- modifyList(default_susie_params, susie_params)

  cat("Running SuSiE...\n")
  susie_fit <- susie(X, y, L,
                     scaled_prior_variance = susie_params$scaled_prior_variance,
                     residual_variance = susie_params$residual_variance,
                     prior_weights = prior_weights,
                     null_weight = susie_params$null_weight,
                     standardize = standardize,
                     estimate_residual_variance = susie_params$estimate_residual_variance,
                     estimate_prior_variance = susie_params$estimate_prior_variance)
  info <- param_setup(var(y), susie_fit, overlapping,
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
                      ss = 0,
                      X_input = X,
                      y_input = y)

  cat("Summarizing results...\n")
  output <- get_summarization(susie_params, info, results,
                              pir_threshold, r2_threshold,
                              prior_weights, snp_names,
                              coverage, twas_weight)

  cat("Done!\n")
  return(output)
}