#' DAP-S fine mapping with infinitesimal effects
#' @param bhat Estimated effect size
#' @param shat Standard error of bhat
#' @param z Vector of z-scores
#' @param var_y Variance of y
#' @param n Sample size
#' @param L Number of modeled causal effects
#' @param LD LD matrix
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
#' @examples
#' # Set random seed for reproducibility
#' set.seed(2025)
#' n <- 5000  # Number of samples
#' p <- 500   # Number of SNPs
#' MAF <- 0.1  # Minor allele frequency
#' Ltrue <- 5  # Number of true causal variants
#' ssq <- 0.01  # Effect size variance
#' sigmasq <- 1  # Error variance
#'
#' # Generate genotype matrix X
#' X <- matrix(rbinom(n*p, 2, MAF), nrow=n, ncol=p)
#' X <- scale(X, center=TRUE, scale=FALSE)
#'
#' # Generate sparse effects
#' b <- rep(0, p)
#' inds <- sample(1:p, size=Ltrue, replace=FALSE)
#' b[inds] <- rnorm(Ltrue) * sqrt(ssq)
#' order_idx <- order(inds)
#' cat('True effects:', inds[order_idx], '\n')
#' cat('Effect sizes:', b[inds[order_idx]], '\n')
#'
#' # Generate infinitesimal effects
#' tausq <- 1e-3
#' infinitesimal <- X %*% rnorm(p) * sqrt(tausq)
#' effects <- X %*% b + infinitesimal
#' y <- effects + rnorm(n) * sqrt(sigmasq)
#' y <- scale(y, center = TRUE, scale = FALSE)
#' cat('Total fraction of variance explained by SNPs:', var(as.vector(effects))/var(y), '\n')
#'
#' # Generate summary statistics
#' res = susieR::univariate_regression(X, y)
#' suff_stat <- susieR::compute_suff_stat(X, as.vector(y), standardize = FALSE)
#' LD <- cov2cor(suff_stat$XtX)
#' z <- res$betahat/res$sebetahat
#'
#' # Work on original scale
#' rst <- dap_inf(bhat = res$betahat, shat = res$sebetahat, var_y = var(y), n=n, L=5, LD = LD, twas_weight=TRUE)
#' rst$sigmasq
#' rst$tausq
#'
#' # SuSiE-inf
#' output <- susie_inf(bhat = res$betahat, shat = res$sebetahat, var_y = var(y), n=n, L=5, LD = LD, null_weight = (1-1/p)^p)
#' 
#' # PIP comparison
#' plot(output$spip, rst$pip, xlab = "PIP of SuSiE-inf", ylab = "PIP of DAP-S")
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#' 
#' # TWAS weight
#' plot(rowSums(output$mu * output$PIP)[1:p], rst$twas_weights,
#'     xlab = "Coefficients (sparse) of SuSiE-inf",
#'     ylab = "Coefficients (sparse) of DAP-S") + abline(a=0,b=1)
#'plot(rowSums(output$mu * output$PIP)[1:p] + output$alpha[1:p], rst$twas_weights,
#'     xlab = "Coefficients (sparse+infinitesimal) of SuSiE-inf",
#'     ylab = "Coefficients (sparse) of DAP-S") + abline(a=0,b=1)
dap_inf <- function(bhat = NULL, shat = NULL, z = NULL, var_y = NULL,
                    n, LD = NULL, L = 10,
                    prior_weights = NULL, standardize = TRUE,
                    use_susie_variance_estimate = TRUE,
                    grid = c(0.04, 0.16, 0.64),
                    overlapping = TRUE,
                    pir_threshold = 1e-6,
                    r2_threshold = 0.25,
                    coverage = NULL,
                    method = "moments",
                    twas_weight = FALSE) {

  if (!is.null(bhat)) {
    p <- length(bhat)
  } else if (!is.null(z)) {
    p <- length(z)
  }

  if (!is.null(colnames(LD))) {
    snp_names <- colnames(LD)
  } else {
    snp_names <- paste0("SNP_", 1:p)
  }

  if (is.null(prior_weights)) prior_weights <- rep(1 / p, p)
  if (is.null(coverage)) coverage <- 2

  cat("Running SuSiE-inf...\n")
  susie_fit <- susie_inf(bhat = bhat, shat = shat, z = z, var_y = var_y,
                      n = n, LD = LD, L = L, method = method,
                      pi = prior_weights, null_weight = prod(1 - prior_weights))
  fit <- NULL
  fit$alpha <- t(susie_fit$PIP)
  fit$V <- susie_fit$ssq
  info <- param_setup(var_y, fit, overlapping,
                      use_susie_variance_estimate, grid)

  susie_params <- list(
    null_weight = prod(1 - prior_weights)
  )

  cat("Running DAP-S fine-mapping with infinitesimal effects...\n")
  results <- dap_main(matrix = info$mat,
                      pir_threshold = pir_threshold,
                      prior_weights = prior_weights,
                      phi2_mat = info$phi2_mat,
                      r2_threshold = r2_threshold,
                      coverage = coverage,
                      overlapping = overlapping,
                      twas_weight = twas_weight,
                      snp_names = snp_names,
                      ss = 2, # meaning infinitesimal effect
                      V_input = susie_fit$V,
                      Dsq_input = susie_fit$Dsq,
                      var_input = susie_fit$var,
                      XtOmegay_input = susie_fit$XtOmegay)

  cat("Summarizing results...\n")
  output <- get_summarization(susie_params, info, results,
                              pir_threshold, r2_threshold,
                              prior_weights, snp_names,
                              coverage, twas_weight)
  output$sigmasq <- susie_fit$sigmasq
  output$tausq <- susie_fit$tausq
  cat("Done!\n")
  return(output)
}
