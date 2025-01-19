#' DAP-S fine-mapping with individual-level data. The parameters are directly obtained from SuSiE results. # nolint: line_length_linter.
#' @param X Genotype data
#' @param y Phenotype data
#' @param L Number of causal variants
#' @param scaled_prior_variance Scaled prior effect size variance. The prior variance, divided by \code{var(y)}; that is, the prior variance of each non-zero element of b is \code{var(y) * scaled_prior_variance}. So it is equal to phi^2. # nolint: line_length_linter.
#' @param residual_variance Residual variance. If \code{estimate_residual_variance = TRUE}, this value provides the initial estimate of the residual variance. # nolint: line_length_linter.
#' @param prior_weights Prior weights
#' @param null_weight Null weight. Probability of no effects, it should be between 0 and 1. # nolint: line_length_linter.
#' @param estimate_residual_variance If \code{estimate_residual_variance = TRUE}, the residual variance is estimated, using \code{residual_variance} as an initial value. If \code{estimate_residual_variance = FALSE}, the residual variance is fixed to the value supplied by \code{residual_variance}. # nolint: line_length_linter.
#' @param estimate_prior_variance If \code{estimate_prior_variance = TRUE}, the prior variance is estimated (this is a separate parameter for each of the L effects). If provided, \code{scaled_prior_variance} is then used as an initial value for the optimization. When \code{estimate_prior_variance = FALSE}, the prior variance for each of the L effects is determined by the value supplied to \code{scaled_prior_variance}. # nolint: line_length_linter.
#' @param r2_threshold Genotype R2 threshold for LD
#' @param coverage Coverage of credible sets. When not set, it outputs signal clusters; otherwise, it outputs credible sets at the specified coverage level # nolint: line_length_linter.
#' @param use_susie_variance_estimate If \code{use_susie_variance_estimate = TRUE}, the SuSiE variance estimate (V) is used to compute the Bayes factor. If \code{use_susie_variance_estimate = FALSE}, the Bayes factor is calculated by averaging over a grid of values provided in \code{phi2_vec}. # nolint: line_length_linter.
#' @param phi2_vec Scaled prior effect size variance vector. If \code{use_susie_variance_estimate = FALSE}, this vector is used to compute the Bayes factor. # nolint: line_length_linter.
#' @param greedy If \code{greedy = TRUE}, the greedy algorithm is used to find the best model, instead of conducting Pseudo-Importance Resampling. # nolint: line_length_linter.
#' @param pir_threshold Threshold for PIR when \code{greedy = FALSE}. # nolint: line_length_linter.
#' @import Rfast
#' @importFrom susieR susie
#' @importFrom dplyr %>% arrange
#' @return Fine mapping results
#'
#' @export
dap_susie <- function(X, y, L = min(10, ncol(X)),
                scaled_prior_variance = 0.2,
                residual_variance = NULL,
                prior_weights = NULL,
                null_weight = NULL,
                estimate_residual_variance = TRUE,
                estimate_prior_variance = TRUE,
                r2_threshold = 0.25,
                coverage = NULL,
                use_susie_variance_estimate = FALSE,
                phi2_vec = c(0.04, 0.16, 0.64),
                greedy = FALSE,
                pir_threshold = 1e-6) {

  # Process inputs
  cat("Processing inputs...\n")
  X <- scale(X, scale = FALSE)
  y <- scale(y, scale = FALSE)
  n <- length(y)
  p <- ncol(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("SNP_", 1:p)

  ## Initialize parameters
  if (is.null(null_weight)) null_weight <- (1-1/p)^p
  if (is.null(prior_weights)) prior_weights <- rep(1/p, p)
  if (is.null(coverage)) coverage <- 2

  # Run SuSiE
  susie_fit <- susie(X, y, L, max_iter = 1000, coverage = 0.95,
               scaled_prior_variance = scaled_prior_variance,
               residual_variance = residual_variance,
               null_weight = null_weight, prior_weights = prior_weights,
               estimate_residual_variance = estimate_residual_variance,
               estimate_prior_variance = estimate_prior_variance)
  matrix <- t(susie_fit$alpha)
  for (col in 1:sum(susie_fit$V!=0)) {
    if (matrix[p+1, col] == max(matrix[, col])) break
  }
  mat <- as.matrix(matrix[,1:col]) # We only keep at most one column with null as top among all non-zero effects.
  mat <- process_matrix(mat)  # we need to pre-process the mat, so that each SNP only appears in the effect where it has highest rank.

  # then we do PIR (or greedy) to get models
  results <- dap_susie_main(X, y, mat, prior_weights, r2_threshold, coverage, use_susie_variance_estimate, phi2_vec, greedy, pir_threshold)


  # then we do Posterior probability (either using SuSiE estimates, or using a grid of values)
  # then we do credible set or signal cluster

  # Create model results
  result_df <- data.frame(
    Model = results$model_config,
    Posterior_Prob = results$posterior_prob,
    Log10_Posterior_Score = results$log10_posterior_score,
    Log10_BF = results$log10_BF,
    Log10_Prior = results$log10_prior,
    stringsAsFactors = FALSE) %>% arrange(desc(Log10_Posterior_Score))
  params <- list(
    log10_nc = results$log10_nc,
    pir_threshold = pir_threshold,
    r2_threshold = r2_threshold,
    phi2 = scaled_prior_variance,
    phi2_vec = phi2_vec,
    prior_weights = prior_weights,
    null_weight = null_weight)
  snp <- data.frame(
    SNP = colnames(X),
    PIP = results$pip,
    cluster_index = -1,
    in_cluster_order = Inf,
    stringsAsFactors = FALSE)

  # Make sure each SNP is included in the cluster where it has higher order
  if (!is.null(results$signal_cluster)) {
    for (i in seq_along(results$signal_cluster$clusters)) {
      for (order in seq_along(results$signal_cluster$clusters[[i]])) {
        snp_name <- results$signal_cluster$clusters[[i]][order]
        idx <- which(snp$SNP == snp_name)

        if (snp$cluster_index[idx] != -1) {
          # Only reassign if new order is better than current order
          if (order < snp$in_cluster_order[idx]) {
            message(sprintf("SNP %s reassigned from cluster %d to cluster %d (order improved from %d to %d)", snp_name, snp$cluster_index[idx], i, snp$cluster_order[idx], order))
            snp$cluster_index[idx] <- i
            snp$in_cluster_order[idx] <- order
          }
        } else {
          # First assignment
          snp$cluster_index[idx] <- i
          snp$in_cluster_order[idx] <- order
        }
      }
    }
  }

  # Return results
  cat("Done!\n")
  return(list(alpha = mat,
              variants = snp %>% arrange(desc(PIP)),
              models = result_df,
              sets = results$signal_cluster,
              elements = results$element_cluster,
              params = params))
}


#' Process matrix to keep only columns where the SNP ranks highest.
#' @param mat Input matrix
#' @import Rfast
#' @return Matrix with only first maximum value kept in each row
#' @noRd
process_matrix <- function(mat) {
  ranks <- Rfast::colRanks(mat, descending = TRUE)   # returns matrix with ranks by columns
  mask <- ranks == Rfast::rowMins(ranks, value = TRUE)
  mask <- col(mask) == apply(mask, 1, which.max)
  return(mat * mask)
}