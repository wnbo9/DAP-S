#' DAP fine mapping
#' @param X Genotype data
#' @param y Phenotype data
#' @param L Number of causal variants
#' @param prior_weights Prior weights
#' @param null_weight Null weight
#' @param residual_tau Residual variance
#' @param threshold Threshold for proposal density
#' @param phi2 Scaled prior effect size variance
#' @param phi2_vec Scaled prior effect size variance vector
#' @param r2_threshold Genotype R2 threshold for LD
#' @param coverage Coverage of credible sets. When not set, it will output signal clusters; otherwise, it outputs credible sets at the specified coverage level
#'
#' @import Rfast
#' @importFrom susieR susie
#' @importFrom dplyr %>% arrange
#' @return Fine mapping results
#'
#' @export
dap <- function(X, y,
                L = min(10, ncol(X)),
                prior_weights = NULL,
                null_weight = NULL,
                residual_tau = NULL,
                threshold = 1e-6,
                phi2 = 0.36,
                phi2_vec = NULL,
                r2_threshold = 0.25,
                coverage = NULL) {

  # Process inputs
  cat("Processing inputs...\n")
  X <- scale(X, scale = FALSE)
  y <- scale(y, scale = FALSE)
  n <- length(y)
  p <- ncol(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("SNP_",1:p)


  ## Initialize parameters
  if (is.null(null_weight)) null_weight <- (1-1/p)^p
  if (is.null(prior_weights)) prior_weights <- rep(1/p, p)
  if (is.null(phi2_vec)) phi2_vec <- c(0.04, 0.16, 0.64)
  if (is.null(residual_tau)) residual_tau <- 1/var(y)
  if (is.null(coverage)) coverage <- 10

  # Run SuSiE
  susie_fit <- susie(X, y, L, max_iter = 1000, coverage = 0.95, standardize = FALSE,
               null_weight = null_weight, prior_weights = prior_weights,
               estimate_residual_variance = FALSE, residual_variance = 1/residual_tau,
               estimate_prior_variance = FALSE, scaled_prior_variance = phi2)
  matrix <- t(susie_fit$alpha)

  # Run PIR
  results <- dap_main(X, y, matrix, threshold, prior_weights, phi2_vec, r2_threshold, coverage)


  # Create model results
  result_df <- data.frame(
    Model = results$model_config,
    Posterior_Prob = results$posterior_prob,
    Log10_Posterior_Score = results$log10_posterior_score,
    Log10_BF = results$log10_BF,
    Log10_Prior = results$log10_prior,
    stringsAsFactors = FALSE
  ) %>% arrange(desc(Log10_Posterior_Score))
  params <- list(
        log10_nc = results$log10_nc,
        threshold = threshold,
        r2_threshold = r2_threshold,
        phi2 = phi2,
        phi2_vec = phi2_vec,
        prior_weights = prior_weights,
        null_weight = null_weight,
        residual_tau = residual_tau
  )
  snp <- data.frame(
      SNP = colnames(X),
      PIP = results$pip,
      cluster_index = -1,
      in_cluster_order = Inf,
      stringsAsFactors = FALSE
  )
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
  return(list(
        alpha = matrix,
        variants = snp %>% arrange(desc(PIP)),
        models = result_df,
        sets = results$signal_cluster,
        params = params
  ))
}
