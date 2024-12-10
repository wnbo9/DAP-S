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
#' @param option Option for DAP. 1: default SuSiE; 2: SuSiE with residual variance estimation; 3: SuSiE with prior variance estimation
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
                coverage = NULL,
                option = 1) {

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
  if (is.null(phi2_vec)) phi2_vec <- c(0.04, 0.16, 0.36, 0.64)
  if (is.null(residual_tau)) residual_tau <- 1/var(y)
  if (is.null(coverage)) coverage <- 10

  # Run SuSiE
  if (option == 1){
      susie_fit <- susie(X, y, L, max_iter = 1000, coverage = 0.95, standardize = FALSE,
               null_weight = null_weight, prior_weights = prior_weights)
  } else {
      susie_fit <- susie(X, y, L, max_iter = 1000, coverage = 0.95, standardize = FALSE,
               null_weight = null_weight, prior_weights = prior_weights,
               estimate_residual_variance = FALSE, residual_variance = 1/residual_tau,
               estimate_prior_variance = FALSE, scaled_prior_variance = phi2)
  }
  matrix <- t(susie_fit$alpha)

  # Run PIR
  results <- dap_main(X, y, matrix, threshold, prior_weights, phi2_vec, r2_threshold, coverage)


  # Create results dataframe
  result_df <- data.frame(
    Model = results$model_config,
    Posterior_Prob = results$posterior_prob,
    Log10_Posterior_Score = results$log10_posterior_score,
    stringsAsFactors = FALSE
  ) %>% arrange(desc(Log10_Posterior_Score))
  params <- list(
        log10_nc = results$log10_nc,
        threshold = threshold,
        phi2 = phi2,
        phi2_vec = phi2_vec,
        prior_weights = prior_weights,
        null_weight = null_weight,
        residual_tau = residual_tau
  )


  # Return results
  cat("Done!\n")
  return(list(
        alpha = matrix,
        pip = results$pip,
        models = result_df,
        sets = results$signal_cluster,
        params = params
  ))
}
