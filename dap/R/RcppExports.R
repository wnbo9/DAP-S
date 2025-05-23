# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compute log10 posterior scores for multiple model configurations
#' @param cmfg_matrix An m*p matrix of model configurations
#' @param pi_vec A vector of prior probabilities
#' @param phi2_mat A matrix of phi2 values
#' @param ss Flag indicating use of summary statistics (0 = use X,y; 1 = use summary statistics)
#' @param X_input An n*p matrix of genotype data (required when ss=0)
#' @param y_input An n-vector of phenotype data (required when ss=0)
#' @param XtX_input A p*p matrix of precomputed X'X (required when ss=1)
#' @param Xty_input A p*1 vector of precomputed X'y (required when ss=1)
#' @param yty_input A scalar representing y'y (required when ss=1)
#' @param n_input Sample size (required when ss=1)
#' @param twas_weight A boolean indicating whether to compute TWAS weights
#' @return An m-vector of log10 posterior scores
#' @export
compute_log10_posterior <- function(cmfg_matrix, pi_vec, phi2_mat, ss = 0L, X_input = NULL, y_input = NULL, XtX_input = NULL, Xty_input = NULL, yty_input = NULL, n_input = NULL, V_input = NULL, Dsq_input = NULL, var_input = NULL, XtOmegay_input = NULL, twas_weight = FALSE) {
    .Call(`_dap_compute_log10_posterior`, cmfg_matrix, pi_vec, phi2_mat, ss, X_input, y_input, XtX_input, Xty_input, yty_input, n_input, V_input, Dsq_input, var_input, XtOmegay_input, twas_weight)
}

#' Implementation of DAP-S algorithm in C++
#' @param X Genotype matrix
#' @param y Phenotype vector
#' @param matrix Proposal density matrix from SuSiE
#' @param pir_threshold Threshold for PIR
#' @param prior_weights Vector of prior probabilities
#' @param phi2_mat Matrix of scaled prior effect size variances
#' @param r2_threshold Threshold for LD
#' @param coverage Coverage for credible set
#' @param overlapping If TRUE. enforce overlapping clusters
#' @param twas_weight If TRUE, return TWAS weights
#' @param snp_names SNP names
#'
#' @return A list containing:
#' \itemize{
#'   \item model_config - Model configurations
#'   \item posterior_prob - Posterior probabilities
#'   \item log10_posterior_score - Log10 posterior scores
#'   \item log10_nc - Log10 normalizing constant
#'   \item pip - Posterior inclusion probabilities
#'   \item signal_cluster - Signal clusters
#' }
#' @export
dap_main <- function(matrix, pir_threshold, prior_weights, phi2_mat, r2_threshold, coverage, overlapping, twas_weight, snp_names, ss = 0L, X_input = NULL, y_input = NULL, XtX_input = NULL, Xty_input = NULL, yty_input = NULL, n_input = NULL, V_input = NULL, Dsq_input = NULL, var_input = NULL, XtOmegay_input = NULL) {
    .Call(`_dap_dap_main`, matrix, pir_threshold, prior_weights, phi2_mat, r2_threshold, coverage, overlapping, twas_weight, snp_names, ss, X_input, y_input, XtX_input, Xty_input, yty_input, n_input, V_input, Dsq_input, var_input, XtOmegay_input)
}

#' Implementation of updating DAP-S results algorithm in C++
#' @param X Genotype matrix
#' @param dap_result DAP-S results
#' @param prior_weights Vector of prior probabilities
#' @param r2_threshold Threshold for LD
#' @param coverage Coverage for credible set
#'
#' @return A list containing:
#' \itemize{
#'   \item posterior_prob - Posterior probabilities
#'   \item log10_posterior_score - Log10 posterior scores
#'   \item log10_nc - Log10 normalizing constant
#'   \item pip - Posterior inclusion probabilities
#'   \item signal_cluster - Signal clusters
#' }
dap_update_main <- function(X, dap_result, prior_weights, r2_threshold, coverage) {
    .Call(`_dap_dap_update_main`, X, dap_result, prior_weights, r2_threshold, coverage)
}

