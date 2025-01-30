#' DAP-S fine-mapping with individual-level data.
#'
#' The parameters are directly obtained from SuSiE results. Most of the commonly
#' used SuSiE parameters (e.g., \code{scaled_prior_variance}, \code{residual_variance},
#' etc.) can be controlled via \code{default_susie_params} and \code{susie_params}.
#'
#' @param X A matrix of genotype data (\eqn{n} x \eqn{p}), where \eqn{n} is the
#'   number of subjects and \eqn{p} is the number of variants.
#' @param y A numeric vector of phenotypes (length \eqn{n}).
#' @param L The number of potential causal variants (or “effects”) to include in
#'   the model. Defaults to \code{min(10, ncol(X))}.
#'
#' @param default_susie_params A named list of default SuSiE parameters (e.g.,
#'   \code{scaled_prior_variance}, \code{residual_variance}, \code{prior_weights},
#'   \code{null_weight}, \code{standardize}, \code{estimate_residual_variance},
#'   \code{estimate_prior_variance}). These specify the “factory” SuSiE settings.
#'   You may modify any of these if you want to change the defaults.
#'
#' @param susie_params A named list of parameters that override elements in
#'   \code{default_susie_params}. Only those specified here will be changed; the
#'   rest remain as in \code{default_susie_params}.
#'
#' @param scaled_prior_variance SuSiE input: scaled prior effect size variance.
#'   The prior variance is \code{var(y) * scaled_prior_variance}. Equivalently,
#'   it corresponds to \eqn{phi^2} in the BVSR with a \eqn{D_2} prior.
#'
#' @param residual_variance SuSiE input: residual variance. If
#'   \code{estimate_residual_variance = TRUE}, this value is used as the initial
#'   guess; otherwise, the residual variance is fixed to this value.
#'
#' @param prior_weights SuSiE input: optional prior weights on each variant;
#'   a length-\eqn{p} vector giving the prior probability each variant is causal.
#'
#' @param null_weight SuSiE input: a probability between 0 and 1 that no variant
#'   is causal.
#'
#' @param standardize SuSiE input: whether to standardize the genotype matrix.
#'   If \code{TRUE}, each column of \code{X} is standardized (mean 0, variance 1).
#'
#' @param estimate_residual_variance SuSiE input: if \code{TRUE}, the residual
#'   variance is estimated internally (using \code{residual_variance} as an initial
#'   value). If \code{FALSE}, it is held fixed.
#'
#' @param estimate_prior_variance SuSiE input: if \code{TRUE}, the prior variance
#'   is estimated (one parameter for all L effects) starting from
#'   \code{scaled_prior_variance}. If \code{FALSE}, the prior variance is fixed to
#'   \code{scaled_prior_variance}.
#'
#' @param r2_threshold A genotype \eqn{R^2} threshold for constructing signal
#'   clusters in linkage disequilibrium (LD). Default is \eqn{0.25}.
#'
#' @param coverage Coverage level for constructing credible sets (e.g., \eqn{0.95}).
#'   If specified, the function produces credible sets covering the specified
#'   probability. If not set, it returns signal clusters.
#'
#' @param use_susie_variance_estimate If \code{TRUE}, uses the SuSiE-estimated
#'   variance to compute the Bayes factor. If \code{FALSE}, a grid of values
#'   given by \code{grid} is used, and the Bayes factor is averaged across
#'   that grid.
#'
#' @param grid A numeric vector of scaled prior variances used only if
#'   \code{use_susie_variance_estimate = FALSE}. Defaults to
#'   \eqn{c(0.04, 0.16, 0.64)}.
#'
#' @param exclusive If \code{TRUE}, DAP-S enforces mutually exclusive signals
#'   (each variant can only be allocated to one signal). If \code{FALSE}, the
#'   same variant may be included in multiple signals (i.e., “splitting”).
#'
#' @param pir_threshold Pseudo-Importance Resampling (PIR) density threshold.
#'   Models with PIR density above this threshold are considered candidates.
#'   Defaults to \eqn{10^{-6}}.
#'
#' @import Rfast
#' @importFrom susieR susie
#' @importFrom dplyr %>% arrange desc
#'
#' @return A list containing the fine-mapping results, including posterior
#'   inclusion probabilities, credible sets (if \code{coverage} is specified), or
#'   signal clusters (if \code{coverage} is \code{NULL}). Typically, it includes
#'   the SuSiE fit object extended with DAP-S features.
#'
#' @export
dap <- function(X, y, L = min(10, ncol(X)),
                default_susie_params = list(
                  scaled_prior_variance = 0.2,
                  residual_variance = NULL,
                  prior_weights = NULL,
                  null_weight = NULL,
                  standardize = TRUE,
                  estimate_residual_variance = TRUE,
                  estimate_prior_variance = TRUE
                ),
                susie_params = list(),
                r2_threshold = 0.25,
                coverage = NULL,
                use_susie_variance_estimate = TRUE,
                grid = c(0.04, 0.16, 0.64),
                exclusive = TRUE,
                pir_threshold = 1e-6) {

  # Process inputs
  cat("Processing inputs...\n")
  params <- modifyList(default_susie_params, susie_params)
  X <- scale(X, scale = params$standardize)
  y <- scale(y, scale = FALSE)
  p <- ncol(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("SNP_", 1:p)

  ## Initialize parameters
  cat("Initializing parameters...\n")
  if (is.null(params$prior_weights)) {
    params$prior_weights <- rep(1 / p, p)
  }
  if (is.null(params$null_weight)) {
    params$null_weight <- prod(1 - params$prior_weights)
  }
  if (is.null(coverage)) {
    coverage <- 2
  }

  # Run SuSiE
  cat(ifelse(length(susie_params) == 0,
             "Running SuSiE with default settings...\n", "Running SuSiE...\n"))
  susie_fit <- susie(X, y, L, max_iter = 1000, coverage = 0.95,
                     scaled_prior_variance = params$scaled_prior_variance,
                     residual_variance = params$residual_variance,
                     prior_weights = params$prior_weights,
                     null_weight = params$null_weight,
                     standardize = params$standardize,
                     estimate_residual_variance = params$estimate_residual_variance,
                     estimate_prior_variance = params$estimate_prior_variance)
  # Extract matrix of sampling density from SuSiE results
  mat <- get_mat(susie_fit, exclusive)
  if (use_susie_variance_estimate) { # use SuSiE-estimated variance
    phi2_mat <- matrix(susie_fit$V[1:ncol(mat)] / as.vector(var(y)),
                       nrow = 1)
  } else { # use grid of values
    phi2_mat <- matrix(rep(grid, each = ncol(mat)),
                       nrow = length(grid), ncol = ncol(mat), byrow = TRUE)
  }

  cat("Running DAP-S fine-mapping...\n")
  # then we do PIR to get models
  results <- dap_main(X, y, mat, params$prior_weights,
                      r2_threshold, coverage, phi2_mat,
                      exclusive, pir_threshold)

  # summarize results
  cat("Summarizing results...\n")
  output <- get_summarization(results, mat,
                              pir_threshold, r2_threshold,
                              params$scaled_prior_variance, phi2_mat,
                              params$prior_weights, params$null_weight, colnames(X))

  # Return results
  cat("Done!\n")
  return(output)
}