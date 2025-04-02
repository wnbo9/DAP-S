#' DAP-S fine-mapping using summary statistics data
#' @param z The z-score matrix
#' @param R The correlation matrix
#' @param n The sample size
#' @param bhat The effect size estimates
#' @param shat
#' @param var_y The variance of the phenotype
#' @param estimate_residual_variance Whether to estimate the residual variance
#' @return The residual sum of squares
#' @export
#'
dap_rss <- function(z, R, n, bhat, shat, var_y,
                    estimate_residual_variance = FALSE,
                    prior_variance = 50) {

  if (missing(z)) {
    p <- length(bhat)
    z <- bhat / shat
  } else {
    p <- length(z)
  }

  # When n is provided, compute the PVE-adjusted z-scores
  if (!missing(n)) {
    adj <- (n-1) / (z^2 + n-2)
    z <- sqrt(adj) * z
  }


  if (missing(n)) {
    # When n is not provided, use unadjusted z-scores
    rst <- dap_suff_stat(XtX = R, Xty = z, n = 2, yty = 1,
                         susie_params = list(scaled_prior_variance = prior_variance,
                                            estimate_residual_variance = estimate_residual_variance),
                         standardize = FALSE)
  } else {
    # When n is provided, use PVE-adjusted z-scores
    if (!missing(shat) & !missing(var_y)) {
      # var_y, shat (and bhat) are provided, so the effects are on the *original scale*.
      XtXdiag <- var_y * adj/(shat^2)
      XtX <- t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
      XtX <- (XtX + t(XtX))/2
      Xty <- z * sqrt(adj) * var_y / shat
    } else {

      # the effects are on the standardized X, y scale
      XtX <- (n-1) * R
      Xty <- sqrt(n-1) * z
      var_y <- 1
    }

    rst <- dap_suff_stat(XtX = XtX, Xty = Xty, n = n, yty = (n-1) * var_y,
                        susie_params = list(estimate_residual_variance = estimate_residual_variance),
                        standardize = FALSE)
  }

  return(rst)
}
