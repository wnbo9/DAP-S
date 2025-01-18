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
#' @param coverage Coverage of credible sets. When not set, it outputs signal clusters; otherwise, it outputs credible sets at the specified coverage level
#'
#' @import Rfast
#' @importFrom susieR susie_suff_stat
#' @importFrom dplyr %>% arrange
#' @return Fine mapping results
#'
#' @export
dap_suff_stat <- function(XtX, Xty, yty, n,
                L = min(10, ncol(X))) {
  print(class(X))
}
