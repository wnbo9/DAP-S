#' Process matrix to keep only columns where the SNP ranks highest.
#' @param mat Input matrix
#' @import Rfast
#' @return Matrix with only first maximum value kept in each row
#' @noRd
process_matrix <- function(mat) {
  ranks <- Rfast::colRanks(mat, descending = TRUE)
  # returns matrix with ranks by columns
  mask <- ranks == Rfast::rowMins(ranks, value = TRUE)
  mask <- col(mask) == apply(mask, 1, which.max)
  return(mat * mask)
}

#' Extract matrix from SuSiE fit
#' @param susie_fit SuSiE fit
#' @param exclusive Whether to create mutually exclusive signal clusters
#' @return Matrix of sampling
#' @noRd
get_mat <- function(susie_fit, exclusive) {
  matrix <- t(susie_fit$alpha)
  # We only keep at most one column with null as top among all non-zero effects.
  # Note SuSiE can output two columns with similar ranks, we keep the first.
  mat <- as.matrix(matrix[, 1:sum(susie_fit$V != 0)])
  max_indices <- apply(mat, 2, which.max)
  # Keep first occurrence of each maximum position
  keep_cols <- !duplicated(max_indices)
  mat <- mat[, keep_cols, drop = FALSE]
  # If exclusive, we need to process the mat, so that each SNP only appears
  # in the effect where it has highest rank; otherwise, we keep all values.
  if (exclusive) {
    mat <- process_matrix(mat)
  }
  return(mat)
}


#' Get summarization of fine mapping results
#' @param results Fine mapping results
#' @param mat Matrix of sampling density
#' @param pir_threshold Threshold for PIR
#' @param r2_threshold Threshold for R2
#' @param scaled_prior_variance Scaled prior variance
#' @param phi2_mat Vector of prior variances
#' @param prior_weights Prior weights
#' @param null_weight Null weight
#' @param snp_names SNP names
#' @return Summarized results
#' @noRd
get_summarization <- function(results, mat,
                              pir_threshold, r2_threshold,
                              scaled_prior_variance, phi2_mat,
                              prior_weights, null_weight, snp_names) {
  # Create model results dataframe
  p <- length(snp_names)

  result_df <- data.frame(
    model_config = results$model_config,
    posterior_prob = results$posterior_prob,
    log10_posterior_score = results$log10_posterior_score,
    log10_BF = results$log10_BF,
    log10_prior = results$log10_prior,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(log10_posterior_score))

  # Store parameters
  params <- list(
    log10_nc = results$log10_nc,
    pir_threshold = pir_threshold,
    r2_threshold = r2_threshold,
    phi2 = scaled_prior_variance,
    phi2_mat = phi2_mat,
    prior_weights = prior_weights,
    null_weight = null_weight
  )

  # Initialize SNP dataframe
  snp <- data.frame(
    snp = snp_names,
    index = 1:p,
    pip = results$pip,
    cluster_index = -1,
    duplicate = 0,
    stringsAsFactors = FALSE
  )

  sc <- results$signal_cluster
  pip_matrix <- sc$pip_matrix
  snp_index <- sc$snp_index
  for(i in unique(unlist(snp_index))) {
    snp$cluster_index[i] <- which.max(pip_matrix[i, ])
    snp$duplicate[i] <- sum(pip_matrix[i, ] > 0) > 1
  }
  snp <- cbind(snp, pip_matrix)
  sc$pip_matrix <- NULL


  # Return compiled results
  list(
    alpha = mat,
    variants = snp,
    models = result_df,
    sets = sc,
    elements = results$element_cluster,
    params = params
  )
}