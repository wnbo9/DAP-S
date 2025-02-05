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
  mat * mask
}

#' Initialize parameters for fine mapping from SuSiE fit
#' @param fit SuSiE fit
#' @param exclusive Whether to create mutually exclusive signal clusters
#' @param use_susie_variance_estimate Use SuSiE-estimated prior variance
#' @param grid Grid of prior variances
#' @return List of parameters: matrix and scaled prior variance matrix
#' @noRd
param_setup <- function(fit, exclusive, use_susie_variance_estimate, grid) {
  # Extract the non-zero columns
  matrix <- t(fit$alpha)
  cols <- which(fit$V != 0)
  mat <- as.matrix(matrix[, cols])
  V <- fit$V[cols]

  # Remove duplicated columns
  max_indices <- apply(mat, 2, which.max)
  keep_cols <- !duplicated(max_indices)
  mat <- mat[, keep_cols, drop = FALSE]
  V <- V[keep_cols]

  # Process matrix if exclusive
  if (exclusive) {
    mat <- process_matrix(mat)
  }

  # Set up scaled prior variance matrix
  if (use_susie_variance_estimate) {
    phi2_mat <- matrix(V / as.vector(var(y)), nrow = 1)
  } else {
    #if (max(V) > 2) { # shift the upper bound of the grid
    #  grid[length(grid)] <- 1.28
    #}
    phi2_mat <- matrix(grid, ncol = ncol(mat), byrow = TRUE)
  }

  list(mat = mat, phi2_mat = phi2_mat)
}




#' Get summarization of fine mapping results
#' @param susie_params SuSiE parameters
#' @param info Information from SuSiE fit
#' @param results Fine mapping results
#' @param pir_threshold Threshold for PIR
#' @param r2_threshold Threshold for R2
#' @param prior_weights Prior weights
#' @param snp_names SNP names
#' @return Summarized results
#' @noRd
get_summarization <- function(susie_params, info, results,
                              pir_threshold, r2_threshold,
                              prior_weights, snp_names) {

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
    phi2 = susie_params$scaled_prior_variance,
    phi2_mat = info$phi2_mat,
    prior_weights = prior_weights,
    null_weight = susie_params$null_weight
  )

  # Initialize SNP dataframe
  snp_df <- data.frame(
    snp = snp_names,
    index = 1:p,
    marginal_pip = results$pip,
    cluster_index = "",
    stringsAsFactors = FALSE
  )

  sc <- results$signal_cluster
  snp_index <- sc$snp_index
  all_snps <- sort(unique(unlist(snp_index)))
  for (snp in all_snps) {
    clusters <- which(sapply(snp_index, function(x) snp %in% x))
    snp_df$cluster_index[snp_df$index == snp] <- paste(clusters, collapse = ",")
  }
  snp_df$cluster_index[snp_df$cluster_index == ""] <- "-1"
  snp_df$duplicate <- grepl(",", snp_df$cluster_index)

  new_df <- snp_df
  colnames(new_df)[colnames(new_df) == "marginal_pip"] <- "pip"
  dup_indices <- which(new_df$cluster_index != "-1")
  all_new_rows <- list()
  for (i in dup_indices) {
    clusters <- as.numeric(strsplit(new_df$cluster_index[i], ",")[[1]])
    new_rows <- lapply(clusters, function(cluster) {
      row <- new_df[i, ]
      row$cluster_index <- cluster
      effect_value <- results$effect_pip[i, sc$sc_index[cluster]]
      row$pip <- effect_value
      return(row)
    })
    all_new_rows <- c(all_new_rows, new_rows)
  }

  new_df <- new_df[-dup_indices, ]
  if (length(all_new_rows) > 0) {
    new_rows_df <- do.call(rbind, all_new_rows)
    new_df <- rbind(new_df, new_rows_df)
  }
  new_df$cluster_index <- as.numeric(new_df$cluster_index)

  # Return compiled results
  list(
    alpha = info$mat,
    alpha_dap = results$effect_pip,
    variants = snp_df,
    variants_coloc = new_df,
    models = result_df,
    sets = sc,
    elements = results$element_cluster,
    params = params
  )
}