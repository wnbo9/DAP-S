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
#' @param y Phenotype vector
#' @param fit SuSiE fit
#' @param overlapping Whether to create overlapping signal clusters
#' @param use_susie_variance_estimate Use SuSiE-estimated prior variance
#' @param grid Grid of prior variances
#' @return List of parameters: matrix and scaled prior variance matrix
#' @noRd
param_setup <- function(y, fit, overlapping,
                        use_susie_variance_estimate, grid) {
  # Extract the non-zero columns
  matrix <- t(fit$alpha)
  cols <- which(fit$V != 0)
  if (length(cols) == 0) {
    mat <- as.matrix(matrix[, 1])
    V <- 0.2
  } else {
    mat <- as.matrix(matrix[, cols])
    V <- fit$V[cols]
  }

  # Remove duplicated columns caused by model overfitting
  max_indices <- apply(mat, 2, which.max)
  keep_cols <- !duplicated(max_indices) | max_indices == nrow(mat)
  mat <- mat[, keep_cols, drop = FALSE]
  V <- V[keep_cols]

  # Process matrix if not overlapping
  if (!overlapping) {
    mat <- process_matrix(mat)
  }

  # Set up scaled prior variance matrix
  if (use_susie_variance_estimate && max(V / as.vector(var(y))) > 1) {
    #phi2_mat <- matrix(V / as.vector(var(y)), nrow = 1)
    grid[length(grid)] <- 1.28 # shift the upper bound of the grid
  }
  phi2_mat <- matrix(rep(grid, ncol(mat)), ncol = ncol(mat), byrow = FALSE)

  list(mat = mat, phi2_mat = phi2_mat, overlapping = overlapping)
}




#' Get summarization of fine mapping results
#' @param susie_params SuSiE parameters
#' @param info Information from SuSiE fit
#' @param results Fine mapping results
#' @param pir_threshold Threshold for PIR
#' @param r2_threshold Threshold for R2
#' @param prior_weights Prior weights
#' @param snp_names SNP names
#' @param coverage Coverage for credible set
#' @param twas_weight Whether to compute TWAS weights
#' @return Summarized results
#' @noRd
get_summarization <- function(susie_params, info, results,
                              pir_threshold, r2_threshold,
                              prior_weights, snp_names,
                              coverage, twas_weight) {

  # Create model results dataframe
  p <- length(snp_names)

  model_df <- data.frame(
    model_config = results$model_config,
    posterior_prob = results$posterior_prob,
    log10_posterior_score = results$log10_posterior_score,
    log10_BF = results$log10_BF,
    log10_prior = results$log10_prior,
    stringsAsFactors = FALSE
  )

  # Store parameters
  params <- list(
    log10_nc = results$log10_nc,
    pir_threshold = pir_threshold,
    r2_threshold = r2_threshold,
    phi2 = susie_params$scaled_prior_variance,
    phi2_mat = info$phi2_mat,
    null_weight = susie_params$null_weight,
    twas_weight = twas_weight,
    overlapping = info$overlapping,
    coverage = coverage
  )

  sc <- results$signal_cluster
  cs_id <- list()
  snp_id <- list()
  pip_values <- list()

  # Loop through each cluster
  for (i in seq_along(sc$snp_index)) {
    current_cs_id <- sc$sc_index[i]  # Get cluster ID
    current_snps <- sc$snp_index[[i]]  # Get SNPs in this
    current_pips <- results$effect_pip[current_snps, current_cs_id]
    # Add to lists
    cs_id[[i]] <- rep(current_cs_id, length(current_snps))
    snp_id[[i]] <- current_snps
    pip_values[[i]] <- current_pips
  }

  # Combine into dataframe
  varaint_df <- data.frame(
    SC_ID = unlist(cs_id),
    SNP_ID = unlist(snp_id),
    PIP = unlist(pip_values),
    SNP_Name = snp_names[unlist(snp_id)],
    stringsAsFactors = FALSE
  )

  # Return compiled results
  list(
    alpha = info$mat,
    alpha_dap = results$effect_pip,
    pip = results$pip,
    prior = prior_weights,
    snp_names = snp_names,
    variants = varaint_df,
    models = model_df,
    model_combo = results$model_combo,
    sets = sc,
    params = params,
    reg_weights = results$reg_weights,
    twas_weights = results$twas_weights
  )
}