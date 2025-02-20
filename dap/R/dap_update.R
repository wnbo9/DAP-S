#' Update DAP-S fine-mapping results.
#' @param X Genotype matrix (n x p)
#' @param dap_result DAP-S result from \code{dap}
#' @param prior_weights New prior weights of each SNP being causal
#' @param r2_threshold Threshold for R2 when constructing credible sets
#' @param coverage Coverage for credible set
#' @import Rfast
#' @importFrom susieR susie
#' @importFrom dplyr %>% arrange desc
#' @return List of fine-mapping results.
#' @export
dap_update <- function(X, dap_result,
                       prior_weights = NULL,
                       r2_threshold = 0.25,
                       coverage = NULL) {

  cat("Processing inputs...\n")

  # Validate inputs
  if (is.null(dap_result)) {
    stop("DAP result object is required")
  }

  # Use original coverage if not specified
  if (is.null(coverage)) {
    coverage <- dap_result$params$coverage
  }

  # Use original prior weights if not specified
  if (is.null(prior_weights)) {
    prior_weights <- dap_result$params$prior_weights
  }

  susie_params <- list(
    scaled_prior_variance = 0.2,
    residual_variance = NULL,
    null_weight = prod(1 - prior_weights)
  )
  p <- length(dap_result$snp_names)

  cat("Updating DAP-S fine-mapping results...\n")
  results <- dap_update_main(X, dap_result,
                             prior_weights,
                             r2_threshold,
                             coverage)

  cat("Summarizing results...\n")
  new_df <- data.frame(
    snp = dap_result$snp_names,
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
    new_df$cluster_index[new_df$index == snp] <- paste(clusters, collapse = ",")
  }
  new_df$cluster_index[new_df$cluster_index == ""] <- "-1"
  new_df$duplicate <- grepl(",", new_df$cluster_index)

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


  dap_result$pip <- results$pip
  dap_result$variants_coloc <- new_df
  dap_result$models$posterior_prob <- results$posterior_prob
  dap_result$models$log10_posterior_score <- results$log10_posterior_score
  dap_result$models$log10_prior <- results$log10_prior
  dap_result$sets <- results$signal_cluster
  dap_result$params$log10_nc <- results$log10_nc
  dap_result$params$prior_weights <- prior_weights
  dap_result$reg_weights <- results$reg_weights
  dap_result$twas_weights <- results$twas_weights

  cat("Done!\n")
  return(dap_result)
}