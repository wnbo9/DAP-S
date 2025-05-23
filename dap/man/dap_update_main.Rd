% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dap_update_main}
\alias{dap_update_main}
\title{Implementation of updating DAP-S results algorithm in C++}
\usage{
dap_update_main(X, dap_result, prior_weights, r2_threshold, coverage)
}
\arguments{
\item{X}{Genotype matrix}

\item{dap_result}{DAP-S results}

\item{prior_weights}{Vector of prior probabilities}

\item{r2_threshold}{Threshold for LD}

\item{coverage}{Coverage for credible set}
}
\value{
A list containing:
\itemize{
  \item posterior_prob - Posterior probabilities
  \item log10_posterior_score - Log10 posterior scores
  \item log10_nc - Log10 normalizing constant
  \item pip - Posterior inclusion probabilities
  \item signal_cluster - Signal clusters
}
}
\description{
Implementation of updating DAP-S results algorithm in C++
}
