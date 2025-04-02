#' Method of Moments (MoM) estimation for sigma^2 and tau^2
#'
#' @param PIP Matrix of posterior inclusion probabilities
#' @param mu Matrix of posterior means conditional on causal
#' @param omega Matrix of posterior precisions conditional on causal
#' @param sigmasq Current estimate of sigma^2
#' @param tausq Current estimate of tau^2
#' @param n Sample size
#' @param V Matrix of eigenvectors of X'X
#' @param Dsq Eigenvalues of X'X
#' @param VtXty V'X'y
#' @param Xty X'y
#' @param yty y'y
#' @param est_sigmasq Whether to estimate sigma^2
#' @param est_tausq Whether to estimate tau^2
#' @param verbose Whether to print update information
#'
#' @return List containing updated sigmasq and tausq
#' @keywords internal
MoM <- function(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
                est_sigmasq, est_tausq, verbose) {
  # Get dimensions
  p <- nrow(mu)
  L <- ncol(mu)

  # Compute A
  A <- matrix(c(n, sum(Dsq), sum(Dsq), sum(Dsq^2)), nrow = 2)

  # Compute diag(V'MV)
  b <- rowSums(mu * PIP)
  Vtb <- t(V) %*% b
  diagVtMV <- Vtb^2
  tmpD <- numeric(p)

  for (l in 1:L) {
    bl <- mu[, l] * PIP[, l]
    Vtbl <- t(V) %*% bl
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + PIP[, l] * (mu[, l]^2 + 1/omega[, l])
  }

  diagVtMV <- diagVtMV + colSums(V^2 * tmpD)

  # Compute x
  x <- numeric(2)
  x[1] <- yty - 2 * sum(b * Xty) + sum(Dsq * diagVtMV)
  x[2] <- sum(Xty^2) - 2 * sum(Vtb * VtXty * Dsq) + sum(Dsq^2 * diagVtMV)

  if (est_tausq) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigmasq <- sol[1]
      tausq <- sol[2]
    } else {
      sigmasq <- x[1]/n
      tausq <- 0
    }

    if (verbose) {
      cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigmasq, tausq))
    }
  } else if (est_sigmasq) {
    sigmasq <- (x[1] - A[1, 2] * tausq) / n
    if (verbose) {
      cat(sprintf("Update sigma^2 to %f\n", sigmasq))
    }
  }

  return(list(sigmasq = sigmasq, tausq = tausq))
}

#' Maximum Likelihood Estimation (MLE) for sigma^2 and tau^2
#'
#' @param PIP Matrix of posterior inclusion probabilities
#' @param mu Matrix of posterior means conditional on causal
#' @param omega Matrix of posterior precisions conditional on causal
#' @param sigmasq Current estimate of sigma^2
#' @param tausq Current estimate of tau^2
#' @param n Sample size
#' @param V Matrix of eigenvectors of X'X
#' @param Dsq Eigenvalues of X'X
#' @param VtXty V'X'y
#' @param yty y'y
#' @param est_sigmasq Whether to estimate sigma^2
#' @param est_tausq Whether to estimate tau^2
#' @param sigmasq_range Range for sigma^2 optimization
#' @param tausq_range Range for tau^2 optimization
#' @param it Current iteration
#' @param verbose Whether to print update information
#'
#' @return List containing updated sigmasq and tausq
#' @keywords internal
MLE <- function(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, yty,
                est_sigmasq, est_tausq, sigmasq_range, tausq_range, it, verbose) {
  # Get dimensions
  p <- nrow(mu)
  L <- ncol(mu)

  # Set default ranges if NULL
  if (is.null(sigmasq_range)) sigmasq_range <- c(0.2 * yty/n, 1.2 * yty/n)
  if (is.null(tausq_range)) tausq_range <- c(1e-12, 1.2 * yty/(n * p))

  # Compute diag(V'MV)
  b <- rowSums(mu * PIP)
  Vtb <- t(V) %*% b
  diagVtMV <- Vtb^2
  tmpD <- numeric(p)

  for (l in 1:L) {
    bl <- mu[, l] * PIP[, l]
    Vtbl <- t(V) %*% bl
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + PIP[, l] * (mu[, l]^2 + 1/omega[, l])
  }

  diagVtMV <- diagVtMV + colSums(V^2 * tmpD)

  # Define negative ELBO as function of x = (sigma_e^2, sigma_g^2)
  f <- function(x) {
    sigmasq <- x[1]
    tausq <- x[2]

    0.5 * (n - p) * log(sigmasq) + 0.5/sigmasq * yty +
      sum(0.5 * log(tausq * Dsq + sigmasq) -
            0.5 * tausq/sigmasq * VtXty^2/(tausq * Dsq + sigmasq) -
            Vtb * VtXty/(tausq * Dsq + sigmasq) +
            0.5 * Dsq/(tausq * Dsq + sigmasq) * diagVtMV)
  }

  if (est_tausq) {
    res <- optim(c(sigmasq, tausq), f, method = "L-BFGS-B",
                 lower = c(sigmasq_range[1], tausq_range[1]),
                 upper = c(sigmasq_range[2], tausq_range[2]))

    if (res$convergence == 0) {
      sigmasq <- res$par[1]
      tausq <- res$par[2]

      if (verbose) {
        cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigmasq, tausq))
      }
    } else {
      warning(sprintf("sigma^2 and tau^2 update for iteration %d failed to converge; keeping previous parameters", it))
    }
  } else if (est_sigmasq) {
    g <- function(x) f(c(x, tausq))

    res <- optim(sigmasq, g, method = "L-BFGS-B",
                 lower = sigmasq_range[1], upper = sigmasq_range[2])

    if (res$convergence == 0) {
      sigmasq <- res$par

      if (verbose) {
        cat(sprintf("Update sigma^2 to %f\n", sigmasq))
      }
    } else {
      warning(sprintf("sigma^2 update for iteration %d failed to converge; keeping previous parameters", it))
    }
  }

  return(list(sigmasq = sigmasq, tausq = tausq))
}

#' SuSiE with random effects
#'
#' @param z Vector of z-scores (equal to X'y/sqrt(n))
#' @param bhat Beta hat
#' @param shat Standard error hat
#' @param var_y Variance of y, equal to ||y||^2/(n-1)
#' @param n Sample size
#' @param L Number of modeled causal effects
#' @param LD LD matrix (equal to X'X/n)
#' @param V Precomputed p x p matrix of eigenvectors of X'X
#' @param Dsq Precomputed length-p vector of eigenvalues of X'X
#' @param est_ssq Estimate prior effect size variances s^2 using MLE
#' @param ssq Length-L initialization s^2 for each effect
#' @param ssq_range Lower and upper bounds for each s^2, if estimated
#' @param prior_weights Length-p vector of prior causal probability for each SNP; must sum to 1
#' @param null_weight A value specify the probability of no effect; should be between 0 and 1.
#' @param est_sigmasq Estimate variance sigma^2
#' @param est_tausq Estimate both variances sigma^2 and tau^2
#' @param sigmasq Initial value for sigma^2
#' @param tausq Initial value for tau^2
#' @param method One of c('moments','MLE')
#' @param sigmasq_range Lower and upper bounds for sigma^2, if estimated using MLE
#' @param tausq_range Lower and upper bounds for tau^2, if estimated using MLE
#' @param PIP p x L initializations of PIPs
#' @param mu p x L initializations of mu
#' @param maxiter Maximum number of SuSiE iterations
#' @param PIP_tol Convergence threshold for PIP difference between iterations
#' @param verbose Whether to print update information
#'
#' @return List with results
#' @export
susie_inf <- function(z, bhat, shat, var_y, n, L, LD = NULL, V = NULL, Dsq = NULL,
                      est_ssq = TRUE, ssq = NULL, ssq_range = c(0, 1), prior_weights = NULL, null_weight = 0,
                      est_sigmasq = TRUE, est_tausq = TRUE, sigmasq = 1, tausq = 0,
                      method = "moments", sigmasq_range = NULL, tausq_range = NULL,
                      PIP = NULL, mu = NULL, maxiter = 100, PIP_tol = 1e-3, verbose = TRUE) {

  if (missing(z)) {
    p <- length(bhat)
    z = bhat / shat
  } else {
    p <- length(z)
  }

  if (is.numeric(null_weight) && null_weight == 0) null_weight <- NULL
  # Initialize prior causal probabilities
  if (is.null(null_weight)) {
    if (is.null(prior_weights)) {
      logpi <- rep(log(1.0/p), p)
    } else {
      logpi <- rep(-Inf, p)
      inds <- which(prior_weights > 0)
      logpi[inds] <- log(prior_weights[inds])
    }
  } else {
    if (is.null(prior_weights)) {
      logpi <- log(c(rep(1.0/p * (1-null_weight), p), null_weight))
    } else {
      prior_weights <- c(prior_weights * (1-null_weight), null_weight)
      logpi <- rep(-Inf, p+1)
      inds <- which(prior_weights > 0)
      logpi[inds] <- log(prior_weights[inds])
    }
    LD <- cbind(rbind(LD,0), 0)
    z <- c(z, 0)
    p <- p+1
  }

  adj = (n-1)/(z^2 + n - 2)
  z = sqrt(adj) * z

  # Precompute V, D^2 in the SVD X = UDV', and V'X'y and y'y
  if ((is.null(V) || is.null(Dsq)) && is.null(LD)) {
    stop("Missing LD")
  } else if (is.null(V) || is.null(Dsq)) {
    if (!missing(shat) & !missing(var_y)) {
      XtXdiag <- suppressWarnings(var_y * adj/(shat^2))
      XtX <- t(LD * sqrt(XtXdiag)) * sqrt(XtXdiag)
      XtX <- (XtX + t(XtX))/2
      R <- XtX/(n-1)
      eig <- eigen(R, symmetric = TRUE)
    } else {
      eig <- eigen(LD, symmetric = TRUE)
    }
    V <- eig$vectors
    Dsq <- pmax( (n-1) * eig$values, 0 )
  } else {
    Dsq <- pmax(Dsq, 0)
  }

  if (!missing(shat) & !missing(var_y)) {
    # on original scale
    Xty = suppressWarnings(z * sqrt(adj) * var_y / shat)
    VtXty <- t(V) %*% Xty
    yty <- (n-1) * var_y
  } else {
    # on standardized x and y
    Xty <- sqrt(n-1) * z
    VtXty <- t(V) %*% Xty
    yty <- (n-1)
  }


  # Initialize s_l^2, PIP_j, mu_j, omega_j
  if (is.null(ssq)) ssq <- rep(0.2, L)
  if (is.null(PIP)) PIP <- matrix(1/p, nrow = p, ncol = L)
  if (is.null(mu)) mu <- matrix(0, nrow = p, ncol = L)

  # Initialize diagonal variances, diag(X' Omega X), X' Omega y
  var <- tausq * Dsq + sigmasq
  diagXtOmegaX <- rowSums(sweep(V^2, 2, Dsq/var, "*"))
  XtOmegay <- V %*% (VtXty/var)

  lbf_variable <- matrix(0, nrow = p, ncol = L)
  lbf <- numeric(L)
  omega <- outer(diagXtOmegaX, 1/ssq, "+")


  # Main SuSiE iteration loop
  for (it in 1:maxiter) {
    if (verbose) cat(sprintf("Iteration %d\n", it))
    PIP_prev <- PIP

    # Single effect regression for each effect l = 1,...,L
    for (l in 1:L) {
      # Compute X' Omega r_l for residual r_l
      b <- rowSums(mu * PIP) - mu[, l] * PIP[, l]
      XtOmegaXb <- V %*% (t(V) %*% b * Dsq/var)
      XtOmegar <- XtOmegay - XtOmegaXb

      if (est_ssq) {
        # Update prior variance ssq[l]
        f <- function(x) {
          -log(sum(exp(-0.5 * log(1 + x * diagXtOmegaX) +
                         x * XtOmegar^2/(2 * (1 + x * diagXtOmegaX)) +
                         logpi)))
        }

        res <- optim(par = mean(ssq_range), fn = f, method = "Brent", lower = ssq_range[1], upper = ssq_range[2])

        if (res$convergence == 0) {
          ssq[l] <- res$par

          if (verbose) {
            cat(sprintf("Update s^2 for effect %d to %f\n", l, ssq[l]))
          }
        } else {
          cat(sprintf("WARNING: s^2 update for iteration %d, effect %d failed to converge; keeping previous parameters\n", it, l))
        }
      }

      # Update omega, mu, and PIP
      omega[, l] <- diagXtOmegaX + 1/ssq[l]
      mu[, l] <- XtOmegar/omega[, l]
      lbf_variable[, l] <- XtOmegar^2/(2 * omega[, l]) - 0.5 * log(omega[, l] * ssq[l])
      logPIP <- lbf_variable[, l] + logpi
      lbf[l] <- log(sum(exp(logPIP)))
      PIP[, l] <- exp(logPIP - lbf[l])
    }

    # Update variance components
    if (est_sigmasq || est_tausq) {
      if (method == "moments") {
        result <- MoM(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
                      est_sigmasq, est_tausq, verbose)
        sigmasq <- result$sigmasq
        tausq <- result$tausq
      } else if (method == "MLE") {
        result <- MLE(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, yty,
                      est_sigmasq, est_tausq, sigmasq_range, tausq_range, it, verbose)
        sigmasq <- result$sigmasq
        tausq <- result$tausq
      } else {
        stop("Unsupported variance estimation method")
      }

      # Update X' Omega X, X' Omega y
      var <- tausq * Dsq + sigmasq
      diagXtOmegaX <- rowSums(sweep(V^2, 2, Dsq/var, "*"))
      XtOmegay <- V %*% (VtXty/var)
    }

    # Determine convergence from PIP differences
    PIP_diff <- max(abs(PIP_prev - PIP))
    if (verbose) cat(sprintf("Maximum change in PIP: %f\n", PIP_diff))
    if (PIP_diff < PIP_tol) {
      if (verbose) cat("CONVERGED\n")
      break
    }
  }

  # Compute posterior means of b and alpha
  b <- rowSums(mu * PIP)
  XtOmegaXb <- V %*% (t(V) %*% b * Dsq/var)
  XtOmegar <- XtOmegay - XtOmegaXb
  alpha <- tausq * XtOmegar

  return(list(
    PIP = PIP,
    spip = 1 - apply(1 - PIP, 1, prod),
    mu = mu,
    omega = omega,
    lbf = lbf,
    lbf_variable = lbf_variable,
    ssq = ssq,
    sigmasq = sigmasq,
    tausq = tausq,
    alpha = alpha
  ))
}

#' Compute credible sets from single-effect PIPs
#'
#' @param PIP p x L PIP matrix output by susie
#' @param coverage Coverage level for each credible set
#' @param purity Sets with minimum absolute correlation < purity are removed
#' @param LD LD matrix (equal to X'X/n)
#' @param V Precomputed p x p matrix of eigenvectors of X'X
#' @param Dsq Precomputed length-p vector of eigenvalues of X'X
#' @param n Sample size
#' @param dedup Remove duplicate CS's
#'
#' @return List of variable indices corresponding to credible sets
#' @export
cred <- function(PIP, coverage = 0.9, purity = 0.5, LD = NULL, V = NULL, Dsq = NULL,
                 n = NULL, dedup = TRUE) {

  if ((is.null(V) || is.null(Dsq) || is.null(n)) && is.null(LD)) {
    stop("Missing inputs for purity filtering")
  }

  # Compute credible sets
  cred_list <- list()

  for (l in 1:ncol(PIP)) {
    sortinds <- order(PIP[, l], decreasing = TRUE)
    cumsum_pip <- cumsum(PIP[sortinds, l])
    ind <- min(which(cumsum_pip >= coverage))
    credset <- sortinds[1:(ind)]

    # Filter by purity
    if (length(credset) == 1) {
      cred_list[[length(cred_list) + 1]] <- credset
      next
    }

    if (length(credset) < 100) {
      rows <- credset
    } else {
      set.seed(123)
      rows <- sample(credset, size = 100, replace = FALSE)
    }

    if (!is.null(LD)) {
      LDloc <- LD[rows, rows, drop = FALSE]
    } else {
      LDloc <- t(V[rows, , drop = FALSE] * Dsq) %*% V[rows, , drop = FALSE] / n
    }

    if (min(abs(LDloc)) > purity) {
      cred_list[[length(cred_list) + 1]] <- sort(credset)
    }
  }

  if (dedup) {
    # Remove duplicates while maintaining order
    cred_list_str <- sapply(cred_list, function(x) paste(x, collapse = ","))
    unique_indices <- !duplicated(cred_list_str)
    cred_list <- cred_list[unique_indices]
  }

  return(cred_list)
}
