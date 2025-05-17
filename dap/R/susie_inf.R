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
MoM <- function(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
                est_sigmasq, est_tausq, verbose) {
  # Get dimensions
  p <- nrow(mu)
  L <- ncol(mu)

  # Compute A
  A <- matrix(c(n, sum(Dsq), sum(Dsq), sum(Dsq^2)), nrow = 2)
  b <- rowSums(mu * PIP)
  Vtb <- t(V) %*% b
  diagVtMV <- Vtb^2
  tmpD <- numeric(p)

  for (l in 1:L) {
    bl <- mu[, l] * PIP[, l]
    Vtbl <- t(V) %*% bl
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
  }

  diagVtMV <- diagVtMV + rowSums(sweep(t(V)^2, 2, tmpD, `*`))

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

  diagVtMV <- diagVtMV + rowSums(sweep(t(V)^2, 2, tmpD, `*`))

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
                 upper = c(sigmasq_range[2], tausq_range[2]),
                 control = list(ndeps = rep(1e-8, 2), lmm = 10,
                                pgtol = 1e-5, factr = 1e7, maxit = 15000))

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

    res <- optim(sigmasq, g, method = "Brent",
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
#' @param bhat Vector of effect size estimates
#' @param shat Vector of standard errors
#' @param z Vector of z-scores
#' @param var_y Variance of y
#' @param n Sample size
#' @param L Number of modeled causal effects
#' @param LD LD matrix
#' @param V Precomputed p x p matrix of eigenvectors of X'X
#' @param Dsq Precomputed length-p vector of eigenvalues of X'X
#' @param null_weight Null weight
#' @param est_ssq Estimate prior effect size variances s^2 using MLE
#' @param ssq Length-L initialization s^2 for each effect
#' @param ssq_range Lower and upper bounds for each s^2, if estimated
#' @param pi Length-p vector of prior causal probability for each SNP; must sum to 1
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
#'
#' @examples
#' # Set random seed for reproducibility
#' set.seed(2025)
#' n <- 5000  # Number of samples
#' p <- 500   # Number of SNPs
#' MAF <- 0.1  # Minor allele frequency
#' Ltrue <- 5  # Number of true causal variants
#' ssq <- 0.01  # Effect size variance
#' sigmasq <- 1  # Error variance
#'
#' # Generate genotype matrix X
#' X <- matrix(rbinom(n*p, 2, MAF), nrow=n, ncol=p)
#' X <- scale(X, center=TRUE, scale=FALSE)
#'
#' # Generate sparse effects
#' b <- rep(0, p)
#' inds <- sample(1:p, size=Ltrue, replace=FALSE)
#' b[inds] <- rnorm(Ltrue) * sqrt(ssq)
#' order_idx <- order(inds)
#' cat('True effects:', inds[order_idx], '\n')
#' cat('Effect sizes:', b[inds[order_idx]], '\n')
#'
#' # Generate infinitesimal effects
#' tausq <- 1e-3
#' infinitesimal <- X %*% rnorm(p) * sqrt(tausq)
#' effects <- X %*% b + infinitesimal
#' y <- effects + rnorm(n) * sqrt(sigmasq)
#' y <- scale(y, center = TRUE, scale = FALSE)
#' cat('Total fraction of variance explained by SNPs:', var(as.vector(effects))/var(y), '\n')
#'
#' # Generate summary statistics
#' res = susieR::univariate_regression(X, y)
#' suff_stat <- susieR::compute_suff_stat(X, as.vector(y), standardize = FALSE)
#' LD <- cov2cor(suff_stat$XtX)
#' z <- res$betahat/res$sebetahat
#'
#' # Work on original scale without estimating infinitesimal effects
#' output <- susie_inf(bhat = res$betahat, shat = res$sebetahat, var_y = var(y), est_tausq = FALSE, n=n, L=5, LD = LD, null_weight = 0.3)
#' output$sigmasq
#' output$tausq
#'
#' # Work on original scale with estimating infinitesimal effects
#' output2 <- susie_inf(bhat = res$betahat, shat = res$sebetahat, var_y = var(y), n=n, L=5, LD = LD, null_weight = 0.3)
#' output2$sigmasq
#' output2$tausq
#'
#' plot(output2$spip, output$spip, xlab = "PIP of SuSiE-inf", ylab = "PIP of SuSiE")
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' # Work on standardized scale
#' output3 <- susie_inf(z = z, n = n, L = 5, LD = LD, null_weight = 0.3)
#' output3$sigmasq
#' output3$tausq
susie_inf <- function(bhat = NULL, shat = NULL, z = NULL, var_y = NULL, n, L = 10,
                      LD = NULL, V = NULL, Dsq = NULL, null_weight = 0,
                      est_ssq = TRUE, ssq = NULL, ssq_range = c(0, 1), pi = NULL,
                      est_sigmasq = TRUE, est_tausq = TRUE, sigmasq = 1, tausq = 0,
                      method = "moments", sigmasq_range = NULL, tausq_range = NULL,
                      PIP = NULL, mu = NULL, maxiter = 100, PIP_tol = 1e-3, verbose = FALSE) {

  suppressWarnings({
    if (is.null(z)) z <- bhat / shat
    p <- length(z)

    if (is.numeric(null_weight) && null_weight == 0) null_weight <- NULL

    if (is.null(null_weight)) {
      if (is.null(pi)) {
        logpi <- rep(log(1.0/p), p)
      } else {
        logpi <- rep(-Inf, p)
        inds <- which(pi > 0)
        logpi[inds] <- log(pi[inds])
      }
    } else {
      if (is.null(pi)) {
        logpi <- log(c(rep(1.0/p * (1-null_weight), p), null_weight))
      } else {
        pi <- c(pi * (1-null_weight), null_weight)
        logpi <- rep(-Inf, p+1)
        inds <- which(pi > 0)
        logpi[inds] <- log(pi[inds])
      }
      LD <- cbind(rbind(LD,0), 0)
      z <- c(z, 0)
      p <- p+1
    }

    adj <- (n-1)/(z^2 + n - 2)
    z   <- sqrt(adj) * z


    # Precompute V, D^2 in the SVD X = UDV', and V'X'y and y'y
    if ((is.null(V) || is.null(Dsq)) && is.null(LD)) {
      stop("Missing LD")
    } else if (is.null(V) || is.null(Dsq)) {
      if (!is.null(shat) && !is.null(var_y)) {
        if (!is.null(null_weight)) shat <- c(shat, Inf)
        XtXdiag <- var_y * adj/(shat^2)
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

    if (!is.null(shat) && !is.null(var_y)) {
      Xty <- z * sqrt(adj) * var_y / shat
      VtXty <- t(V) %*% Xty
      yty <- (n-1) * var_y
    } else {
      # on standardized x and y
      Xty <- sqrt(n-1) * z
      VtXty <- t(V) %*% Xty
      yty <- (n-1)
    }

    # Initialize diagonal variances, diag(X' Omega X), X' Omega y
    var <- tausq * Dsq + sigmasq
    diagXtOmegaX <- rowSums(sweep(V^2, 2, Dsq/var, "*"))
    XtOmegay <- V %*% (VtXty/var)

    # Initialize s_l^2, PIP_j, mu_j, omega_j
    if (is.null(ssq)) ssq <- rep(0.2, L)
    if (is.null(PIP)) PIP <- matrix(1/p, nrow = p, ncol = L)
    if (is.null(mu)) mu <- matrix(0, nrow = p, ncol = L)

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
          f <- function(logx) {
            x <- exp(logx)
            -log(sum(exp(-0.5 * log(1 + x * diagXtOmegaX) +
                           x * XtOmegar^2/(2 * (1 + x * diagXtOmegaX)) +
                           logpi)))
          }
          res <- optim(par = log(max(c(ssq[l], ssq_range[1], 1e-10))),
                       fn = f,
                       method = "Brent",
                       lower = log(max((ssq_range[1]), 1e-10)),
                       upper = log(ssq_range[2]))

          f0 <- -log(sum(exp(logpi)))  # equivalent to f(0)
          if (f0 < res$value) {
            ssq[l] <- 0
          } else {
            ssq[l] <- exp(res$par)
          }

          if (verbose) {
            cat(sprintf("Update s^2 for effect %d to %f\n", l, ssq[l]))
          }
        }

        # Update omega, mu, and PIP
        omega[, l] <- diagXtOmegaX + 1/ssq[l]
        mu[, l] <- XtOmegar/omega[, l]

        if (ssq[l] > 0) {
          lbf_variable[, l] <- XtOmegar^2/(2 * omega[, l]) - 0.5 * log(omega[, l] * ssq[l])
        } else {
          lbf_variable[, l] <- 0
        }

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
    spip <- 1 - apply(1 - PIP, 1, prod)
    if (!is.null(null_weight)) spip <- spip[-p]

    return(list(
      PIP = PIP,
      spip = spip,
      mu = mu,
      omega = omega,
      lbf = lbf,
      lbf_variable = lbf_variable,
      ssq = ssq,
      sigmasq = sigmasq,
      tausq = tausq,
      alpha = alpha,
      V = V,
      Dsq = Dsq,
      var = var,
      XtOmegay = XtOmegay
    ))
  })
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
  
  cred_list <- list()
  pip_sums <- c()  # store total PIP for each credible set
  
  for (l in 1:ncol(PIP)) {
    sortinds <- order(PIP[, l], decreasing = TRUE)
    cumsum_pip <- cumsum(PIP[sortinds, l])
    ind_candidates <- which(cumsum_pip >= coverage)
    
    if (length(ind_candidates) == 0) {
      #warning(sprintf("No credible set found for column %d (max cumulative PIP = %.3f)", 
      #                l, max(cumsum_pip, na.rm = TRUE)))
      next
    }
    
    ind <- min(ind_candidates)
    credset <- sortinds[1:ind]
    
    # Filter by purity
    if (length(credset) == 1) {
      cred_list[[length(cred_list) + 1]] <- credset
      pip_sums <- c(pip_sums, sum(PIP[credset, l], na.rm = TRUE))
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
      pip_sums <- c(pip_sums, sum(PIP[credset, l], na.rm = TRUE))
    }
  }
  
  if (dedup) {
    cred_list_str <- sapply(cred_list, function(x) paste(x, collapse = ","))
    unique_indices <- !duplicated(cred_list_str)
    cred_list <- cred_list[unique_indices]
    pip_sums <- pip_sums[unique_indices]
  }
  
  return(list(sets = cred_list, coverage = pip_sums))
}