library(pracma)

#' Convert interval probabilities (length K-1) to grid-point weights (length K)
#'
#' We view the K-1 intervals as (Theta[k], Theta[k+1]].
#' A simple, mass-preserving projection is:
#'   w[1]   = interval[1] / 2
#'   w[K]   = interval[K-1] / 2
#'   w[k]   = (interval[k-1] + interval[k]) / 2,  k = 2,...,K-1
#'
#' This preserves total mass: sum(w) == sum(interval).
interval_to_grid_weights <- function(interval_prob) {
  K_minus_1 <- length(interval_prob)
  K <- K_minus_1 + 1L
  w <- numeric(K)

  if (K_minus_1 == 1L) {
    # Degenerate case: just split mass equally between the two endpoints
    w[1] <- interval_prob[1] / 2
    w[2] <- interval_prob[1] / 2
  } else {
    w[1] <- interval_prob[1] / 2
    w[K] <- interval_prob[K_minus_1] / 2
    if (K > 2L) {
      for (k in 2:(K - 1L)) {
        w[k] <- (interval_prob[k - 1L] + interval_prob[k]) / 2
      }
    }
  }
  # Normalize just in case of numerical noise
  w / sum(w)
}


#' Compute subinterval counts N_{l,i} from beta vector
#' 
#' @param beta Vector of coefficients (length p)
#' @param endpoints Vector of interval endpoints (length 2^L + 1)
#' @param L Number of levels in Polya tree
#' @return List with N matrix (L x max_intervals_per_level) containing counts
compute_interval_counts <- function(beta, endpoints, L) {
  # Pre-allocate matrix: row l has 2^l intervals
  max_intervals <- 2^L
  N <- matrix(0, nrow = L, ncol = max_intervals)
  
  # For each level l = 1, ..., L
  for (l in 1:L) {
    n_intervals <- 2^l
    # Compute interval boundaries for this level
    # Level l divides [a_min, a_max] into 2^l equal parts
    for (i in 1:n_intervals) {
      # Interval I_{l,i} = (endpoints[(i-1)*2^(L-l) + 1], endpoints[i*2^(L-l) + 1]]
      left_idx <- (i - 1) * 2^(L - l) + 1
      right_idx <- i * 2^(L - l) + 1
      
      left_bound <- endpoints[left_idx]
      right_bound <- endpoints[right_idx]
      
      # Count beta values in this interval (left, right]
      if (i == n_intervals) {
        N[l, i] <- sum(beta >= left_bound & beta <= right_bound)
      } else {
        N[l, i] <- sum(beta >= left_bound & beta < right_bound)
      }

    }
  }
  
  return(N)
}


#' Sample phi parameters from their conditional posterior (conjugate update)
#'
#' According to Eq. (7) in the paper:
#' phi_{l,i} | N ~ Beta(1 + N_{l, 2*i-1}, 1 + N_{l, 2*i})
#'
#' @param N Matrix of interval counts (L x max_intervals)
#' @param L Number of levels
#' @return List of phi values organized by level
sample_phi_given_beta <- function(N, L) {
  phi <- vector("list", L)
  
  for (l in 1:L) {
    n_parent_intervals <- 2^(l - 1)
    phi[[l]] <- numeric(n_parent_intervals)

    level_param = 4

    
    for (i in 1:n_parent_intervals) {
      # For parent interval i at level l-1, 
      # we have child intervals 2*i-1 and 2*i at level l
      left_child <- 2 * i - 1
      right_child <- 2 * i
      
      # Beta parameters: alpha = 1 + N_{l, left_child}, beta = 1 + N_{l, right_child}
      alpha <- level_param + N[l, left_child]
      beta_param <- level_param + N[l, right_child]
      
      # Sample from Beta distribution
      phi[[l]][i] <- rbeta(1, alpha, beta_param)
    }
  }
  
  return(phi)
}


#' Posterior means E[phi | N] for Rao–Blackwellization
#'
#' For each node (l, i) we have the posterior:
#'   phi_{l,i} | N ~ Beta(1 + N_{l,2*i-1}, 1 + N_{l,2*i})
#' so
#'   E[phi_{l,i} | N] = (1 + N_{l,2*i-1}) / (2 + N_{l,2*i-1} + N_{l,2*i})
phi_posterior_means <- function(N, L) {
  phi_mean <- vector("list", L)
  
  for (l in 1:L) {
    n_parent_intervals <- 2^(l - 1)
    phi_mean[[l]] <- numeric(n_parent_intervals)
    
    level_param = 4

    for (i in 1:n_parent_intervals) {
      left_child  <- 2 * i - 1
      right_child <- 2 * i
      
      alpha      <- level_param + N[l, left_child]
      beta_param <- level_param + N[l, right_child]
      
      phi_mean[[l]][i] <- alpha / (alpha + beta_param)
    }
  }
  
  phi_mean
}



#' Convert phi parameters to interval probabilities pi at level L
#'
#' The interval probabilities are products along the tree:
#' pi_{1,1} = phi_{1,1}, pi_{1,2} = 1 - phi_{1,1}
#' pi_{l, 2*i-1} = phi_{l,i} * pi_{l-1, i}
#' pi_{l, 2*i} = (1 - phi_{l,i}) * pi_{l-1, i}
#'
#' @param phi List of phi values by level
#' @param L Number of levels
#' @return Vector of interval probabilities at level L (length 2^L)
phi_to_interval_probs <- function(phi, L) {
  # Start with level 1
  pi_prev <- c(phi[[1]][1], 1 - phi[[1]][1])
  
  if (L == 1) {
    return(pi_prev)
  }
  
  # Iterate through levels 2 to L
  for (l in 2:L) {
    n_intervals <- 2^l
    pi_curr <- numeric(n_intervals)
    
    for (i in 1:(2^(l-1))) {
      # Children of interval i at level l-1
      left_child <- 2 * i - 1
      right_child <- 2 * i
      
      pi_curr[left_child] <- phi[[l]][i] * pi_prev[i]
      pi_curr[right_child] <- (1 - phi[[l]][i]) * pi_prev[i]
    }
    
    pi_prev <- pi_curr
  }
  
  return(pi_prev)
}


#' Find which interval a value belongs to
#' @param x A single value
#' @param endpoints Vector of interval endpoints (length K)
#' @return Index i such that x is in interval (endpoints[i], endpoints[i+1]]
find_interval_index <- function(x, endpoints) {
  # findInterval returns i such that endpoints[i] < x <= endpoints[i+1]
  idx <- findInterval(x, endpoints, rightmost.closed = TRUE)
  # Handle edge case: if x <= min(endpoints), place in first interval
  if (idx == 0) idx <- 1
  return(idx)
}



#' Compute log probability of an interval under current Polya tree
#' @param interval_idx Interval index (1 to 2^L)
#' @param phi List of phi values by level
#' @param L Number of levels
#' @return Log probability
log_interval_prob <- function(interval_idx, phi, L) {
  log_prob <- 0.0
  current_idx <- interval_idx
  
  # Small epsilon to avoid log(0)
  eps <- 1e-10
  
  # Traverse from level L up to level 1
  for (l in L:1) {
    # Determine parent interval at level l-1
    parent_idx <- ceiling(current_idx / 2)
    
    # Is current_idx the left (odd) or right (even) child?
    is_left_child <- (current_idx %% 2 == 1)
    
    if (is_left_child) {
      # P(left | parent) = phi[l][parent_idx]
      p <- phi[[l]][parent_idx]
      log_prob <- log_prob + log(p)
    } else {
      # P(right | parent) = 1 - phi[l][parent_idx]
      p <- 1 - phi[[l]][parent_idx]
      log_prob <- log_prob + log(p)
    }
    
    current_idx <- parent_idx
  }
  
  return(log_prob)
}


#' Compute log posterior weights for each interval for a single beta_i
#'
#' @param mean  Conditional Gaussian mean m for beta_i
#' @param sd    Conditional Gaussian sd for beta_i
#' @param phi   List of Polya tree phi parameters
#' @param L     Number of Polya tree levels
#' @param endpoints Grid endpoints (length K = 2^L + 1)
#' @return      Vector of log-weights of length (K - 1)
log_interval_posterior_weights <- function(mean, sd, phi, L, endpoints) {
  num_intervals <- length(endpoints) - 1L
  log_w <- rep(-Inf, num_intervals)
  
  for (r in 1:num_intervals) {
    lower <- endpoints[r]
    upper <- endpoints[r + 1L]
    
    # Likelihood mass of N(mean, sd^2) inside (lower, upper]
    a <- (lower - mean) / sd
    b <- (upper - mean) / sd
    mass <- pnorm(b) - pnorm(a)
    
    if (mass > 0) {
      log_prior_r <- log_interval_prob(r, phi, L)
      log_w[r] <- log_prior_r + log(mass)
    }
    # if mass == 0, log_w[r] stays -Inf => prob = 0
  }
  
  log_w
}


# Fast vectorized version: prior given by log_pi, likelihood mass via CDF at endpoints
fast_log_interval_posterior_weights <- function(mean, sd, log_pi, endpoints) {
  # endpoints: length K = 2^L + 1
  # log_pi:   length K-1 = 2^L  (log prior per interval)

  # Standardize endpoints once
  z    <- (endpoints - mean) / sd
  cdf  <- pnorm(z)

  # Mass in each interval (endpoints[r], endpoints[r+1]]
  mass <- diff(cdf)            # length K-1

  # Avoid log(0)
  mass[mass <= 0] <- .Machine$double.xmin

  # Combine prior and likelihood mass
  log_pi + log(mass)
}


#' Truncated normal sampling using inverse CDF method
#' @param n Number of samples
#' @param mean Mean of normal distribution
#' @param sd Standard deviation
#' @param lower Lower bound
#' @param upper Upper bound
#' @return Samples from truncated normal
rtruncnorm <- function(n = 1, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  a <- (lower - mean) / sd
  b <- (upper - mean) / sd
  
  # CDF at bounds
  Fa <- pnorm(a)
  Fb <- pnorm(b)
  
  # Sample uniform in [Fa, Fb], then apply inverse CDF
  u <- runif(n, Fa, Fb)
  z <- qnorm(u)
  
  mean + sd * z
}


#' Main Polya Tree estimation function with MH-within-Gibbs sampler
#'
#' @param X Design matrix (n x p)
#' @param y Response vector (length n)
#' @param sigma Known noise standard deviation
#' @param grid Grid points for density estimation (length K)
#' @param K Number of grid points (should equal 2^L + 1)
#' @param L Number of Polya tree levels
#' @param max.iter Maximum number of MCMC iterations
#' @param verbose Print progress
#' @param w.true True density for diagnostics (optional)
#' @param theta.true True betas for MSE diagnostics (optional)
#' @param print.iter Print every print.iter iterations
polya_tree_skeleton <- function(X, y, sigma, grid, K = 65, L = 6,
                                max.iter = 60000, 
                                verbose = TRUE, w.true = NULL, theta.true = NULL, 
                                print.iter = 100) {
  
  # Dimensions
  n <- nrow(X)
  p <- ncol(X)
  burn_in <- 200

  
  # Verify K = 2^L + 1
  expected_K <- 2^L + 1
  if (K != expected_K) {
    warning(sprintf("K=%d but 2^L + 1 = %d. Using K=%d", K, expected_K, expected_K))
    K <- expected_K
  }
  
  # Initialize with ridge solution
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  # beta <- as.numeric(solve(XtX, Xty))
  # beta <- as.numeric(solve(XtX + diag(1e-6, p), Xty))
  set initial beta to vector of 0s
  beta <- rep(0, p)


  # Polya tree endpoints
  if (length(grid) != K) {
    stop(sprintf("grid must have length K=%d", K))
  }
  endpoints <- grid
  num_intervals <- 2^L  # K - 1 intervals
  
  # Pre-compute hessians for Gaussian likelihood: -d²logL/dβ²ᵢ = ||X[:,i]||² / σ²
  sigma_sq <- sigma^2
  hessians <- -colSums(X^2) / sigma_sq
  
  # Current state
  Xbeta <- as.numeric(X %*% beta)
  residuals <- y - Xbeta  
  
  # Storage for MCMC samples
  w_samples_train <- matrix(0, nrow = max.iter, ncol = K)
  beta_samples_train <- matrix(0, nrow = max.iter, ncol = p)
  
  # Initialize phi for first iteration: start with a flat tree (like Python)
  phi <- vector("list", L)
  for (l in 1:L) {
    n_parent_intervals <- 2^(l - 1)
    phi[[l]] <- rep(0.5, n_parent_intervals)
  }
  
  
  # ============================================================================
  # MAIN MCMC LOOP
  # ============================================================================
  for (iter in 1:max.iter) {
    
        # Compute interval priors pi_r under current phi for this whole sweep
    pi_L <- phi_to_interval_probs(phi, L)       # length = num_intervals = 2^L
    eps  <- 1e-12
    log_pi <- log(pi_L + eps)

    # --- Step 1: Update each beta[i] via exact Gibbs over intervals ---
    for (i in 1:p) {
        Xi <- X[, i]

        grad <- sum(residuals * Xi) / sigma_sq
        hess <- hessians[i]
        mean_i <- beta[i] - grad / hess
        sd_i   <- 1 / sqrt(-hess)

        # Now use a faster, vectorized version that takes log_pi
        log_w <- fast_log_interval_posterior_weights(mean_i, sd_i, log_pi, endpoints)

        max_log_w <- max(log_w)
        w_int <- exp(log_w - max_log_w)
        w_int[!is.finite(w_int)] <- 0
        w_int <- w_int / sum(w_int)

        r_star <- sample.int(num_intervals, size = 1L, prob = w_int)

        lower_bound <- endpoints[r_star]
        upper_bound <- endpoints[r_star + 1L]
        beta_star   <- rtruncnorm(1, mean = mean_i, sd = sd_i,
                                  lower = lower_bound, upper = upper_bound)

        delta_beta <- beta_star - beta[i]
        if (delta_beta != 0) {
          delta_pred <- Xi * delta_beta
          Xbeta      <- Xbeta + delta_pred
          residuals  <- residuals - delta_pred
        }

        beta[i]      <- beta_star
    }

    # --- Step 2: Update Polya tree parameters phi | beta (Gibbs step) ---
    N <- compute_interval_counts(beta, endpoints, L)
    phi <- sample_phi_given_beta(N, L)

    # --- Step 3 (Rao–Blackwellized): use E[w | beta] instead of w(phi) ---
    # Posterior means of phi given N
    phi_mean <- phi_posterior_means(N, L)
    
    # Expected interval probabilities E[pi_L | N]
    pi_L_rb <- phi_to_interval_probs(phi_mean, L)
    
    # Project to grid weights: this is w_RB = E[w | beta]
    w_rb <- interval_to_grid_weights(pi_L_rb)
    
    # --- Save samples (Rao–Blackwellized w) ---
    w_samples_train[iter, ]     <- w_rb
    beta_samples_train[iter, ]  <- beta
    
    # --- Print progress ---

    if (verbose && (iter %% print.iter == 0)) {
      cat(sprintf("Iteration %d/%d\n", iter, max.iter))
      
      # Compute diagnostics

      mse_inst <- mean((beta - theta.true)^2)
      cat(sprintf("  Instantaneous MSE: %.4f\n", mse_inst))

      if (iter > burn_in) {
        w_guess <- colMeans(w_samples_train[(burn_in + 1):iter, , drop = FALSE])
        theta_guess <- colMeans(beta_samples_train[(burn_in + 1):iter, , drop = FALSE])
        
        if (!is.null(w.true)) {
          err <- sum(abs(w_guess - w.true)) / 2
          cat(sprintf("  TV distance: %.4f\n", err))
        }
        
        if (!is.null(theta.true)) {
          mse <- mean((theta_guess - theta.true)^2)
          cat(sprintf("  MSE of beta: %.4f\n", mse))
        }
      }
    }
  }
  
  # Return results
  w_guess <- colMeans(w_samples_train[(burn_in + 1):iter, , drop = FALSE])
  theta_guess <- colMeans(beta_samples_train[(burn_in + 1):iter, , drop = FALSE])



  list(
    w = w_guess,
    theta.mean = theta_guess,
    w_samples = w_samples_train,
    theta_samples = beta_samples_train,
    grid = grid,
    sigma = sigma,
    K = K,
    L = L
  )
}