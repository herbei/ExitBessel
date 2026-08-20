simulate_q_I_f <- function(x, r, s, n = 1L) {
  .check_f_parameters(x = x, r = r)
  .check_positive_scalar(s, "s")
  n <- .check_nonnegative_integer(n)

  if (n == 0L) {
    return(numeric(0L))
  }

  d_minus <- x - r
  rinv_gamma_truncated_upper(
    n = n,
    shape = 1 / 2,
    beta = d_minus^2 / 2,
    upper = s
  )
}

# -----------------------------------------------------------------------------

.simulate_q_P_f_one <- function(x, r, s) {
  d_minus <- x - r
  d_plus <- x + r
  log_ratio <- log1p(2 * r / d_minus)

  repeat {
    A <- d_minus * exp(log_ratio * runif(1))
    G <- rgamma(1, shape = 3 / 2, rate = 1)
    T <- A^2 / (2 * G)

    if (T <= s && runif(1) <= 1 - T / A^2) {
      return(T)
    }
  }
}

# -----------------------------------------------------------------------------

simulate_q_P_f <- function(x, r, s, n = 1L) {
  .check_f_parameters(x = x, r = r)
  .check_positive_scalar(s, "s")
  n <- .check_nonnegative_integer(n)

  d_minus <- x - r
  d_plus <- x + r

  if (s >= d_minus^2) {
    stop("s must satisfy 0 < s < (x - r)^2 for the first-pair sampler.",
         call. = FALSE)
  }

  if (n == 0L) {
    return(numeric(0L))
  }
  if (n == 1L) {
    return(.simulate_q_P_f_one(x = x, r = r, s = s))
  }

  draws <- numeric(n)
  n_accepted <- 0L
  accept_prob <- D_P_f(s = s, r = r, x = x) / log1p(2 * r / d_minus)
  if (!is.finite(accept_prob) || accept_prob <= 0) {
    accept_prob <- .Machine$double.eps
  } else {
    accept_prob <- min(accept_prob, 1)
  }
  log_ratio <- log1p(2 * r / d_minus)

  while (n_accepted < n) {
    n_needed <- n - n_accepted
    batch_size <- max(1000L, ceiling(1.2 * n_needed / accept_prob))
    batch_size <- min(batch_size, 1000000L)

    A <- d_minus * exp(log_ratio * runif(batch_size))
    G <- rgamma(batch_size, shape = 3 / 2, rate = 1)
    T <- A^2 / (2 * G)
    keep <- T <= s & runif(batch_size) <= 1 - T / A^2

    if (any(keep)) {
      accepted <- T[keep]
      n_take <- min(length(accepted), n_needed)
      idx <- seq.int(n_accepted + 1L, n_accepted + n_take)
      draws[idx] <- accepted[seq_len(n_take)]
      n_accepted <- n_accepted + n_take
    }
  }

  draws
}

# -----------------------------------------------------------------------------

.retrospective_accept_f <- function(threshold, y, s, x, r, max_refinements) {
  if (y <= s) {
    N <- 0L
    refinements <- 0L

    repeat {
      bounds <- f_image_bounds(t = y, x = x, r = r, N = N)
      lower <- bounds$lower
      upper <- bounds$upper

      if (threshold <= lower) {
        return(TRUE)
      }
      if (threshold > upper) {
        return(FALSE)
      }

      N <- N + 1L
      refinements <- refinements + 1L
      if (refinements > max_refinements) {
        stop("image retrospective comparison did not resolve.", call. = FALSE)
      }
    }
  }

  N <- N0_f_E(lambda_x(x) * y)
  refinements <- 0L

  repeat {
    bounds <- f_eigen_bounds(t = y, x = x, r = r, N = N)
    lower <- bounds$lower
    upper <- bounds$upper

    if (threshold <= lower) {
      return(TRUE)
    }
    if (threshold > upper) {
      return(FALSE)
    }

    N <- N + 1L
    refinements <- refinements + 1L
    if (refinements > max_refinements) {
      stop("spectral retrospective comparison did not resolve.", call. = FALSE)
    }
  }
}

# -----------------------------------------------------------------------------

simulate_f <- function(x, r, n = 1L, max_refinements = 100000L) {
  .check_f_parameters(x = x, r = r)
  n <- .check_nonnegative_integer(n)
  max_refinements <- .check_nonnegative_integer(max_refinements, "max_refinements")

  if (n == 0L) {
    return(numeric(0L))
  }

  envelope <- select_envelope_f(x = x, r = r)
  method <- envelope$method
  s <- envelope$s
  d_minus <- x - r
  lambda <- lambda_x(x)
  CE <- envelope$CE
  DE <- D_E_f(s = s, x = x)

  if (method == "I") {
    C_star <- x / r
    D_star <- D_I_f(s = s, r = r, x = x)
  } else {
    one_minus_gamma <- one_minus_Gamma_f(s = s, r = r, x = x)
    C_star <- x / (r * one_minus_gamma)
    D_star <- D_P_f(s = s, r = r, x = x)
  }

  M <- C_star * D_star + CE * DE
  prob_left <- C_star * D_star / M
  if (!is.finite(M) || M <= 0 || !is.finite(prob_left) ||
      prob_left < 0 || prob_left > 1) {
    stop("failed to construct a valid Algorithm 4 proposal envelope.",
         call. = FALSE)
  }

  draws <- numeric(n)

  for (i in seq_len(n)) {
    accepted <- FALSE

    while (!accepted) {
      if (runif(1) <= prob_left) {
        if (method == "I") {
          Y <- simulate_q_I_f(x = x, r = r, s = s, n = 1L)
          H <- C_star * levy_density(Y, a = d_minus)
        } else {
          Y <- simulate_q_P_f(x = x, r = r, s = s, n = 1L)
          H <- C_star * levy_density_pair_diff(Y, a_lower = d_minus, r = r)
        }
      } else {
        Y <- s + rexp(1, rate = 1) / lambda
        H <- CE * exp_density_rate(Y, rate = lambda)
      }

      threshold <- runif(1) * H
      accepted <- .retrospective_accept_f(
        threshold = threshold,
        y = Y,
        s = s,
        x = x,
        r = r,
        max_refinements = max_refinements
      )

      if (accepted) {
        draws[i] <- Y
      }
    }
  }

  draws
}

# -----------------------------------------------------------------------------

simulate_exit_bessel <- function(x, r = 0, n = 1L, max_refinements = 100000L) {
  .check_positive_scalar(x, "x")
  .check_finite_scalar(r, "r")
  if (r < 0 || r >= x) {
    stop("r must satisfy 0 <= r < x.", call. = FALSE)
  }

  if (r == 0) {
    return(simulate_g(x = x, n = n))
  }

  simulate_f(x = x, r = r, n = n, max_refinements = max_refinements)
}
