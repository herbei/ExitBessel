simulate_q_I_g <- function(x, s, n = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("x must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s <= 0 || s >= x^2) {
    stop("s must be a finite scalar with 0 < s < x^2.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n) ||
      n > .Machine$integer.max) {
    stop("n must be a non-negative integer scalar.", call. = FALSE)
  }

  n <- as.integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }

  draws <- numeric(n)
  n_accepted <- 0L
  x2 <- x^2

  log_accept_prob <- log(x) + 0.5 * log(2 / (pi * s)) - x2 / (2 * s)
  accept_prob <- exp(log_accept_prob)

  while (n_accepted < n) {
    n_needed <- n - n_accepted
    batch_size <- max(1000L, ceiling(1.2 * n_needed / max(accept_prob, .Machine$double.eps)))
    batch_size <- min(batch_size, 1000000L)

    g <- rgamma(batch_size, shape = 3 / 2, rate = 1)
    t <- x2 / (2 * g)
    keep <- t <= s & runif(batch_size) <= 1 - t / x2

    if (any(keep)) {
      accepted <- t[keep]
      n_take <- min(length(accepted), n_needed)
      idx <- seq.int(n_accepted + 1L, n_accepted + n_take)
      draws[idx] <- accepted[seq_len(n_take)]
      n_accepted <- n_accepted + n_take
    }
  }

  draws
}

# -----------------------------------------------------------------------------

simulate_g <- function(x, n = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("x must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n) ||
      n > .Machine$integer.max) {
    stop("n must be a non-negative integer scalar.", call. = FALSE)
  }

  n <- as.integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }

  u_star_g <- 0.403723
  lambda <- pi^2 / (2 * x^2)
  s <- x^2 * u_star_g
  tE <- 2 * x^2 * log(2) / pi^2

  if (s < tE) {
    stop("the hard-coded switching point must satisfy s >= 2 * x^2 * log(2) / pi^2.",
         call. = FALSE)
  }

  C_I <- C_I_g(s = s, x = x)
  C_E <- C_E_g()
  D_I <- D_I_g(s = s, x = x)
  D_E <- D_E_g(s = s, x = x)
  M <- C_I * D_I + C_E * D_E
  prob_I <- C_I * D_I / M

  draws <- numeric(n)

  for (i in seq_len(n)) {
    accepted <- FALSE
    while (!accepted) {
      B <- runif(1) <= prob_I

      if (B) {
        Y <- simulate_q_I_g(x = x, s = s, n = 1L)
        A <- C_I * .b_g_I_term(t = Y, x = x, k = 0)
      } else {
        E <- rexp(1, rate = 1)
        Y <- s + E / lambda
        f_E <- lambda * exp(-lambda * Y)
        A <- C_E * f_E
      }

      threshold <- runif(1) * A
      N <- 0L
      decided <- FALSE

      while (!decided) {
        if (Y <= s) {
          lower <- LN_g_I(t = Y, x = x, N = N)
          upper <- UN_g_I(t = Y, x = x, N = N)
        } else {
          lower <- LN_g_E(t = Y, x = x, N = N)
          upper <- UN_g_E(t = Y, x = x, N = N)
        }

        if (threshold <= lower) {
          draws[i] <- Y
          accepted <- TRUE
          decided <- TRUE
        } else if (threshold > upper) {
          decided <- TRUE
        } else {
          N <- N + 1L
        }
      }
    }
  }

  draws
}
