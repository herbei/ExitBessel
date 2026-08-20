f_image_bounds <- function(t, x, r, N = 0L) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }
  .check_f_parameters(x = x, r = r)
  N <- .check_nonnegative_integer(N, "N")

  d_minus <- x - r
  d_plus <- x + r
  out_lower <- rep(NA_real_, length(t))
  out_upper <- rep(NA_real_, length(t))
  ok <- is.finite(t) & t > 0 & t <= d_plus^2

  if (any(ok)) {
    for (idx in which(ok)) {
      tt <- t[idx]
      pair_sum <- 0

      for (k in 0:N) {
        a_lower <- d_minus + 2 * k * x
        pair_sum <- pair_sum + levy_density_pair_diff(tt, a_lower = a_lower, r = r)
      }

      lower <- (x / r) * pair_sum
      upper <- lower + (x / r) * levy_density(tt, a = d_plus + 2 * N * x)

      out_lower[idx] <- lower
      out_upper[idx] <- upper
    }
  }

  list(lower = out_lower, upper = out_upper)
}

# -----------------------------------------------------------------------------

LN_f_I <- function(t, x, r, N = 0L) {
  f_image_bounds(t = t, x = x, r = r, N = N)$lower
}

# -----------------------------------------------------------------------------

UN_f_I <- function(t, x, r, N = 0L) {
  f_image_bounds(t = t, x = x, r = r, N = N)$upper
}

# -----------------------------------------------------------------------------

f_eigen_partial_sum <- function(t, x, r, N = 0L) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }
  .check_f_parameters(x = x, r = r)
  N <- .check_nonnegative_integer(N, "N")

  out <- rep(NA_real_, length(t))
  ok <- is.finite(t) & t > 0

  if (!any(ok)) {
    return(out)
  }

  if (N == 0L) {
    out[ok] <- 0
    return(out)
  }

  tt <- t[ok]
  n <- seq_len(N)
  lambda <- pi^2 / (2 * x^2)
  weights <- (-1)^(n + 1) * n * sin(n * pi * r / x)
  exp_mat <- exp(-outer(n^2, lambda * tt, "*"))
  out[ok] <- as.numeric((pi / (x * r)) * (weights %*% exp_mat))

  out
}

# -----------------------------------------------------------------------------

N0_f_E <- function(a) {
  .check_positive_scalar(a, "a")

  N <- 0L
  repeat {
    q1 <- ((N + 2) / (N + 1)) * exp(-a * (2 * N + 3))
    q2 <- ((N + 2) / (N + 1))^2 * exp(-a * (2 * N + 3))

    if (q1 < 1 && q2 < 1) {
      return(N)
    }

    N <- N + 1L
    if (N > .Machine$integer.max - 1L) {
      stop("failed to find a finite spectral starting index.", call. = FALSE)
    }
  }
}

# -----------------------------------------------------------------------------

.f_eigen_remainder_bound <- function(a, x, r, N) {
  q1 <- ((N + 2) / (N + 1)) * exp(-a * (2 * N + 3))
  q2 <- ((N + 2) / (N + 1))^2 * exp(-a * (2 * N + 3))

  B1 <- pi * (N + 1) * exp(-a * (N + 1)^2) / (x * r * (1 - q1))
  B2 <- pi^2 * (N + 1)^2 * exp(-a * (N + 1)^2) / (x^2 * (1 - q2))

  min(B1, B2)
}

# -----------------------------------------------------------------------------

f_eigen_bounds <- function(t, x, r, N = NULL) {
  .check_finite_scalar(t, "t")
  if (t <= 0) {
    stop("t must be positive.", call. = FALSE)
  }
  .check_f_parameters(x = x, r = r)

  a <- lambda_x(x) * t
  N0 <- N0_f_E(a)
  if (is.null(N)) {
    N <- N0
  } else {
    N <- .check_nonnegative_integer(N, "N")
    N <- max(N, N0)
  }

  if (N > 0L) {
    n <- seq_len(N)
    lambda <- pi^2 / (2 * x^2)
    terms <- (pi / (x * r)) * (-1)^(n + 1) * n *
      sin(n * pi * r / x) * exp(-lambda * t * n^2)
    partial <- cumsum(terms)
  } else {
    partial <- numeric(0L)
  }

  lower <- -Inf
  upper <- Inf

  for (k in seq.int(N0, N)) {
    Sk <- if (k == 0L) 0 else partial[k]
    Bk <- .f_eigen_remainder_bound(a = a, x = x, r = r, N = k)
    lower <- max(lower, Sk - Bk)
    upper <- min(upper, Sk + Bk)
  }

  list(lower = lower, upper = upper, N0 = N0, N = N)
}

# -----------------------------------------------------------------------------

LN_f_E <- function(t, x, r, N = NULL) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }

  vapply(
    t,
    function(tt) f_eigen_bounds(t = tt, x = x, r = r, N = N)$lower,
    numeric(1L)
  )
}

# -----------------------------------------------------------------------------

UN_f_E <- function(t, x, r, N = NULL) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }

  vapply(
    t,
    function(tt) f_eigen_bounds(t = tt, x = x, r = r, N = N)$upper,
    numeric(1L)
  )
}

# -----------------------------------------------------------------------------

f_density_spectral <- function(t, x, r, N = 200L) {
  f_eigen_partial_sum(t = t, x = x, r = r, N = N)
}
