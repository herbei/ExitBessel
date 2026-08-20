log_Gamma_f <- function(s, r, x) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r < 0) {
    stop("r must be a finite non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= r) {
    stop("x must be a finite scalar with x > r.", call. = FALSE)
  }
  d_minus <- x - r
  if (any(!is.finite(s) | s <= 0 | s >= d_minus^2)) {
    stop("s must contain only finite values with 0 < s < (x - r)^2.",
         call. = FALSE)
  }

  log((3 * x - r)^2 - s) - log(d_minus^2 - s) -
    2 * x * (2 * x - r) / s
}

# -----------------------------------------------------------------------------

myGamma <- function(s, r, x) {
  exp(log_Gamma_f(s = s, r = r, x = x))
}

# -----------------------------------------------------------------------------

one_minus_Gamma_f <- function(s, r, x) {
  -expm1(log_Gamma_f(s = s, r = r, x = x))
}

# -----------------------------------------------------------------------------

C_I_g <- function(s, x) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("x must be a finite positive scalar.", call. = FALSE)
  }
  if (any(!is.finite(s) | s <= 0 | s >= x^2)) {
    stop("s must contain only finite values with 0 < s < x^2.", call. = FALSE)
  }

  one_minus_gamma_g <- one_minus_Gamma_f(s = s, r = 0, x = x)
  if (any(one_minus_gamma_g <= 0)) {
    stop("C_I_g requires myGamma(s, 0, x) < 1.", call. = FALSE)
  }

  1 / one_minus_gamma_g
}

# -----------------------------------------------------------------------------

C_E_g <- function(s = NULL) {
  if (is.null(s)) {
    return(2)
  }
  if (!is.numeric(s)) {
    stop("s must be numeric when supplied.", call. = FALSE)
  }

  rep(2, length(s))
}

# -----------------------------------------------------------------------------

D_I_g <- function(s, x) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("x must be a finite positive scalar.", call. = FALSE)
  }
  if (any(!is.finite(s) | s <= 0 | s >= x^2)) {
    stop("s must contain only finite values with 0 < s < x^2.", call. = FALSE)
  }

  2 * x * sqrt(2 / (pi * s)) * exp(-x^2 / (2 * s))
}

# -----------------------------------------------------------------------------

D_E_g <- function(s, x) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("x must be a finite positive scalar.", call. = FALSE)
  }
  if (any(!is.finite(s) | s < 0)) {
    stop("s must contain only finite non-negative values.", call. = FALSE)
  }

  lambda <- pi^2 / (2 * x^2)
  exp(-lambda * s)
}

# -----------------------------------------------------------------------------

.b_g_I_term <- function(t, x, k) {
  a <- (2 * k + 1)^2 * x^2
  scaled <- a / t
  2 * x / (sqrt(2 * pi) * t^(3 / 2)) * (scaled - 1) * exp(-scaled / 2)
}

# -----------------------------------------------------------------------------

g_I_partial_sum <- function(t, x, N = 0L) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("x must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 0 || N != floor(N) ||
      N > .Machine$integer.max) {
    stop("N must be a non-negative integer scalar.", call. = FALSE)
  }

  N <- as.integer(N)
  out <- rep(NA_real_, length(t))
  ok <- is.finite(t) & t > 0 & t < x^2

  if (any(ok)) {
    tt <- t[ok]
    value <- numeric(length(tt))
    for (k in 0:N) {
      value <- value + .b_g_I_term(tt, x, k)
    }
    out[ok] <- value
  }

  out
}

# -----------------------------------------------------------------------------

LN_g_I <- function(t, x, N = 0L) {
  g_I_partial_sum(t = t, x = x, N = N)
}

# -----------------------------------------------------------------------------

UN_g_I <- function(t, x, N = 0L) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("x must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 0 || N != floor(N) ||
      N > .Machine$integer.max) {
    stop("N must be a non-negative integer scalar.", call. = FALSE)
  }

  N <- as.integer(N)
  out <- rep(NA_real_, length(t))
  ok <- is.finite(t) & t > 0 & t < x^2

  if (any(ok)) {
    tt <- t[ok]
    one_minus_gamma_g <- one_minus_Gamma_f(tt, r = 0, x = x)
    bound_ok <- one_minus_gamma_g > 0

    if (any(bound_ok)) {
      tt_bound <- tt[bound_ok]
      value <- LN_g_I(tt_bound, x = x, N = N) +
        .b_g_I_term(tt_bound, x, N + 1) / one_minus_gamma_g[bound_ok]
      out[which(ok)[bound_ok]] <- value
    }
  }

  out
}

# -----------------------------------------------------------------------------

.g_E_partial_sum <- function(t, x, n_max) {
  lambda <- pi^2 / (2 * x^2)
  value <- numeric(length(t))

  for (n in seq_len(n_max)) {
    value <- value + (-1)^(n + 1) * (pi^2 / x^2) * n^2 * exp(-lambda * t * n^2)
  }

  value
}

# -----------------------------------------------------------------------------

LN_g_E <- function(t, x, N = 0L) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("x must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 0 || N != floor(N) ||
      N > .Machine$integer.max) {
    stop("N must be a non-negative integer scalar.", call. = FALSE)
  }

  N <- as.integer(N)
  tE <- 2 * x^2 * log(2) / pi^2
  out <- rep(NA_real_, length(t))
  ok <- is.finite(t) & t >= tE

  if (any(ok)) {
    out[ok] <- .g_E_partial_sum(t[ok], x, 2 * N + 2)
  }

  out
}

# -----------------------------------------------------------------------------

UN_g_E <- function(t, x, N = 0L) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("x must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 0 || N != floor(N) ||
      N > .Machine$integer.max) {
    stop("N must be a non-negative integer scalar.", call. = FALSE)
  }

  N <- as.integer(N)
  tE <- 2 * x^2 * log(2) / pi^2
  out <- rep(NA_real_, length(t))
  ok <- is.finite(t) & t >= tE

  if (any(ok)) {
    out[ok] <- .g_E_partial_sum(t[ok], x, 2 * N + 1)
  }

  out
}
