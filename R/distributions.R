lambda_x <- function(x) {
  .check_positive_scalar(x, "x")
  pi^2 / (2 * x^2)
}

# -----------------------------------------------------------------------------

tE_x <- function(x) {
  .check_positive_scalar(x, "x")
  2 * x^2 * log(2) / pi^2
}

# -----------------------------------------------------------------------------

uE_value <- function() {
  2 * log(2) / pi^2
}

# -----------------------------------------------------------------------------

levy_density <- function(t, a) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }
  .check_positive_scalar(a, "a")

  out <- rep(0, length(t))
  out[is.na(t)] <- NA_real_
  ok <- is.finite(t) & t > 0

  if (any(ok)) {
    tt <- t[ok]
    out[ok] <- a / (sqrt(2 * pi) * tt^(3 / 2)) * exp(-a^2 / (2 * tt))
  }

  out
}

# -----------------------------------------------------------------------------

levy_log_cdf <- function(s, a) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  .check_positive_scalar(a, "a")

  out <- rep(-Inf, length(s))
  out[is.na(s)] <- NA_real_
  out[is.infinite(s) & s > 0] <- 0
  ok <- is.finite(s) & s > 0

  if (any(ok)) {
    out[ok] <- log(2) + pnorm(-a / sqrt(s[ok]), log.p = TRUE)
  }

  out
}

# -----------------------------------------------------------------------------

levy_cdf <- function(s, a) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  .check_positive_scalar(a, "a")

  exp(levy_log_cdf(s = s, a = a))
}

# -----------------------------------------------------------------------------

.normal_cdf_diff <- function(lower, upper) {
  out <- pnorm(upper) - pnorm(lower)

  bad <- is.finite(lower) & is.finite(upper) & lower < upper & out <= 0
  if (any(bad)) {
    log_upper <- pnorm(upper[bad], log.p = TRUE)
    log_lower <- pnorm(lower[bad], log.p = TRUE)
    ratio <- exp(log_lower - log_upper)
    out[bad] <- exp(log_upper + log1p(-ratio))
  }

  out
}

# -----------------------------------------------------------------------------

levy_cdf_diff <- function(s, a_lower, a_upper) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  .check_positive_scalar(a_lower, "a_lower")
  .check_positive_scalar(a_upper, "a_upper")
  if (a_lower >= a_upper) {
    stop("a_lower must be smaller than a_upper.", call. = FALSE)
  }

  log_lower <- levy_log_cdf(s = s, a = a_lower)
  log_upper <- levy_log_cdf(s = s, a = a_upper)
  out <- rep(0, length(s))
  out[is.na(s)] <- NA_real_

  ok <- is.finite(log_lower) & is.finite(log_upper)
  if (any(ok)) {
    log_ratio <- pmin(log_upper[ok] - log_lower[ok], 0)
    out[ok] <- exp(log_lower[ok]) * (-expm1(log_ratio))
  }

  out
}

# -----------------------------------------------------------------------------

levy_density_pair_diff <- function(t, a_lower, r) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }
  .check_positive_scalar(a_lower, "a_lower")
  .check_positive_scalar(r, "r")

  out <- rep(0, length(t))
  out[is.na(t)] <- NA_real_
  ok <- is.finite(t) & t > 0

  if (any(ok)) {
    tt <- t[ok]
    ell_lower <- levy_density(tt, a = a_lower)
    log_ratio <- log1p(2 * r / a_lower) - 2 * r * (a_lower + r) / tt
    out[ok] <- -ell_lower * expm1(log_ratio)
  }

  out
}

# -----------------------------------------------------------------------------

exp_density_rate <- function(t, rate) {
  if (!is.numeric(t)) {
    stop("t must be numeric.", call. = FALSE)
  }
  .check_positive_scalar(rate, "rate")

  out <- rep(0, length(t))
  out[is.na(t)] <- NA_real_
  ok <- is.finite(t) & t >= 0

  if (any(ok)) {
    out[ok] <- rate * exp(-rate * t[ok])
  }

  out
}

# -----------------------------------------------------------------------------

rinv_gamma_truncated_upper <- function(n, shape, beta, upper) {
  n <- .check_nonnegative_integer(n)
  .check_positive_scalar(shape, "shape")
  .check_positive_scalar(beta, "beta")
  .check_positive_scalar(upper, "upper")

  if (n == 0L) {
    return(numeric(0L))
  }

  p0 <- pgamma(beta / upper, shape = shape, rate = 1, lower.tail = TRUE)
  if (p0 >= 1) {
    stop("upper is too small for the truncated inverse-gamma sampler.",
         call. = FALSE)
  }

  u <- runif(n, min = p0, max = 1)
  g <- qgamma(u, shape = shape, rate = 1, lower.tail = TRUE)

  beta / g
}
