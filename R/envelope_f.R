.KAPPA_1_F <- 1.262
.KAPPA_2_F <- 1.536

# -----------------------------------------------------------------------------

Gamma_f <- function(s, r, x) {
  myGamma(s = s, r = r, x = x)
}

# -----------------------------------------------------------------------------

C_E_f <- function(r, x) {
  .check_f_parameters(x = x, r = r)

  rho <- r / x
  min(
    2 * .KAPPA_1_F / (pi * rho),
    2 * .KAPPA_2_F * min(1, (1 - rho) / rho)
  )
}

# -----------------------------------------------------------------------------

D_I_f <- function(s, r, x) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  .check_f_parameters(x = x, r = r)

  levy_cdf(s = s, a = x - r)
}

# -----------------------------------------------------------------------------

D_P_f <- function(s, r, x) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  .check_f_parameters(x = x, r = r)

  levy_cdf_diff(s = s, a_lower = x - r, a_upper = x + r)
}

# -----------------------------------------------------------------------------

D_E_f <- function(s, x) {
  if (!is.numeric(s)) {
    stop("s must be numeric.", call. = FALSE)
  }
  .check_positive_scalar(x, "x")

  exp(-lambda_x(x) * s)
}

# -----------------------------------------------------------------------------

.Gamma_f_u <- function(u, rho) {
  Gamma_f(s = u, r = rho, x = 1)
}

# -----------------------------------------------------------------------------

.log_Gamma_f_u <- function(u, rho) {
  log_Gamma_f(s = u, r = rho, x = 1)
}

# -----------------------------------------------------------------------------

M_I_f_u <- function(u, rho) {
  .check_positive_scalar(rho, "rho")
  if (rho >= 1) {
    stop("rho must satisfy 0 < rho < 1.", call. = FALSE)
  }

  if (!is.numeric(u)) {
    stop("u must be numeric.", call. = FALSE)
  }

  out <- rep(Inf, length(u))
  ok <- is.finite(u) & u >= uE_value() & u <= (1 + rho)^2

  if (any(ok)) {
    CE <- C_E_f(r = rho, x = 1)
    out[ok] <- (1 / rho) * levy_cdf(s = u[ok], a = 1 - rho) +
      CE * exp(-pi^2 * u[ok] / 2)
  }

  out
}

# -----------------------------------------------------------------------------

M_P_f_u <- function(u, rho) {
  .check_positive_scalar(rho, "rho")
  if (rho >= 1) {
    stop("rho must satisfy 0 < rho < 1.", call. = FALSE)
  }

  if (!is.numeric(u)) {
    stop("u must be numeric.", call. = FALSE)
  }

  out <- rep(Inf, length(u))
  upper <- (1 - rho)^2
  ok <- is.finite(u) & u >= uE_value() & u < upper

  if (any(ok)) {
    CE <- C_E_f(r = rho, x = 1)
    for (idx in which(ok)) {
      log_gamma <- .log_Gamma_f_u(u[idx], rho = rho)
      if (is.finite(log_gamma) && log_gamma < 0) {
        one_minus_gamma <- -expm1(log_gamma)
        D_pair <- levy_cdf_diff(s = u[idx], a_lower = 1 - rho, a_upper = 1 + rho)
        out[idx] <- D_pair / (rho * one_minus_gamma) +
          CE * exp(-pi^2 * u[idx] / 2)
      }
    }
  }

  out
}

# -----------------------------------------------------------------------------

.uP_interval_f <- function(rho) {
  lower <- uE_value()
  upper_domain <- (1 - rho)^2

  if (lower >= upper_domain) {
    return(NULL)
  }

  if (.log_Gamma_f_u(lower, rho = rho) >= 0) {
    return(NULL)
  }

  upper_probe <- upper_domain - max(sqrt(.Machine$double.eps) * upper_domain,
                                    .Machine$double.eps)
  if (upper_probe <= lower) {
    return(NULL)
  }

  root <- tryCatch(
    uniroot(
      f = function(u) .log_Gamma_f_u(u, rho = rho),
      lower = lower,
      upper = upper_probe,
      tol = sqrt(.Machine$double.eps)
    )$root,
    error = function(e) NA_real_
  )

  if (!is.finite(root) || root <= lower) {
    return(NULL)
  }

  c(lower, root)
}

# -----------------------------------------------------------------------------

.optimize_with_endpoints <- function(fn, lower, upper) {
  if (!is.finite(lower) || !is.finite(upper) || lower > upper) {
    stop("invalid optimization interval.", call. = FALSE)
  }

  opt <- tryCatch(
    optimize(f = fn, interval = c(lower, upper), tol = sqrt(.Machine$double.eps)),
    error = function(e) NULL
  )

  candidates <- c(lower, upper)
  if (!is.null(opt)) {
    candidates <- c(candidates, opt$minimum)
  }

  values <- vapply(candidates, function(u) fn(u), numeric(1L))
  finite <- is.finite(values)
  if (!any(finite)) {
    stop("optimization failed to find a finite envelope mass.", call. = FALSE)
  }

  candidates <- candidates[finite]
  values <- values[finite]
  idx <- which.min(values)

  list(u = candidates[idx], mass = values[idx])
}

# -----------------------------------------------------------------------------

select_envelope_f <- function(x, r) {
  .check_f_parameters(x = x, r = r)

  rho <- r / x
  uI <- .optimize_with_endpoints(
    fn = function(u) M_I_f_u(u = u, rho = rho),
    lower = uE_value(),
    upper = (1 + rho)^2
  )

  p_interval <- .uP_interval_f(rho = rho)
  uP <- NULL

  if (!is.null(p_interval)) {
    uP <- .optimize_with_endpoints(
      fn = function(u) M_P_f_u(u = u, rho = rho),
      lower = p_interval[1],
      upper = p_interval[2]
    )
  }

  if (is.null(uP) || uI$mass <= uP$mass) {
    method <- "I"
    u_star <- uI$u
    mass <- uI$mass
  } else {
    method <- "P"
    u_star <- uP$u
    mass <- uP$mass
  }

  list(
    method = method,
    u = u_star,
    s = x^2 * u_star,
    mass = mass,
    rho = rho,
    mass_I = uI$mass,
    mass_P = if (is.null(uP)) NA_real_ else uP$mass,
    u_I = uI$u,
    u_P = if (is.null(uP)) NA_real_ else uP$u,
    CE = C_E_f(r = r, x = x)
  )
}
