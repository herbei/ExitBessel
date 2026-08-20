.check_finite_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop(name, " must be a finite numeric scalar.", call. = FALSE)
  }

  invisible(value)
}

# -----------------------------------------------------------------------------

.check_positive_scalar <- function(value, name) {
  .check_finite_scalar(value, name)
  if (value <= 0) {
    stop(name, " must be positive.", call. = FALSE)
  }

  invisible(value)
}

# -----------------------------------------------------------------------------

.check_nonnegative_integer <- function(value, name = "n") {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || value != floor(value) || value > .Machine$integer.max) {
    stop(name, " must be a non-negative integer scalar.", call. = FALSE)
  }

  as.integer(value)
}

# -----------------------------------------------------------------------------

.check_f_parameters <- function(x, r) {
  .check_positive_scalar(x, "x")
  .check_positive_scalar(r, "r")

  if (r >= x) {
    stop("r must satisfy 0 < r < x.", call. = FALSE)
  }

  invisible(TRUE)
}
