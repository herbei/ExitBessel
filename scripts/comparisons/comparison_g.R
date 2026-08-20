set.seed(20260818)

.script_file <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
  } else {
    normalizePath(sys.frames()[[1L]]$ofile, mustWork = TRUE)
  }
})
source(file.path(dirname(dirname(.script_file)), "script_utils.R"))
.project_root <- exit_bessel_project_root(.script_file)

# Source SH26 code into its own environment.
sh26_env <- new.env(parent = globalenv())
source(
  file.path(.project_root, "R", "SH26_functions.R"),
  local = sh26_env
)

# Source the current code into its own environment.
current_env <- new.env(parent = globalenv())
source_exit_bessel(root = .project_root, local = current_env)

# Number of accepted draws to generate for each sampler and each x.
B <- 10000L

# Values of x used in the comparison.
x_values <- c(0.5, 1, 2, 4, 8)

# -----------------------------------------------------------------------------

sh26_expected_proposals <- function(x) {
  # Constants used by SH26 sample_from_g().
  lam <- pi^2 / (2 * x^2)
  c_g <- log(3) * log(4) / (3 * pi^2)
  tstar <- c_g * 2 * x^2 / log(3)
  gamma_star <- 9 * exp(-4 * x^2 / tstar)

  # First proposal component mass.
  D1 <- 2 / (1 - gamma_star)
  N1 <- sh26_env$F_b(tstar, 3 / 2, x^2 / 2)
  tildeD1 <- D1 * N1

  # Second proposal component mass.
  D2 <- 2
  N2 <- exp(-lam * tstar)
  tildeD2 <- D2 * N2

  # The proposal envelope integrates to tildeD1 + tildeD2.
  tildeD1 + tildeD2
}

# -----------------------------------------------------------------------------

current_expected_proposals <- function(x) {
  # Constants used by the current Algorithm 5 implementation.
  u_star_g <- 0.403723
  s <- x^2 * u_star_g

  # The current proposal envelope integrates to C_I D_I + C_E D_E.
  current_env$C_I_g(s = s, x = x) * current_env$D_I_g(s = s, x = x) +
    current_env$C_E_g() * current_env$D_E_g(s = s, x = x)
}

# -----------------------------------------------------------------------------

sh26_sample_g_counted <- function(x) {
  # Constants copied from SH26 sample_from_g().
  lam <- pi^2 / (2 * x^2)
  c_g <- log(3) * log(4) / (3 * pi^2)
  tstar <- c_g * 2 * x^2 / log(3)
  gamma_star <- 9 * exp(-4 * x^2 / tstar)

  D1 <- 2 / (1 - gamma_star)
  N1 <- sh26_env$F_b(tstar, 3 / 2, x^2 / 2)
  tildeD1 <- D1 * N1

  D2 <- 2
  N2 <- exp(-lam * tstar)
  tildeD2 <- D2 * N2

  total_mass <- tildeD1 + tildeD2
  prob_ig <- tildeD1 / total_mass

  proposal_count <- 0L
  bound_iterations <- 0L
  series_terms <- 0L

  repeat {
    # Count each proposed Y.
    proposal_count <- proposal_count + 1L

    # Draw from the SH26 proposal mixture.
    if (runif(1) < prob_ig) {
      Y <- sh26_env$rinv_gamma_trunc(1, 3 / 2, x^2 / 2, tstar)
    } else {
      Y <- sh26_env$rexp_trunc(1, lam, tstar, Inf)
    }

    # Evaluate the SH26 proposal density at Y.
    dq_g <- prob_ig * sh26_env$dinv_gamma_trunc_g(Y, tstar, x) +
      (1 - prob_ig) * sh26_env$dexp_trunc_g(Y, tstar, lam)

    # Draw the rejection-test uniform variable.
    U <- runif(1)
    threshold <- U * total_mass * dq_g

    # SH26 starts at N = 200 and increases by 50.
    N <- 200L

    repeat {
      N <- N + 50L
      bound_iterations <- bound_iterations + 1L

      # Count the number of eigen-series terms evaluated by upper/lower bounds.
      m_tilde <- ceiling(1 / sqrt(pi^2 * Y / (2 * x^2)))
      upper_terms <- max(2 * m_tilde - 1, 2 * N - 1)
      lower_terms <- max(2 * m_tilde, 2 * N)
      series_terms <- series_terms + upper_terms + lower_terms

      # Compute SH26 retrospective upper and lower bounds.
      UB <- sh26_env$upper_bound_g(Y, x, N)
      LB <- sh26_env$lower_bound_g(Y, x, N)

      # Accept when the lower bound is already above the threshold.
      if (threshold <= LB) {
        return(list(
          value = Y,
          proposals = proposal_count,
          qI_attempts = 0L,
          bound_iterations = bound_iterations,
          series_terms = series_terms,
          accept_N = N
        ))
      }

      # Reject this proposal when the upper bound is already below the threshold.
      if (threshold > UB) {
        break
      }
    }
  }
}

# -----------------------------------------------------------------------------

current_q_I_g_one_counted <- function(x, s) {
  # This is Algorithm 5 for one draw, with a counter for Gamma proposals.
  attempts <- 0L
  x2 <- x^2

  repeat {
    attempts <- attempts + 1L
    G <- rgamma(1, shape = 3 / 2, rate = 1)
    T <- x2 / (2 * G)
    U <- runif(1)

    if (T <= s && U <= 1 - T / x2) {
      return(list(value = T, attempts = attempts))
    }
  }
}

# -----------------------------------------------------------------------------

current_sample_g_counted <- function(x) {
  # Constants copied from the current simulate_g().
  u_star_g <- 0.403723
  lambda <- pi^2 / (2 * x^2)
  s <- x^2 * u_star_g

  C_I <- current_env$C_I_g(s = s, x = x)
  C_E <- current_env$C_E_g()
  D_I <- current_env$D_I_g(s = s, x = x)
  D_E <- current_env$D_E_g(s = s, x = x)
  total_mass <- C_I * D_I + C_E * D_E
  prob_I <- C_I * D_I / total_mass

  proposal_count <- 0L
  qI_attempts <- 0L
  bound_iterations <- 0L
  series_terms <- 0L

  repeat {
    # Count each proposed Y from the Algorithm 6 proposal mixture.
    proposal_count <- proposal_count + 1L

    # Draw from the image proposal component.
    if (runif(1) <= prob_I) {
      qI_draw <- current_q_I_g_one_counted(x = x, s = s)
      Y <- qI_draw$value
      qI_attempts <- qI_attempts + qI_draw$attempts
      A <- C_I * get(".b_g_I_term", envir = current_env)(t = Y, x = x, k = 0)
    } else {
      # Draw from the shifted exponential proposal component.
      E <- rexp(1, rate = 1)
      Y <- s + E / lambda
      f_E <- lambda * exp(-lambda * Y)
      A <- C_E * f_E
    }

    # Draw the rejection-test uniform variable.
    U <- runif(1)
    threshold <- U * A

    # The current retrospective comparison starts at N = 0.
    N <- 0L

    repeat {
      bound_iterations <- bound_iterations + 1L

      # Use image bounds to the left of s and eigen bounds to the right of s.
      if (Y <= s) {
        lower <- current_env$LN_g_I(t = Y, x = x, N = N)
        upper <- current_env$UN_g_I(t = Y, x = x, N = N)
        series_terms <- series_terms + (2 * N + 3)
      } else {
        lower <- current_env$LN_g_E(t = Y, x = x, N = N)
        upper <- current_env$UN_g_E(t = Y, x = x, N = N)
        series_terms <- series_terms + (4 * N + 3)
      }

      # Accept when the lower bound is already above the threshold.
      if (threshold <= lower) {
        return(list(
          value = Y,
          proposals = proposal_count,
          qI_attempts = qI_attempts,
          bound_iterations = bound_iterations,
          series_terms = series_terms,
          accept_N = N
        ))
      }

      # Reject this proposal when the upper bound is already below the threshold.
      if (threshold > upper) {
        break
      }

      # Otherwise, refine the approximation.
      N <- N + 1L
    }
  }
}

# -----------------------------------------------------------------------------

benchmark_sampler <- function(sampler_name, sample_fun, expected_fun, x, B) {
  proposals <- 0L
  qI_attempts <- 0L
  bound_iterations <- 0L
  series_terms <- 0L
  accept_N <- numeric(B)

  start_time <- proc.time()[["elapsed"]]

  for (i in seq_len(B)) {
    draw <- sample_fun(x)
    proposals <- proposals + draw$proposals
    qI_attempts <- qI_attempts + draw$qI_attempts
    bound_iterations <- bound_iterations + draw$bound_iterations
    series_terms <- series_terms + draw$series_terms
    accept_N[i] <- draw$accept_N
  }

  elapsed <- proc.time()[["elapsed"]] - start_time

  data.frame(
    sampler = sampler_name,
    x = x,
    accepted_draws = B,
    expected_proposals_per_draw = expected_fun(x),
    empirical_proposals_per_draw = proposals / B,
    avg_qI_attempts_per_draw = qI_attempts / B,
    runtime_seconds = elapsed,
    runtime_per_10000_draws = elapsed * 10000 / B,
    avg_bound_iterations_per_draw = bound_iterations / B,
    avg_series_terms_per_draw = series_terms / B,
    avg_accept_N = mean(accept_N)
  )
}

# -----------------------------------------------------------------------------

all_results <- list()
row_id <- 1L

for (x in x_values) {
  cat("Running SH26 sampler for x =", x, "\n")
  all_results[[row_id]] <- benchmark_sampler(
    sampler_name = "SH26",
    sample_fun = sh26_sample_g_counted,
    expected_fun = sh26_expected_proposals,
    x = x,
    B = B
  )
  row_id <- row_id + 1L

  cat("Running current sampler for x =", x, "\n")
  all_results[[row_id]] <- benchmark_sampler(
    sampler_name = "Current",
    sample_fun = current_sample_g_counted,
    expected_fun = current_expected_proposals,
    x = x,
    B = B
  )
  row_id <- row_id + 1L
}

results <- do.call(rbind, all_results)

print(results)

results_file <- file.path(.project_root, "results", "comparison_g_results.csv")
dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
write.csv(results, file = results_file, row.names = FALSE)
cat("Results written to:", normalizePath(results_file), "\n")
