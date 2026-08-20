set.seed(20260818)

# Instrumented, setup-amortized benchmark for the f samplers.
# Deterministic envelope constants are computed once per (x, rho) before timing;
# the timed loop then reproduces each sampler's proposal and retrospective tests.

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

# Number of accepted draws to generate for each sampler and each rho.
B_env <- Sys.getenv("EXIT_BESSEL_B", "10000")
B <- suppressWarnings(as.integer(B_env))
if (!is.finite(B) || B <= 0L) {
  stop("EXIT_BESSEL_B must be a positive integer.", call. = FALSE)
}

# Scale is fixed because the optimized f envelope depends on rho = r / x.
x_value <- 1
rho_values <- c(0.01, 0.1, 0.4, 0.8, 0.99)

# -----------------------------------------------------------------------------

expected_mean_f <- function(x, r) {
  (x^2 - r^2) / 3
}

# -----------------------------------------------------------------------------

sh26_f_constants <- function(x, r) {
  y <- (x - r) / (2 * x)
  t1_f <- 2 * x^2 * (2 * y + 1) / log(1 + 1 / y)
  t1_f <- (log(2) / pi^2) * t1_f
  lambda <- pi^2 / (2 * x^2)
  t2_f <- log(2) / lambda
  option <- if (t1_f < t2_f) 1L else 2L

  if (option == 2L) {
    t_f <- min(t1_f, t2_f)
    gamma1_star <- sh26_env$gamma_f(t1_f, r, x)
    C1_f <- (x / r) / (1 - gamma1_star)
    N1_f <- sh26_env$F_b(t_f, 1 / 2, 2 * x^2 * y^2)
    tildeC1_f <- C1_f * N1_f

    C2_f <- pi^3 / (6 * x * r * lambda)
    N2_f <- exp(-lambda * t_f)
    tildeC2_f <- C2_f * N2_f

    total_mass <- tildeC1_f + tildeC2_f

    return(list(
      option = "option2",
      y = y,
      t1_f = t1_f,
      t2_f = t2_f,
      t_f = t_f,
      lambda = lambda,
      tildeC1_f = tildeC1_f,
      tildeC2_f = tildeC2_f,
      tildeC_f = 0,
      total_mass = total_mass,
      p1 = tildeC1_f / total_mass,
      p_mid = 0,
      p2 = tildeC2_f / total_mass
    ))
  }

  gamma1_star <- sh26_env$gamma_f(t1_f, r, x)
  C1_f <- (x / r) / (1 - gamma1_star)
  N1_f <- sh26_env$F_b(t1_f, 1 / 2, 2 * x^2 * y^2)
  tildeC1_f <- C1_f * N1_f

  a <- pi^2 * t1_f / (2 * x^2)
  Nf <- ceiling(1 / sqrt(2 * a))
  term1 <- sum(seq_len(Nf) * exp(-a * seq_len(Nf)^2))
  term2 <- (1 / (2 * a)) * exp(-a * Nf^2)
  C_f <- (pi / (x * r)) * (term1 + term2)
  tildeC_f <- C_f * (t2_f - t1_f)

  C2_f <- pi^3 / (6 * x * r * lambda)
  N2_f <- exp(-lambda * t2_f)
  tildeC2_f <- C2_f * N2_f

  total_mass <- tildeC1_f + tildeC_f + tildeC2_f

  list(
    option = "option1",
    y = y,
    t1_f = t1_f,
    t2_f = t2_f,
    t_f = NA_real_,
    lambda = lambda,
    tildeC1_f = tildeC1_f,
    tildeC2_f = tildeC2_f,
    tildeC_f = tildeC_f,
    total_mass = total_mass,
    p1 = tildeC1_f / total_mass,
    p_mid = tildeC_f / total_mass,
    p2 = tildeC2_f / total_mass
  )
}

# -----------------------------------------------------------------------------

sh26_expected_proposals <- function(x, r) {
  sh26_f_constants(x = x, r = r)$total_mass
}

# -----------------------------------------------------------------------------

make_sh26_sample_f_counted <- function(x, r) {
  constants <- sh26_f_constants(x = x, r = r)
  force(x)
  force(r)
  force(constants)

  function(max_refinements = 100000L) {
    proposal_count <- 0L
    bound_iterations <- 0L
    series_terms <- 0L

    repeat {
      proposal_count <- proposal_count + 1L

      coin <- runif(1)
      if (constants$option == "option2") {
        if (coin < constants$p1) {
          Y <- sh26_env$rinv_gamma_trunc(
            1, 1 / 2, 2 * x^2 * constants$y^2, constants$t_f
          )
        } else {
          Y <- sh26_env$rexp_trunc(
            1, constants$lambda, constants$t_f, Inf
          )
        }

        proposal_density <- constants$p1 * sh26_env$dinv_gamma_trunc_f(
          Y, constants$t_f, x, r
        ) + constants$p2 * sh26_env$dexp_trunc_f(
          Y, constants$t_f, constants$lambda
        )
      } else {
        if (coin < constants$p1) {
          Y <- sh26_env$rinv_gamma_trunc(
            1, 1 / 2, 2 * x^2 * constants$y^2, constants$t1_f
          )
        } else if (coin < constants$p1 + constants$p_mid) {
          Y <- runif(1, min = constants$t1_f, max = constants$t2_f)
        } else {
          Y <- sh26_env$rexp_trunc(
            1, constants$lambda, constants$t2_f, Inf
          )
        }

        proposal_density <- constants$p1 * sh26_env$dinv_gamma_trunc_f(
          Y, constants$t1_f, x, r
        ) + constants$p_mid * dunif(
          Y, min = constants$t1_f, max = constants$t2_f
        ) + constants$p2 * sh26_env$dexp_trunc_f(
          Y, constants$t2_f, constants$lambda
        )
      }

      threshold <- runif(1) * constants$total_mass * proposal_density
      N <- 200L
      refinements <- 0L

      repeat {
        N <- N + 50L
        bound_iterations <- bound_iterations + 1L
        refinements <- refinements + 1L

        a <- pi^2 * Y / (2 * x^2)
        m0 <- ceiling(1 / sqrt(2 * a))
        N_max <- max(200L, m0 + N)
        series_terms <- series_terms + 2L * N_max

        lower <- sh26_env$lower_bound_f(Y, x, r, N)
        upper <- sh26_env$upper_bound_f(Y, x, r, N)

        if (threshold <= lower) {
          return(list(
            value = Y,
            proposals = proposal_count,
            qP_attempts = 0L,
            bound_iterations = bound_iterations,
            series_terms = series_terms,
            accept_N = N,
            branch = constants$option
          ))
        }

        if (threshold > upper) {
          break
        }

        if (refinements > max_refinements) {
          stop("SH26 retrospective comparison did not resolve.",
               call. = FALSE)
        }
      }
    }
  }
}

# -----------------------------------------------------------------------------

current_q_P_f_one_counted <- function(x, r, s) {
  d_minus <- x - r
  d_plus <- x + r
  log_ratio <- log1p(2 * r / d_minus)
  attempts <- 0L

  repeat {
    attempts <- attempts + 1L
    A <- d_minus * exp(log_ratio * runif(1))
    G <- rgamma(1, shape = 3 / 2, rate = 1)
    T <- A^2 / (2 * G)

    if (T <= s && runif(1) <= 1 - T / A^2) {
      return(list(value = T, attempts = attempts))
    }
  }
}

# -----------------------------------------------------------------------------

current_f_constants <- function(x, r) {
  envelope <- current_env$select_envelope_f(x = x, r = r)
  method <- envelope$method
  s <- envelope$s
  d_minus <- x - r
  lambda <- current_env$lambda_x(x)
  CE <- envelope$CE
  DE <- current_env$D_E_f(s = s, x = x)

  if (method == "I") {
    C_star <- x / r
    D_star <- current_env$D_I_f(s = s, r = r, x = x)
  } else {
    one_minus_gamma <- current_env$one_minus_Gamma_f(s = s, r = r, x = x)
    C_star <- x / (r * one_minus_gamma)
    D_star <- current_env$D_P_f(s = s, r = r, x = x)
  }

  total_mass <- C_star * D_star + CE * DE
  prob_left <- C_star * D_star / total_mass

  list(
    method = method,
    s = s,
    d_minus = d_minus,
    lambda = lambda,
    CE = CE,
    C_star = C_star,
    total_mass = total_mass,
    prob_left = prob_left
  )
}

# -----------------------------------------------------------------------------

current_expected_proposals <- function(x, r) {
  current_f_constants(x = x, r = r)$total_mass
}

# -----------------------------------------------------------------------------

make_current_sample_f_counted <- function(x, r) {
  constants <- current_f_constants(x = x, r = r)
  force(x)
  force(r)
  force(constants)

  function(max_refinements = 100000L) {
    proposal_count <- 0L
    qP_attempts <- 0L
    bound_iterations <- 0L
    series_terms <- 0L

    repeat {
      proposal_count <- proposal_count + 1L

      if (runif(1) <= constants$prob_left) {
        if (constants$method == "I") {
          Y <- current_env$simulate_q_I_f(
            x = x, r = r, s = constants$s, n = 1L
          )
          envelope_density <- constants$C_star *
            current_env$levy_density(Y, a = constants$d_minus)
        } else {
          qP_draw <- current_q_P_f_one_counted(x = x, r = r, s = constants$s)
          Y <- qP_draw$value
          qP_attempts <- qP_attempts + qP_draw$attempts
          envelope_density <- constants$C_star *
            current_env$levy_density_pair_diff(
              Y, a_lower = constants$d_minus, r = r
            )
        }
      } else {
        Y <- constants$s + rexp(1, rate = 1) / constants$lambda
        envelope_density <- constants$CE *
          current_env$exp_density_rate(Y, rate = constants$lambda)
      }

      threshold <- runif(1) * envelope_density

      if (Y <= constants$s) {
        N <- 0L
        refinements <- 0L

        repeat {
          bounds <- current_env$f_image_bounds(t = Y, x = x, r = r, N = N)
          bound_iterations <- bound_iterations + 1L
          series_terms <- series_terms + 2L * (N + 1L) + 1L

          if (threshold <= bounds$lower) {
            return(list(
              value = Y,
              proposals = proposal_count,
              qP_attempts = qP_attempts,
              bound_iterations = bound_iterations,
              series_terms = series_terms,
              accept_N = N,
              branch = constants$method
            ))
          }

          if (threshold > bounds$upper) {
            break
          }

          N <- N + 1L
          refinements <- refinements + 1L
          if (refinements > max_refinements) {
            stop("image retrospective comparison did not resolve.",
                 call. = FALSE)
          }
        }
      } else {
        N <- current_env$N0_f_E(current_env$lambda_x(x) * Y)
        refinements <- 0L

        repeat {
          bounds <- current_env$f_eigen_bounds(t = Y, x = x, r = r, N = N)
          bound_iterations <- bound_iterations + 1L
          series_terms <- series_terms + max(bounds$N, 1L)

          if (threshold <= bounds$lower) {
            return(list(
              value = Y,
              proposals = proposal_count,
              qP_attempts = qP_attempts,
              bound_iterations = bound_iterations,
              series_terms = series_terms,
              accept_N = bounds$N,
              branch = constants$method
            ))
          }

          if (threshold > bounds$upper) {
            break
          }

          N <- bounds$N + 1L
          refinements <- refinements + 1L
          if (refinements > max_refinements) {
            stop("spectral retrospective comparison did not resolve.",
                 call. = FALSE)
          }
        }
      }
    }
  }
}

# -----------------------------------------------------------------------------

benchmark_sampler <- function(sampler_name, sample_factory, expected_fun, x, rho, B) {
  r <- rho * x
  sample_one <- sample_factory(x = x, r = r)
  values <- numeric(B)
  proposals <- 0L
  qP_attempts <- 0L
  bound_iterations <- 0L
  series_terms <- 0L
  accept_N <- numeric(B)
  branch <- character(B)

  start_time <- proc.time()[["elapsed"]]

  for (i in seq_len(B)) {
    draw <- sample_one()
    values[i] <- draw$value
    proposals <- proposals + draw$proposals
    qP_attempts <- qP_attempts + draw$qP_attempts
    bound_iterations <- bound_iterations + draw$bound_iterations
    series_terms <- series_terms + draw$series_terms
    accept_N[i] <- draw$accept_N
    branch[i] <- draw$branch
  }

  elapsed <- proc.time()[["elapsed"]] - start_time

  data.frame(
    sampler = sampler_name,
    x = x,
    rho = rho,
    r = r,
    accepted_draws = B,
    branch = paste(unique(branch), collapse = ","),
    expected_proposals_per_draw = expected_fun(x = x, r = r),
    empirical_proposals_per_draw = proposals / B,
    avg_qP_attempts_per_draw = qP_attempts / B,
    runtime_seconds = elapsed,
    runtime_per_10000_draws = elapsed * 10000 / B,
    avg_bound_iterations_per_draw = bound_iterations / B,
    avg_series_terms_per_draw = series_terms / B,
    avg_accept_N = mean(accept_N),
    sample_mean = mean(values),
    expected_mean = expected_mean_f(x = x, r = r)
  )
}

# -----------------------------------------------------------------------------

all_results <- list()
row_id <- 1L

for (rho in rho_values) {
  cat("Running SH26 sampler for rho =", rho, "\n")
  all_results[[row_id]] <- benchmark_sampler(
    sampler_name = "SH26",
    sample_factory = make_sh26_sample_f_counted,
    expected_fun = sh26_expected_proposals,
    x = x_value,
    rho = rho,
    B = B
  )
  row_id <- row_id + 1L

  cat("Running current sampler for rho =", rho, "\n")
  all_results[[row_id]] <- benchmark_sampler(
    sampler_name = "Current",
    sample_factory = make_current_sample_f_counted,
    expected_fun = current_expected_proposals,
    x = x_value,
    rho = rho,
    B = B
  )
  row_id <- row_id + 1L
}

results <- do.call(rbind, all_results)

print(results)

results_file <- file.path(.project_root, "results", "comparison_f_results.csv")
dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
write.csv(results, file = results_file, row.names = FALSE)
cat("Results written to:", normalizePath(results_file), "\n")
