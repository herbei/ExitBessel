suppressPackageStartupMessages({
  library(ggplot2)
})

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

set.seed(20260818)

B_env <- Sys.getenv("EXIT_BESSEL_RUNTIME_B", "10000")
B <- suppressWarnings(as.integer(B_env))
if (!is.finite(B) || B <= 0L) {
  stop("EXIT_BESSEL_RUNTIME_B must be a positive integer.", call. = FALSE)
}

x_value <- 1
rho_values <- sort(unique(c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)))

time_sh26_sampler <- function(x, r, B) {
  draws <- numeric(B)
  elapsed <- system.time({
    for (i in seq_len(B)) {
      draws[i] <- sh26_env$sample_from_f(x = x, r = r)
    }
  })[["elapsed"]]

  elapsed * 10000 / B
}

# -----------------------------------------------------------------------------

time_algorithm_4_sampler <- function(x, r, B) {
  elapsed <- system.time({
    draws <- current_env$simulate_f(x = x, r = r, n = B)
  })[["elapsed"]]

  elapsed * 10000 / B
}

# -----------------------------------------------------------------------------

benchmark_one_rho <- function(rho, x, B) {
  r <- rho * x

  gc(verbose = FALSE)
  sh26_seconds <- time_sh26_sampler(x = x, r = r, B = B)

  gc(verbose = FALSE)
  algorithm_4_seconds <- time_algorithm_4_sampler(x = x, r = r, B = B)

  data.frame(
    rho = rho,
    sampler = c("SH26 sample_from_f()", "Algorithm 4"),
    seconds_per_10000_draws = c(sh26_seconds, algorithm_4_seconds)
  )
}

# -----------------------------------------------------------------------------

results_file <- file.path(.project_root, "results", "f_runtime_curve_results.csv")
use_cached_results <- tolower(Sys.getenv(
  "EXIT_BESSEL_RUNTIME_USE_CACHE", "false"
)) %in% c("1", "true", "yes")

if (use_cached_results && file.exists(results_file)) {
  cat("Reading cached results from:", normalizePath(results_file), "\n")
  runtime_results <- read.csv(results_file)
} else {
  runtime_rows <- vector("list", length(rho_values))

  for (i in seq_along(rho_values)) {
    rho <- rho_values[[i]]
    cat("Benchmarking rho =", rho, "\n")
    runtime_rows[[i]] <- benchmark_one_rho(rho = rho, x = x_value, B = B)
  }

  runtime_results <- do.call(rbind, runtime_rows)
  dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
  write.csv(runtime_results, file = results_file, row.names = FALSE)
}

My_Theme <- theme(
  axis.title.x = element_text(size = 11),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 11),
  legend.position = "inside",
  legend.position.inside = c(0.78, 0.74),
  legend.background = element_rect(fill = "white", colour = "grey80"),
  legend.title = element_blank(),
  legend.text = element_text(size = 10),
  panel.grid.minor = element_blank()
)

runtime_plot <- ggplot(
  runtime_results,
  aes(x = rho, y = seconds_per_10000_draws, color = sampler)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(
    values = c(
      "Algorithm 4" = "purple",
      "SH26 sample_from_f()" = "grey25"
    ),
    labels = c(
      "Algorithm 4" = "Algorithm 4",
      "SH26 sample_from_f()" = "SH26"
    )
  ) +
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20", "50")
  ) +
  labs(
    x = expression(rho),
    y = "Seconds per 10,000 draws"
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw() +
  My_Theme

figure_dir <- file.path(.project_root, "Figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

figure_file <- file.path(figure_dir, "Figure03.png")
ggsave(
  filename = figure_file,
  plot = runtime_plot,
  width = 7,
  height = 8 / 3,
  units = "in",
  dpi = 600
)

cat("Results available at:", normalizePath(results_file), "\n")
cat("Figure written to:", normalizePath(figure_file), "\n")
