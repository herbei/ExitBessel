suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
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
source_exit_bessel(root = .project_root)

.scratch_plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(file = .scratch_plot_file)
.scratch_plot_device <- grDevices::dev.cur()
on.exit({
  if (.scratch_plot_device %in% grDevices::dev.list()) {
    grDevices::dev.off(.scratch_plot_device)
  }
  unlink(.scratch_plot_file)
}, add = TRUE)

set.seed(20260818)

My_Theme <- theme(
  axis.title.x = element_text(size = 10),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 10),
  plot.title = element_text(size = 10, hjust = 0.5)
)

add_panel_tags <- function(plot_list) {
  Map(
    function(plot, tag) {
      plot +
        labs(tag = tag) +
        theme(
          plot.tag = element_text(size = 11, face = "bold"),
          plot.tag.position = c(0.02, 0.98)
        )
    },
    plot_list,
    LETTERS[seq_along(plot_list)]
  )
}

make_t_grid <- function(T, s, d_minus, MM = 1200L) {
  t_min <- max(.Machine$double.eps, d_minus^2 * 1e-4, T * 1e-7)
  dense_end <- min(T, max(10 * d_minus^2, 0.05 * T))

  grid <- seq(t_min, T, length.out = MM)
  if (dense_end > t_min && dense_end < T) {
    grid <- c(grid, seq(t_min, dense_end, length.out = MM))
  }

  switch_point <- if (s >= t_min && s <= T) s else numeric(0L)
  sort(unique(c(grid, switch_point)))
}

make_f_values <- function(t_vals, x, r, s, N_image = 20L, N_eigen = 20L) {
  f_vals <- numeric(length(t_vals))
  left <- t_vals <= s
  right <- t_vals > s

  if (any(left)) {
    bounds <- f_image_bounds(t = t_vals[left], x = x, r = r, N = N_image)
    f_vals[left] <- 0.5 * (bounds$lower + bounds$upper)
  }

  if (any(right)) {
    f_vals[right] <- 0.5 * (
      LN_f_E(t = t_vals[right], x = x, r = r, N = N_eigen) +
        UN_f_E(t = t_vals[right], x = x, r = r, N = N_eigen)
    )
  }

  pmax(f_vals, 0)
}

make_hist_df <- function(draws, T, bins = 180L) {
  breaks <- seq(0, T, length.out = bins + 1L)
  bins_for_draws <- findInterval(draws, breaks, rightmost.closed = TRUE)
  bins_for_draws <- pmin(bins_for_draws, bins + 1L)
  counts <- tabulate(bins_for_draws, nbins = bins + 1L)[seq_len(bins)]
  binwidth <- breaks[2L] - breaks[1L]

  data.frame(
    mid = 0.5 * (breaks[-1L] + breaks[-length(breaks)]),
    density = counts / (length(draws) * binwidth),
    width = binwidth
  )
}

# -----------------------------------------------------------------------------

make_positive_kde_df <- function(draws, T, n = 512L) {
  t_min <- max(.Machine$double.eps, T * 1e-7)
  t_vals <- seq(t_min, T, length.out = n)
  log_kde <- density(log(draws), n = 2048L)
  log_t <- log(t_vals)
  inside <- log_t >= min(log_kde$x) & log_t <= max(log_kde$x)
  kde_vals <- numeric(length(t_vals))

  if (any(inside)) {
    kde_vals[inside] <- approx(
      x = log_kde$x,
      y = log_kde$y,
      xout = log_t[inside]
    )$y / t_vals[inside]
  }

  data.frame(t = t_vals, kde = kde_vals)
}

# -----------------------------------------------------------------------------

make_panel_pair <- function(rho, x = 1, B = 50000L, MM = 1200L) {
  r <- rho * x
  lambda <- lambda_x(x)
  d_minus <- x - r
  envelope <- select_envelope_f(x = x, r = r)
  s <- envelope$s

  if (envelope$method == "I") {
    C_star <- x / r
    D_star <- D_I_f(s = s, r = r, x = x)
  } else {
    C_star <- x / (r * one_minus_Gamma_f(s = s, r = r, x = x))
    D_star <- D_P_f(s = s, r = r, x = x)
  }

  draws <- simulate_f(x = x, r = r, n = B)
  density_prob <- if (rho >= 0.95) 0.95 else 0.995
  density_floor <- if (rho >= 0.95) 20 * d_minus^2 else 0
  T_density <- max(
    as.numeric(quantile(draws, density_prob, names = FALSE)),
    density_floor,
    0.02 * s
  )
  T_envelope <- max(1.15 * s, as.numeric(quantile(draws, 0.995, names = FALSE)))

  t_env <- make_t_grid(T = T_envelope, s = s, d_minus = d_minus, MM = MM)
  f_env_vals <- make_f_values(t_vals = t_env, x = x, r = r, s = s)
  f_env_df <- data.frame(t = t_env, f = f_env_vals)

  t_left <- t_env[t_env <= s]
  if (envelope$method == "I") {
    env_left_vals <- C_star * levy_density(t_left, a = d_minus)
  } else {
    env_left_vals <- C_star * levy_density_pair_diff(t_left, a_lower = d_minus, r = r)
  }
  env_left_df <- data.frame(t = t_left, env = env_left_vals)

  t_right <- seq(s, T_envelope, length.out = max(2L, MM / 2L))
  env_right_df <- data.frame(
    t = t_right,
    env = envelope$CE * exp_density_rate(t_right, rate = lambda)
  )

  t_density <- make_t_grid(T = T_density, s = s, d_minus = d_minus, MM = MM)
  f_density_df <- data.frame(
    t = t_density,
    f = make_f_values(t_vals = t_density, x = x, r = r, s = s)
  )

  kde_df <- make_positive_kde_df(draws = draws, T = T_density)
  hist_df <- make_hist_df(draws = draws, T = T_density)

  env_max <- max(f_env_df$f, env_left_df$env, env_right_df$env, na.rm = TRUE)
  dens_max <- max(f_density_df$f, kde_df$kde, hist_df$density, na.rm = TRUE)
  miny_env <- -0.08 * env_max
  miny_density <- -0.08 * dens_max

  switch_y <- max(
    tail(env_left_df$env[is.finite(env_left_df$env)], 1L),
    head(env_right_df$env[is.finite(env_right_df$env)], 1L),
    na.rm = TRUE
  )

  title <- bquote(rho == .(rho) * "," ~~ J^"*" == .(as.name(envelope$method)))

  pp1 <- ggplot() +
    geom_line(data = f_env_df, aes(x = t, y = f), linewidth = 0.5, color = "red") +
    geom_line(data = env_left_df, aes(x = t, y = env), linewidth = 0.5, color = "blue", linetype = 2) +
    geom_line(data = env_right_df, aes(x = t, y = env), linewidth = 0.5, color = "purple", linetype = 2) +
    geom_point(aes(x = s, y = 0), size = 2, color = "red", inherit.aes = FALSE) +
    geom_segment(
      aes(x = s, y = 0, xend = s, yend = switch_y),
      color = "grey",
      linetype = 2
    ) +
    annotate("text", x = s, y = miny_env, label = 's^"*"', color = "black", parse = TRUE, size = 3) +
    labs(title = title, x = "t", y = "") +
    coord_cartesian(xlim = c(0, T_envelope), ylim = c(miny_env, 1.05 * env_max)) +
    theme_bw() +
    My_Theme

  pp2 <- ggplot() +
    geom_col(
      data = hist_df,
      aes(x = mid, y = density),
      width = unique(hist_df$width),
      fill = "steelblue",
      colour = "black",
      linewidth = 0.1,
      alpha = 0.5
    ) +
    geom_line(data = f_density_df, aes(x = t, y = f), linewidth = 0.5, color = "red") +
    geom_line(data = kde_df, aes(x = t, y = kde), linewidth = 0.5, color = "black") +
    labs(title = title, x = "t", y = "f(t; r, x)") +
    coord_cartesian(xlim = c(0, T_density), ylim = c(miny_density, 1.05 * dens_max)) +
    theme_bw() +
    My_Theme

  list(pp1, pp2)
}

rho_values <- c(0.01, 0.4, 0.8, 0.99)

figure_dir <- file.path(.project_root, "Figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

plot_list <- do.call(c, lapply(rho_values, make_panel_pair))
plot_list <- add_panel_tags(plot_list)

plot_f <- arrangeGrob(
  grobs = plot_list,
  ncol = 2
)

ggsave(
  filename = file.path(figure_dir, "Figure02.png"),
  plot = plot_f,
  width = 9,
  height = 12,
  units = "in",
  dpi = 600
)

cat("Figure written to:", normalizePath(file.path(figure_dir, "Figure02.png")), "\n")
