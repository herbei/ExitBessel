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

set.seed(20260817)

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

make_g_values <- function(t_vals, x, s, N_image = 20L, N_eigen = 20L) {
  g_vals <- numeric(length(t_vals))
  left <- t_vals <= s
  right <- t_vals > s

  g_vals[left] <- g_I_partial_sum(t = t_vals[left], x = x, N = N_image)
  g_vals[right] <- 0.5 * (
    LN_g_E(t = t_vals[right], x = x, N = N_eigen) +
      UN_g_E(t = t_vals[right], x = x, N = N_eigen)
  )

  g_vals
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

make_panel_pair <- function(x, T, B = 10000L, MM = 1000L) {
  u_star_g <- 0.403723
  s <- x^2 * u_star_g
  lam <- pi^2 / (2 * x^2)

  t_vals <- seq(0.001, T, length = MM)
  g_vals <- make_g_values(t_vals = t_vals, x = x, s = s)
  gDF <- data.frame(t = t_vals, g = g_vals)

  t_vals1 <- seq(0.001, s, length = MM / 2)
  t_vals2 <- seq(s, T, length = MM / 2)

  env_vals1 <- C_I_g(s = s, x = x) * g_I_partial_sum(t = t_vals1, x = x, N = 0L)
  env_vals2 <- C_E_g() * lam * exp(-lam * t_vals2)

  env1DF <- data.frame(t = t_vals1, env = env_vals1)
  env2DF <- data.frame(t = t_vals2, env = env_vals2)

  Z_final <- data.frame(Z = simulate_g(x = x, n = B))
  kdeDF <- make_positive_kde_df(draws = Z_final$Z, T = T)

  maxx <- T
  maxy <- 1.05 * max(gDF$g, env1DF$env, env2DF$env)
  miny <- -0.08 * maxy

  pp1 <- ggplot() +
    geom_line(data = gDF, aes(x = t, y = g), linewidth = 0.5, color = "red") +
    geom_line(data = env1DF, aes(x = t, y = env), linewidth = 0.5, color = "blue", linetype = 2) +
    geom_line(data = env2DF, aes(x = t, y = env), linewidth = 0.5, color = "purple", linetype = 2) +
    geom_point(aes(x = s, y = 0), size = 2, color = "red", inherit.aes = FALSE) +
    geom_segment(
      aes(x = s, y = 0, xend = s, yend = max(env1DF$env[MM / 2], env2DF$env[1])),
      color = "grey",
      linetype = 2
    ) +
    annotate("text", x = s, y = miny, label = 's[g]^"*"', color = "black", parse = TRUE, size = 3) +
    labs(title = paste("x =", x), x = "t", y = "") +
    coord_cartesian(xlim = c(0, maxx), ylim = c(miny, maxy)) +
    theme_bw() +
    My_Theme

  pp2 <- ggplot() +
    geom_histogram(
      data = Z_final,
      aes(x = Z, y = after_stat(density)),
      bins = 100,
      fill = "steelblue",
      colour = "black",
      linewidth = 0.1,
      alpha = 0.5
    ) +
    geom_line(data = gDF, aes(x = t, y = g), linewidth = 0.5, color = "red") +
    geom_line(data = kdeDF, aes(x = t, y = kde), linewidth = 0.5, color = "black") +
    labs(title = paste("x =", x), x = "t", y = "g(t;x)") +
    coord_cartesian(xlim = c(0, maxx), ylim = c(miny, maxy)) +
    theme_bw() +
    My_Theme

  list(pp1, pp2)
}

xvalues <- c(1, 4)
Tvalues <- c(1.5, 15)

figure_dir <- file.path(.project_root, "Figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

plot_list <- c(
  make_panel_pair(x = xvalues[1], T = Tvalues[1]),
  make_panel_pair(x = xvalues[2], T = Tvalues[2])
)
plot_list <- add_panel_tags(plot_list)

plot_g <- arrangeGrob(
  plot_list[[1]],
  plot_list[[2]],
  plot_list[[3]],
  plot_list[[4]],
  ncol = 2
)

ggsave(
  filename = file.path(figure_dir, "Figure04.png"),
  plot = plot_g,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600
)

cat("Figure written to:", normalizePath(file.path(figure_dir, "Figure04.png")), "\n")
