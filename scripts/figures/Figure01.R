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

make_u_star_grid <- function(rho_values) {
  rows <- lapply(rho_values, function(rho) {
    envelope <- select_envelope_f(x = 1, r = rho)
    data.frame(
      rho = rho,
      u_star = envelope$u,
      method = envelope$method,
      u_I = envelope$u_I,
      u_P = envelope$u_P,
      mass = envelope$mass,
      mass_I = envelope$mass_I,
      mass_P = envelope$mass_P
    )
  })

  out <- do.call(rbind, rows)
  out$method <- factor(out$method, levels = c("P", "I"))
  out$segment <- cumsum(c(1L, diff(as.integer(out$method)) != 0L))

  out
}

rho_values <- seq(0.001, 0.999, length.out = 600L)
u_star_df <- make_u_star_grid(rho_values)

finite_pair <- is.finite(u_star_df$mass_P)
mass_gap <- u_star_df$mass_I - u_star_df$mass_P
switch_idx <- which(
  finite_pair[-length(finite_pair)] &
    finite_pair[-1L] &
    mass_gap[-length(mass_gap)] * mass_gap[-1L] <= 0
)
switch_df <- data.frame()
if (length(switch_idx) > 0L) {
  switch_df <- data.frame(
    rho = vapply(
      switch_idx,
      function(idx) {
        uniroot(
          f = function(rho) {
            envelope <- select_envelope_f(x = 1, r = rho)
            envelope$mass_I - envelope$mass_P
          },
          lower = u_star_df$rho[idx],
          upper = u_star_df$rho[idx + 1L],
          tol = sqrt(.Machine$double.eps)
        )$root
      },
      numeric(1L)
    )
  )
}

My_Theme <- theme(
  axis.title.x = element_text(size = 11),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 11),
  legend.position = "none"
)

plot_u_star <- ggplot(u_star_df, aes(x = rho, y = u_star)) +
  geom_line(aes(color = method, group = segment), linewidth = 0.8) +
  geom_hline(yintercept = uE_value(), color = "grey65", linewidth = 0.35, linetype = 2) +
  geom_vline(
    data = switch_df,
    aes(xintercept = rho),
    color = "grey50",
    linewidth = 0.35,
    linetype = 2
  ) +
  annotate(
    "text",
    x = 0.03,
    y = uE_value(),
    label = "u[E]",
    parse = TRUE,
    color = "grey35",
    hjust = 0,
    vjust = -0.6,
    size = 3
  ) +
  annotate(
    "text",
    x = 0.18,
    y = 0.625,
    label = "J^\"*\" == P",
    parse = TRUE,
    color = "blue",
    size = 4
  ) +
  annotate(
    "text",
    x = 0.76,
    y = 0.625,
    label = "J^\"*\" == I",
    parse = TRUE,
    color = "purple",
    size = 4
  ) +
  scale_color_manual(values = c(P = "blue", I = "purple"), name = "J*") +
  labs(
    x = expression(rho),
    y = expression(u^"*")
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(uE_value() - 0.02, 0.64)) +
  theme_bw() +
  My_Theme

figure_dir <- file.path(.project_root, "Figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(figure_dir, "Figure01.png"),
  plot = plot_u_star,
  width = 7,
  height = 13 / 6,
  units = "in",
  dpi = 600
)

cat("Figure written to:", normalizePath(file.path(figure_dir, "Figure01.png")), "\n")
