source(file.path("src", "load.R"))

p <- make_params(
  b = 0, delta = 1, r = 1, eps = 0.5,
  sigma = 1, mu = 1, u = 1, l = 1
)

xL0 <- -0.5
xk0 <- -0.1
xl0 <- 0.5
xU0 <- 1

sol_opt <- solve_optimal_barriers(
  p, xL0, xk0, xl0, xU0,
  verbose = interactive(),
  record_iterates = FALSE,
  switch_to_newton = FALSE,
  tol_regime_switch = 1e-3
)

sim_time_varying <- simulate_reflected_jd(
  T = 75, dt = 0.25, params = p,
  thresholds = list(
    xL = function(t) sol_opt$x$xL + 0.05 * sin(t / 5),
    xU = function(t) sol_opt$x$xU + 0.1 * cos(t / 5),
    xk = sol_opt$x$xk,
    xl = sol_opt$x$xl
  ),
  seed = 666
)

save_animation_frames_canvas(
  sol = sol_opt, params = p, sim = sim_time_varying,
  out_dir = "frames/erg_sing_control",
  prefix = "erg_sing_control",
  every = 1, digits = 5,
  width_in = 12, height_in = 9, dpi = 80,
  show_ambiguity = FALSE, show_gamma = FALSE, show_barriers = TRUE,
  show_axis_labels = FALSE, show_x_axes = FALSE, show_y_axes = TRUE,
  bty = "n",
  draw_legend_reflected = FALSE, show = FALSE, save = TRUE,
  top_blank = 0.02, bottom_blank = 0.02, heights = c(0.95, 1.5, 1.1)
)

sim_robust <- simulate_reflected_jd(
  T = 75, dt = 0.25, params = p,
  thresholds = sol_opt$x,
  seed = 666
)

save_animation_frames_canvas(
  sol = sol_opt, params = p, sim = sim_robust,
  out_dir = "frames/robust_erg_sing_control",
  prefix = "robust_erg_sing_control",
  every = 1, digits = 5,
  width_in = 12, height_in = 9, dpi = 80,
  show_ambiguity = TRUE, show_gamma = TRUE, show_barriers = TRUE,
  show_axis_labels = FALSE, show_x_axes = FALSE, show_y_axes = TRUE,
  bty = "n",
  draw_legend_reflected = FALSE, show = FALSE, save = TRUE,
  top_blank = 0.02, bottom_blank = 0.02, heights = c(0.95, 1.5, 1.1)
)

sol_iter <- solve_optimal_barriers(
  p, xL0, xk0, xl0, xU0,
  verbose = interactive(),
  switch_to_newton = FALSE,
  tol_regime_switch = 1e-3,
  record_iterates = TRUE,
  record_xtol = 0
)

save_animation_frames_solver_convergence(
  sol_iter, p,
  out_dir = "frames/solver_convergence",
  prefix = "solver_convergence",
  every = 1, digits = 4,
  width_in = 12, height_in = 9, dpi = 150,
  show = FALSE, save = TRUE,
  show_axes = TRUE, bty = "n"
)
