source(file.path("src", "load.R"))

p <- make_params(
  b = 0, delta = 1, r = 1, eps = 0.5,
  sigma = 1, mu = 1, u = 1, l = 1
)

xL0 <- -0.5
xk0 <- -0.1
xl0 <- 0.5
xU0 <- 1

subopt_H <- build_suboptimal_H(p, xL0, xk0, xl0, xU0)
print(diagnose(subopt_H))

sol_opt <- solve_optimal_barriers(
  p, xL0, xk0, xl0, xU0,
  verbose = interactive(),
  switch_to_newton = FALSE,
  tol_regime_switch = 1e-3
)
print(diagnose(sol_opt))

plot_H(
  sol_opt, p,
  show = interactive(), save = TRUE,
  top_blank = 0.075, bottom_blank = 0.05,
  margins = c(3, 2, 1, 1),
  show_x_axis_title = TRUE, show_y_axis_title = FALSE,
  axis_mgp = c(2.2, 0.7, 0), axis_title_cex = 1.4,
  tick_cex = 1
)

sim <- simulate_reflected_jd(params = p, thresholds = sol_opt$x, seed = 123)
plot_reflected_jd(sol_opt, p, sim = sim, show = interactive(), save = TRUE)
plot_controls(sim, show = interactive(), save = TRUE)
plot_reflected_with_controls(
  sol_opt, p, sim,
  top_blank = 0, bottom_blank = 0,
  heights = c(0.7, 1.7, 0.85), draw_legend = FALSE,
  margins = c(3, 2, 1, 1), axis_mgp = c(2.2, 0.7, 0),
  show_x_axis_title = TRUE, x_axis_title = "t",
  show_y_axis_title = FALSE, axis_title_cex = 1.4,
  tick_cex = 1,
  show = interactive(), save = TRUE
)
