source(file.path("src", "load.R"))

p <- make_params(
  b = 0, delta = 1, r = 1, eps = 0.5,
  sigma = 1, mu = 1, u = 1, l = 1
)

sol <- solve_optimal_barriers(
  p, xL0 = -0.5, xk0 = -0.1, xl0 = 0.5, xU0 = 1,
  verbose = interactive(),
  switch_to_newton = FALSE,
  tol_regime_switch = 1e-3
)

cat("Optimal thresholds:\n")
print(unlist(sol$x))
cat("Ergodic value gamma:", format(sol$gamma, digits = 12), "\n")

if (interactive()) {
  plot_H(sol, p, show = TRUE, save = FALSE)
  plot_H_prime(sol, p, show = TRUE, save = FALSE)
}
