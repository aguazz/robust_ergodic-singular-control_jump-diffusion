source(file.path("src", "load.R"))

p <- make_params(
  b = 0, delta = 1, r = 1, eps = 0.5,
  sigma = 1, mu = 1, u = 1, l = 1
)

sol <- solve_optimal_barriers(
  p, xL0 = -0.5, xk0 = -0.1, xl0 = 0.5, xU0 = 1,
  verbose = FALSE,
  switch_to_newton = FALSE,
  tol_regime_switch = 1e-3
)

expected_x <- c(
  xL = -0.9900086,
  xk = -0.2214323,
  xl =  0.5642733,
  xU =  0.7511408
)
expected_gamma <- 2.48011711465

x <- unlist(sol$x)[names(expected_x)]
fmax <- max(abs(sol$nleqslv_primary$fvec))

stopifnot(
  all(abs(x - expected_x) < 1e-6),
  abs(sol$gamma - expected_gamma) < 1e-8,
  is.finite(fmax),
  fmax < 1e-8
)

cat("Smoke check passed.\n")
print(x)
cat("gamma:", format(sol$gamma, digits = 12), "\n")
cat("max residual:", format(fmax, digits = 6), "\n")
