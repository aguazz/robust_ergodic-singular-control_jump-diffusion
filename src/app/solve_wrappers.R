# Solver, simulation, sweep, and cache wrappers for the companion Shiny app.

app_solver_control <- function(options) {
  list(
    ftol = options$ftol,
    xtol = options$xtol,
    maxit = options$maxit,
    stepmax = options$stepmax,
    trace = 0
  )
}

app_solve_solution <- function(params, guesses, options) {
  control <- app_solver_control(options)

  solve_optimal_barriers(
    p = params,
    xL0 = guesses$xL0,
    xk0 = guesses$xk0,
    xl0 = guesses$xl0,
    xU0 = guesses$xU0,
    control = control,
    verbose = FALSE,
    tol_regime_switch = options$tol_regime_switch,
    switch_to_newton = options$switch_to_newton,
    newton_control = control,
    record_iterates = options$record_iterates
  )
}

app_chosen_nleqslv_result <- function(sol) {
  if (identical(sol$chosen_solver, "newton") && !is.null(sol$nleqslv_refined)) {
    return(sol$nleqslv_refined)
  }
  sol$nleqslv_primary
}

app_max_residual <- function(sol) {
  ans <- app_chosen_nleqslv_result(sol)
  if (is.null(ans) || is.null(ans$fvec)) return(NA_real_)
  max(abs(as.numeric(ans$fvec)), na.rm = TRUE)
}

app_solver_message <- function(sol) {
  ans <- app_chosen_nleqslv_result(sol)
  msg <- ans$message %||% ""
  if (!nzchar(msg)) msg <- "No solver message returned."
  msg
}

app_simulate_solution <- function(params, sol, sim_options) {
  simulate_reflected_jd(
    T = sim_options$T,
    dt = sim_options$dt,
    x0 = sim_options$x0,
    params = params,
    thresholds = sol$x,
    seed = sim_options$seed
  )
}

app_run_sweep <- function(params, guesses, options, sweep_param, sweep_values) {
  sweep_param <- app_validate_sweep_param(sweep_param)
  fixed <- as.list(params[c("b", "delta", "r", "eps", "sigma", "mu", "u", "l")])

  solver_args <- list(
    control = app_solver_control(options),
    tol_regime_switch = options$tol_regime_switch,
    switch_to_newton = options$switch_to_newton,
    newton_control = app_solver_control(options)
  )

  comparative_sweeper(
    sweep_param = sweep_param,
    sweep_values = sweep_values,
    b = fixed$b,
    delta = fixed$delta,
    r = fixed$r,
    eps = fixed$eps,
    sigma = fixed$sigma,
    mu = fixed$mu,
    u = fixed$u,
    l = fixed$l,
    xL0 = guesses$xL0,
    xk0 = guesses$xk0,
    xl0 = guesses$xl0,
    xU0 = guesses$xU0,
    save = FALSE,
    solver_args = solver_args,
    verbose = FALSE
  )
}

app_solution_as_initial_guess <- function(sol, fallback) {
  if (is.null(sol) || is.null(sol$x)) return(fallback)
  list(
    xL0 = sol$x$xL,
    xk0 = sol$x$xk,
    xl0 = sol$x$xl,
    xU0 = sol$x$xU
  )
}
