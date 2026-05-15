# Validation helpers for the companion Shiny app.
#
# These functions keep Shiny input checking out of the solver. They return
# plain R objects, so they are also useful from scripts when debugging the app.

app_default_model_params <- function() {
  list(
    b = 0,
    delta = 1,
    r = 1,
    eps = 0.5,
    sigma = 1,
    mu = 1,
    u = 1,
    l = 1
  )
}

app_default_initial_guess <- function() {
  list(xL0 = -0.5, xk0 = -0.1, xl0 = 0.5, xU0 = 1)
}

app_default_solver_options <- function() {
  list(
    ftol = 1e-11,
    xtol = 1e-11,
    maxit = 1e5,
    stepmax = 1,
    tol_regime_switch = 1e-3,
    switch_to_newton = FALSE,
    record_iterates = TRUE
  )
}

app_default_simulation_options <- function() {
  list(T = 8, dt = 0.001, seed = 123, x0 = NULL, max_steps = 50000)
}

app_numeric_scalar <- function(x, name,
                               lower = -Inf, upper = Inf,
                               lower_strict = FALSE,
                               upper_strict = FALSE) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L || !is.finite(x)) {
    stop(sprintf("%s must be a finite number.", name), call. = FALSE)
  }

  lower_ok <- if (isTRUE(lower_strict)) x > lower else x >= lower
  upper_ok <- if (isTRUE(upper_strict)) x < upper else x <= upper

  if (!lower_ok || !upper_ok) {
    lower_op <- if (isTRUE(lower_strict)) ">" else ">="
    upper_op <- if (isTRUE(upper_strict)) "<" else "<="
    if (is.finite(lower) && is.finite(upper)) {
      stop(
        sprintf("%s must satisfy %s %s %g and %s %g.",
                name, name, lower_op, lower, upper_op, upper),
        call. = FALSE
      )
    }
    if (is.finite(lower)) {
      stop(sprintf("%s must be %s %g.", name, lower_op, lower), call. = FALSE)
    }
    if (is.finite(upper)) {
      stop(sprintf("%s must be %s %g.", name, upper_op, upper), call. = FALSE)
    }
  }

  x
}

app_integer_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  x <- app_numeric_scalar(x, name, lower = lower, upper = upper)
  x <- as.integer(round(x))
  if (!is.finite(x)) stop(sprintf("%s must be an integer.", name), call. = FALSE)
  x
}

app_make_params <- function(b, delta, r, eps, sigma, mu, u, l) {
  vals <- list(
    b = app_numeric_scalar(b, "b"),
    delta = app_numeric_scalar(delta, "delta", lower = 0),
    r = app_numeric_scalar(r, "r", lower = 0, lower_strict = TRUE),
    eps = app_numeric_scalar(eps, "epsilon", lower = 0, upper = 1),
    sigma = app_numeric_scalar(sigma, "sigma", lower = 0, lower_strict = TRUE),
    mu = app_numeric_scalar(mu, "mu", lower = 0, lower_strict = TRUE),
    u = app_numeric_scalar(u, "c_U", lower = 0),
    l = app_numeric_scalar(l, "c_D", lower = 0)
  )

  do.call(make_params, vals)
}

app_validate_initial_guess <- function(xL0, xk0, xl0, xU0) {
  guesses <- list(
    xL0 = app_numeric_scalar(xL0, "underline{x}_0"),
    xk0 = app_numeric_scalar(xk0, "xk0"),
    xl0 = app_numeric_scalar(xl0, "xl0"),
    xU0 = app_numeric_scalar(xU0, "overline{x}_0")
  )

  if (!(guesses$xL0 < guesses$xk0 &&
        guesses$xk0 < guesses$xl0 &&
        guesses$xk0 < guesses$xU0)) {
    stop(
      "Initial guesses must satisfy underline{x}_0 < xk0, xk0 < xl0, and xk0 < overline{x}_0.",
      call. = FALSE
    )
  }

  guesses
}

app_validate_solver_options <- function(ftol, xtol, maxit, stepmax,
                                        tol_regime_switch,
                                        switch_to_newton = FALSE,
                                        record_iterates = FALSE) {
  list(
    ftol = app_numeric_scalar(ftol, "ftol", lower = 0, lower_strict = TRUE),
    xtol = app_numeric_scalar(xtol, "xtol", lower = 0, lower_strict = TRUE),
    maxit = app_integer_scalar(maxit, "maxit", lower = 1),
    stepmax = app_numeric_scalar(stepmax, "stepmax", lower = 0, lower_strict = TRUE),
    tol_regime_switch = app_numeric_scalar(
      tol_regime_switch,
      "tol_regime_switch",
      lower = 0,
      lower_strict = TRUE
    ),
    switch_to_newton = isTRUE(switch_to_newton),
    record_iterates = isTRUE(record_iterates)
  )
}

app_validate_simulation_options <- function(T, dt, seed, use_midpoint,
                                            x0 = NULL, max_steps = 50000) {
  T <- app_numeric_scalar(T, "T", lower = 0, lower_strict = TRUE)
  dt <- app_numeric_scalar(dt, "dt", lower = 0, lower_strict = TRUE)
  n_steps <- as.integer(round(T / dt))
  max_steps <- app_integer_scalar(max_steps, "max simulation steps", lower = 1)

  if (n_steps < 1L) {
    stop("T / dt must define at least one simulation step.", call. = FALSE)
  }
  if (n_steps > max_steps) {
    stop(
      sprintf("The simulation would use %s steps; reduce T, increase dt, or raise the limit.",
              format(n_steps, big.mark = ",")),
      call. = FALSE
    )
  }

  seed <- app_integer_scalar(seed, "seed")
  x0 <- if (isTRUE(use_midpoint)) NULL else app_numeric_scalar(x0, "x0")

  list(T = T, dt = dt, seed = seed, x0 = x0, n_steps = n_steps)
}

app_validate_sweep_values <- function(from, to, n, max_points = 80) {
  from <- app_numeric_scalar(from, "sweep lower bound")
  to <- app_numeric_scalar(to, "sweep upper bound")
  n <- app_integer_scalar(n, "sweep grid size", lower = 2, upper = max_points)

  if (from == to) {
    stop("Sweep bounds must be distinct.", call. = FALSE)
  }

  seq(from, to, length.out = n)
}

app_validate_sweep_param <- function(sweep_param) {
  sweep_param <- as.character(sweep_param)[1]
  if (!sweep_param %in% .sweep_known_params) {
    stop(
      sprintf("sweep_param must be one of: %s", paste(.sweep_known_params, collapse = ", ")),
      call. = FALSE
    )
  }
  sweep_param
}

app_filename_tag <- function(x) {
  x <- gsub("\\s+", "_", as.character(x))
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("_+", "_", x)
}
