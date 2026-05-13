# Simulation helpers for reflected jump-diffusion paths and pathwise costs.
#
# Code notation:
#   X      : reflected state path.
#   U      : cumulative upward reflection at the lower barrier xL.
#   L      : cumulative downward reflection at the upper barrier xU.
#            Some plots label this process D_t to emphasize downward control.
#   xL, xU : lower and upper reflecting barriers.
#   xk, xl : drift- and intensity-ambiguity thresholds.
#   dJ     : uncompensated compound-Poisson jump increment.
#
# Loading this file defines functions only.

.threshold_path <- function(value, tt, name) {
  if (is.null(value)) {
    stop(sprintf("thresholds$%s is required.", name), call. = FALSE)
  }

  n <- length(tt)
  if (is.function(value)) {
    out <- tryCatch(value(tt), error = function(e) NULL)
    if (is.null(out) || length(out) != n) {
      out <- vapply(tt, function(t) {
        y <- value(t)
        if (length(y) < 1) {
          stop(sprintf("threshold function '%s' returned length 0.", name), call. = FALSE)
        }
        as.numeric(y[1])
      }, numeric(1))
    }
  } else {
    out <- as.numeric(value)
  }

  out <- as.numeric(out)
  if (length(out) == 1L) {
    out <- rep(out, n)
  }
  if (length(out) != n) {
    stop(sprintf("thresholds$%s must be a scalar, a length-%d vector, or a function of time.",
                 name, n), call. = FALSE)
  }
  if (any(!is.finite(out))) {
    stop(sprintf("thresholds$%s contains non-finite values.", name), call. = FALSE)
  }
  out
}

# Simulate a reflected jump-diffusion path under bang-bang worst-case ambiguity.
simulate_reflected_jd <- function(T=8, dt=0.001, x0=NULL, params, thresholds,
                                  seed = NULL) {
  if (!is.numeric(T) || length(T) != 1L || !is.finite(T) || T <= 0) {
    stop("T must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt <= 0) {
    stop("dt must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.list(thresholds)) {
    stop("thresholds must be a list with xL, xU, xk, and xl.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  n <- as.integer(round(T/dt))
  if (n < 1L) stop("T/dt must define at least one simulation step.", call. = FALSE)
  tt <- seq(0, T, length.out=n+1)

  xL_path <- .threshold_path(thresholds$xL, tt, "xL")
  xU_path <- .threshold_path(thresholds$xU, tt, "xU")
  xk_path <- .threshold_path(thresholds$xk, tt, "xk")
  xl_path <- .threshold_path(thresholds$xl, tt, "xl")
  if (any(xL_path >= xU_path)) {
    stop("thresholds must satisfy xL < xU at every simulated time.", call. = FALSE)
  }

  xL0 <- xL_path[1]
  xU0 <- xU_path[1]
  if (is.null(x0)) x0 <- (xL0 + xU0)/2
  if (!is.numeric(x0) || length(x0) != 1L || !is.finite(x0)) {
    stop("x0 must be a finite scalar when provided.", call. = FALSE)
  }

  EY <- -1/params$mu
  X <- numeric(n+1)
  X[1] <- x0
  U_proc <- L_proc <- numeric(n+1)
  jumped <- logical(n+1)
  dJ <- numeric(n+1)

  for (i in seq_len(n)) {
    # Ambiguity is chosen from the current state and current thresholds.
    xk <- xk_path[i]
    xl <- xl_path[i]

    # Worst-case drift b + sigma*kappa, with kappa in {-delta, +delta}.
    kappa_wc <- if (X[i] <= xk) -params$delta else params$delta
    b_wc <- params$b + params$sigma * kappa_wc

    # Worst-case jump intensity.
    lam_wc <- if (X[i] <= xl) params$r * (1 + params$eps) else params$r * (1 - params$eps)
    lam_wc <- max(0, lam_wc)

    # Downward compound-Poisson jump increment.
    Nj <- rpois(1, lam_wc * dt)
    dJ[i+1] <- if (Nj > 0) -sum(rexp(Nj, rate=params$mu)) else 0

    # The jump increment dJ is uncompensated, while the baseline model drift
    # includes the compensation term for the reference intensity params$r.
    drift_wc <- b_wc - params$r * EY

    dW <- sqrt(dt)*rnorm(1)
    x_star <- X[i] + drift_wc * dt + params$sigma * dW + dJ[i+1]

    xL <- xL_path[i+1]
    xU <- xU_path[i+1]

    U_inc <- L_inc <- 0
    if (x_star < xL) { U_inc <- xL - x_star; X[i+1] <- xL }
    else if (x_star > xU) { L_inc <- x_star - xU; X[i+1] <- xU }
    else X[i+1] <- x_star

    U_proc[i+1] <- U_proc[i] + U_inc
    L_proc[i+1] <- L_proc[i] + L_inc
    jumped[i+1] <- Nj > 0
  }

  list(time=tt, X=X, U=U_proc, L=L_proc, jumped=jumped, dJ=dJ,
       xL=xL_path, xU=xU_path, xk=xk_path, xl=xl_path)
}

# Finite-horizon running-average cost along one simulated path.
# J(t) = (running-cost integral + u*U_t + l*L_t) / t.
finite_horizon_ergodic_cost <- function(sim, params, running = function(x) x^2) {
  required <- c("time", "X", "U", "L")
  missing <- setdiff(required, names(sim))
  if (length(missing)) {
    stop(sprintf("sim is missing required field(s): %s", paste(missing, collapse = ", ")),
         call. = FALSE)
  }

  t <- as.numeric(sim$time)
  X <- as.numeric(sim$X)
  U <- as.numeric(sim$U)
  L <- as.numeric(sim$L)

  len <- lengths(list(t, X, U, L))
  if (any(len != length(t))) {
    stop("sim$time, sim$X, sim$U, and sim$L must have the same length.", call. = FALSE)
  }

  dt <- diff(t)
  if (any(!is.finite(dt)) || any(dt <= 0)) stop("time grid must be increasing")

  g <- as.numeric(running(X))
  if (length(g) != length(X)) {
    stop("running(X) must return one value for each simulated state.", call. = FALSE)
  }
  g[!is.finite(g)] <- NA_real_

  # Trapezoid integral of the running cost.
  integ <- c(0, cumsum(0.5 * (g[-length(g)] + g[-1]) * dt))

  # Total cost including singular controls. U and L are cumulative processes.
  cum_cost <- integ + params$u * U + params$l * L

  avg_cost <- cum_cost / t
  avg_cost[1] <- NA_real_

  data.frame(time = t, cum = cum_cost, avg = avg_cost)
}
