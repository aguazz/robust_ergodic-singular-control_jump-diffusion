# Simulation helpers for reflected jump-diffusion paths and pathwise costs.
# Source split from legacy/animations.R; supports constant or time-dependent thresholds.

# --- simulate reflected path (time-dependent boundaries)
simulate_reflected_jd <- function(T=8, dt=0.001, x0=NULL, params, thresholds,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- as.integer(round(T/dt)); tt <- seq(0, T, length.out=n+1)
  
  # --- NEW: allow time-dependent boundaries ---
  # --- NEW: ambiguity threshold functions (can be constants or functions) ---
  xL_raw <- thresholds$xL
  xU_raw <- thresholds$xU
  xk_raw <- thresholds$xk
  xl_raw <- thresholds$xl
  xk_fun <- if (is.function(xk_raw)) xk_raw else function(t) rep(as.numeric(xk_raw), length(t))
  xl_fun <- if (is.function(xl_raw)) xl_raw else function(t) rep(as.numeric(xl_raw), length(t))
  xL_fun <- if (is.function(xL_raw)) xL_raw else function(t) rep(as.numeric(xL_raw), length(t))
  xU_fun <- if (is.function(xU_raw)) xU_raw else function(t) rep(as.numeric(xU_raw), length(t))
  
  xL0 <- xL_fun(tt[1])[1]; xU0 <- xU_fun(tt[1])[1]
  if (is.null(x0)) x0 <- (xL0 + xU0)/2
  
  EY <- -1/params$mu
  X <- numeric(n+1); X[1] <- x0; U_proc <- L_proc <- numeric(n+1)
  jumped <- logical(n+1); dJ <- numeric(n+1)
  
  # (optional but handy) store boundaries over time
  xL_path <- numeric(n+1); xU_path <- numeric(n+1)
  xL_path[1] <- xL0; xU_path[1] <- xU0
  
  for (i in 1:n) {
    # ambiguity thresholds at current time (use tt[i] since we condition on X[i])
    xk <- xk_fun(tt[i])[1]
    xl <- xl_fun(tt[i])[1]
    
    # worst-case drift (bang-bang at xk)
    kappa_wc <- if (X[i] <= xk) params$b - params$sigma*params$delta else params$b + params$sigma*params$delta
    
    # worst-case intensity (bang-bang at xl)
    lam_wc <- if (X[i] <= xl) params$r * (1 + params$eps) else params$r * (1 - params$eps)
    lam_wc <- max(0, lam_wc)  # safety
    
    # jumps under worst-case intensity
    Nj <- rpois(1, lam_wc * dt)
    dJ[i+1] <- if (Nj > 0) -sum(rexp(Nj, rate=params$mu)) else 0
    
    # drift term consistent with your compensated simulation convention
    # drift_wc <- kappa_wc - lam_wc * EY
    # NEW (correct with uncompensated dJ)
    drift_wc <- kappa_wc - params$r * EY
    
    dW <- sqrt(dt)*rnorm(1)
    x_star <- X[i] + drift_wc * dt + params$sigma * dW + dJ[i+1]
    
    # --- NEW: boundary values at current step time ---
    xL <- xL_fun(tt[i+1])[1]
    xU <- xU_fun(tt[i+1])[1]
    xL_path[i+1] <- xL; xU_path[i+1] <- xU
    
    U_inc <- L_inc <- 0
    if (x_star < xL) { U_inc <- xL - x_star; X[i+1] <- xL }
    else if (x_star > xU) { L_inc <- x_star - xU; X[i+1] <- xU }
    else X[i+1] <- x_star
    
    U_proc[i+1] <- U_proc[i] + U_inc; L_proc[i+1] <- L_proc[i] + L_inc
    jumped[i+1] <- Nj > 0
  }
  
  list(time=tt, X=X, U=U_proc, L=L_proc, jumped=jumped, dJ=dJ,
       xL=xL_path, xU=xU_path)  # NEW (won't break old code)
}

# --- finite-horizon ergodic cost along ONE simulated path ---
# J(t) = (1/t) [ ∫_0^t running(X_s) ds + u U_t + l D_t ]
finite_horizon_ergodic_cost <- function(sim, params, running = function(x) x^2) {
  t <- as.numeric(sim$time)
  X <- as.numeric(sim$X)
  U <- as.numeric(sim$U)
  L <- as.numeric(sim$L)
  
  dt <- diff(t)
  if (any(!is.finite(dt)) || any(dt <= 0)) stop("time grid must be increasing")
  
  g <- running(X)
  g[!is.finite(g)] <- NA_real_
  
  # trapezoid integral of running cost
  integ <- c(0, cumsum(0.5 * (g[-length(g)] + g[-1]) * dt))
  # total cost including singular controls (U,D are cumulative)
  cum_cost <- integ + params$u * U + params$l * L
  
  avg_cost <- cum_cost / t
  avg_cost[1] <- NA_real_  # avoid 0/0 at t=0
  
  data.frame(time = t, cum = cum_cost, avg = avg_cost)
}
