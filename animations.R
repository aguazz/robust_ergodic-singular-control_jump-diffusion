# ========================= Frame-by-frame plotting helpers =========================
# Saves FOUR separate frame sequences:
#   1) Reflected state path  X̄_t  (with thresholds)
#   2) Lower control         D_t   (pushes down)
#   3) Upper control         U_t   (pushes up)
#   4) Finite-horizon ergodic cost  J_t = (1/t)[∫_0^t X_s^2 ds + u U_t + l D_t]

# --- colors ---
cols <- c(
  firebrick3      = "#CF2120",  # RGB(207, 33, 32)
  deepskyblue3    = "#009BCE",  # RGB(0, 155, 206)
  darkorange3     = "#CE6600",  # RGB(206, 102, 0)
  darkolivegreen3 = "#A3CF5B"   # RGB(163, 207, 91)
)


# --- mask future points so axes stay fixed across frames ---
.mask_after <- function(v, i) {
  v <- as.numeric(v)
  if (i < length(v)) v[(i+1):length(v)] <- NA_real_
  v
}

# --- PNG-only saver (fast; stable size for animations) ---
render_and_save_png <- function(
    fname_base, plotfun,
    dir="frames",
    width_in = 6.67, height_in = 4.67, dpi = 200,
    pointsize = NULL, family = NULL,
    match_current = TRUE,
    show = interactive(), save = TRUE
) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  if (match_current && dev.cur() != 1L) {
    sz <- dev.size("in")
    width_in  <- sz[1]; height_in <- sz[2]
    if (is.null(pointsize)) pointsize <- par("ps")
    if (is.null(family))    family    <- par("family")
  }
  if (is.null(pointsize)) pointsize <- 12
  fam <- if (is.null(family) || !nzchar(family)) "sans" else family
  
  if (isTRUE(show)) {
    op <- par(family = fam); on.exit(par(op), add=TRUE)
    plotfun()
  }
  if (!isTRUE(save)) return(invisible())
  
  png(file.path(dir, paste0(fname_base, ".png")),
      width = width_in, height = height_in, units = "in",
      res = dpi, pointsize = pointsize,
      type = getOption("bitmapType", "cairo"), antialias = "subpixel",
      bg = "transparent")
  op <- par(family = fam); on.exit(par(op), add=TRUE)
  plotfun()
  dev.off()
  invisible()
}

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

# Residuals of optimality equations with Regime-2 support (xl may exceed xU)
opt_conditions_residuals <- function(
    z, p, tol_build=1e-10, tol_regime_switch = 1e-3,
    return_sol = FALSE,
    recorder = NULL               # <-- NEW
) {
  xs  <- x_from_z(z)
  sol <- build_suboptimal_H(
    p,
    as.numeric(xs["xL"]), as.numeric(xs["xk"]),
    as.numeric(xs["xl"]), as.numeric(xs["xU"]),
    tol_regime_switch = tol_regime_switch
  )
  
  xL <- xs["xL"]; xk <- xs["xk"]; xl <- xs["xl"]; xU <- xs["xU"]
  R1 <- sol$pieces$region1; R2 <- sol$pieces$region2
  
  H1p <- function(x) gp_eval(R1$lam$lam_minus,x,R1$anchor_minus)*R1$u_minus +
    gp_eval(R1$lam$lam_plus ,x,R1$anchor_plus )*R1$u_plus  +
    qprime_eval(R1$q, x)
  
  H2p <- function(x) gp_eval(R2$lam$lam_minus,x,R2$anchor_minus)*R2$u_minus +
    gp_eval(R2$lam$lam_plus ,x,R2$anchor_plus )*R2$u_plus  +
    qprime_eval(R2$q, x)
  
  if (!isTRUE(sol$regime2)) {
    R3 <- sol$pieces$region3
    H3p <- function(x) gp_eval(R3$lam$lam_minus,x,R3$anchor_minus)*R3$u_minus +
      gp_eval(R3$lam$lam_plus ,x,R3$anchor_plus )*R3$u_plus  +
      qprime_eval(R3$q, x)
    
    r <- c(
      H1p(xL),
      H3p(xU),
      H1p(xk) - H2p(xk),
      H2p(xl) - H3p(xl)
    )
  } else {
    r <- c(
      H1p(xL),
      H2p(xU),
      H1p(xk) - H2p(xk),
      IH_at_xlambda_stable(sol, xla = xl)
    )
  }
  
  # ---------- NEW: record here (side-effect), sol is guaranteed available ----------
  if (is.function(recorder)) {
    # never let recording crash the solver
    try(recorder(z, r, sol), silent = TRUE)
  }
  
  if (isTRUE(return_sol)) attr(r, "sol") <- sol
  r
}

# ------------------------ Outer solver: Broyden → Newton ------------------------
# Switch to a pure Newton method (with FD Jacobian) when close to a solution.
# Nearness is detected by max|f| <= switch_ftol (and after at least switch_iter_min iterations).
# Uses the stable inner builder from Secs. 5.1–5.3.

solve_optimal_barriers <- function(
    p, xL0, xk0, xl0, xU0,
    method = "Broyden",
    control = list(ftol=1e-11, xtol=1e-11, maxit=1e6, stepmax=1, trace=0),
    tol_build = 1e-10,
    verbose = FALSE,
    tol_regime_switch = 1e-3,
    # --- hybrid controls ---
    switch_to_newton = FALSE,
    switch_ftol = 5e-2,          # threshold to trigger Newton (on max|f|)
    switch_iter_min = 2L,        # don't switch on the very first step
    newton_control = list(ftol=1e-12, xtol=1e-12, maxit=1e5, stepmax=1, trace=0),
    fd_step = 1e-6,               # relative step for finite-difference Jacobian
    record_iterates = FALSE,
    record_xtol = 1e-12
) {
  if (!requireNamespace("nleqslv", quietly=TRUE))
    stop("Please install.packages('nleqslv')")
  if (!(xL0 < xk0 && xk0 < xl0 && xk0 < xU0))
    stop("Initial guess must satisfy xL0 < xk0, and xk0 < xl0, xk0 < xU0 (xl0 may be > xU0).")
  
  # ---- reparam z <-> x (keeps barriers ordered) ----
  z0 <- z_from_x(xL0, xk0, xl0, xU0)
  
  # -------------------- iterate recorder (NEW) --------------------
  rec_env <- new.env(parent = emptyenv())
  rec_env$hist <- list()
  rec_env$last_z <- NULL
  rec_env$record_on <- TRUE
  rec_env$stage <- "broyden"
  
  # .c1_jumps_from_sol <- function(sol) {
  #   xs <- sol$x
  #   xL <- xs$xL; xk <- xs$xk; xl <- xs$xl; xU <- xs$xU
  #   
  #   # epsilon chosen relative to smallest gap
  #   gaps <- c(xk - xL, xl - xk, xU - min(xl, xU))
  #   gmin <- max(1e-8, min(gaps[gaps > 0], na.rm=TRUE))
  #   eps  <- min(1e-6 * (1 + max(abs(c(xL,xk,xl,xU)))), 0.25*gmin)
  #   
  #   Hp <- sol$Hp
  #   c(
  #     Hp_xL  = as.numeric(Hp(xL + eps)),
  #     dHp_xk = as.numeric(Hp(xk - eps) - Hp(xk + eps)),
  #     dHp_xl = as.numeric(Hp(xl - eps) - Hp(xl + eps)),
  #     Hp_xU  = as.numeric(Hp(xU - eps))
  #   )
  # }
  
  .record_step <- function(z, r, sol) {
    if (!isTRUE(record_iterates)) return(invisible())
    
    # FORCE a copy (do NOT use as.numeric(z) alone)
    znum <- c(unname(as.double(z)))   # c(...) forces a new vector
    
    if (!is.null(rec_env$last_z)) {
      dz <- suppressWarnings(max(abs(znum - rec_env$last_z)))
      if (is.finite(dz) && dz <= record_xtol) return(invisible())
    }
    
    rec_env$last_z <- znum   # already a fresh copy
    
    # old
    # jumps <- .c1_jumps_from_sol(sol)
    
    # new: store the true residuals that define fmax
    jumps <- c(
      Hp_xL  = as.numeric(r[1]),
      dHp_xk = as.numeric(r[3]),
      dHp_xl = as.numeric(r[4]),  # note: in regime2 this is NOT a derivative jump
      Hp_xU  = as.numeric(r[2])
    )
    
    xs <- sol$x
    xL <- as.numeric(xs$xL); xk <- as.numeric(xs$xk)
    xl <- as.numeric(xs$xl); xU <- as.numeric(xs$xU)
    
    rec_env$hist[[length(rec_env$hist) + 1L]] <- list(
      stage = rec_env$stage,
      z = znum,
      xL = xL, xk = xk, xl = xl, xU = xU,
      gamma = as.numeric(sol$gamma),
      fmax  = max(abs(as.numeric(r))),
      Hp_xL = jumps["Hp_xL"], dHp_xk = jumps["dHp_xk"],
      dHp_xl = jumps["dHp_xl"], Hp_xU = jumps["Hp_xU"]
    )
    invisible()
  }
  
  # Residuals (Eqs. (52)–(55) in stable form)
  .fn_eval <- function(z, do_record = TRUE) {
    
    rec_cb <- if (isTRUE(record_iterates) && isTRUE(do_record)) {
      function(z_loc, r_loc, sol_loc) {
        .record_step(z_loc, as.numeric(r_loc), sol_loc)
      }
    } else NULL
    
    r <- opt_conditions_residuals(
      z, p,
      tol_build = tol_build,
      tol_regime_switch = tol_regime_switch,
      return_sol = FALSE,         # <-- no longer needed
      recorder = rec_cb           # <-- key change
    )
    
    r_num <- as.numeric(r)
    
    if (isTRUE(verbose)) {
      xs <- x_from_z(z)
      cat(sprintf(
        "xL=%.6g  xk=%.6g  xl=%.6g  xU=%.6g  |  r=[% .3e, % .3e, % .3e, % .3e]\n",
        xs["xL"], xs["xk"], xs["xl"], xs["xU"], r_num[1], r_num[2], r_num[3], r_num[4]
      ))
    }
    r_num
  }
  
  fn       <- function(z) .fn_eval(z, do_record = TRUE)
  fn_norec <- function(z) .fn_eval(z, do_record = FALSE)
  
  # Stable, central-difference Jacobian (O(h^2))
  jac_fd <- function(z, fz = NULL) {
    n <- length(z); J <- matrix(0.0, n, n)
    if (is.null(fz)) fz <- fn_norec(z)
    
    for (j in 1:n) {
      h <- fd_step * (1 + abs(z[j]))
      zp <- z; zm <- z
      zp[j] <- z[j] + h;  zm[j] <- z[j] - h
      fp <- fn_norec(zp); fm <- fn_norec(zm)
      J[, j] <- (fp - fm) / (2*h)
    }
    J
  }
  
  # -------------------- Phase 1: robust Broyden --------------------
  ans_broyden <- nleqslv::nleqslv(
    x = z0, fn = fn, method = method,
    global = "dbldog", xscalm = "auto", control = control
  )
  
  if (isTRUE(record_iterates)) {
    # record the final point returned by nleqslv (in case it wasn't evaluated last)
    .fn_eval(ans_broyden$x, do_record = TRUE)
  }
  
  if (isTRUE(record_iterates)) {
    z_end  <- as.numeric(ans_broyden$x)
    xs_end <- x_from_z(z_end)
    sol_end <- build_suboptimal_H(
      p,
      as.numeric(xs_end["xL"]), as.numeric(xs_end["xk"]),
      as.numeric(xs_end["xl"]), as.numeric(xs_end["xU"]),
      tol_regime_switch = tol_regime_switch
    )
    .record_step(z_end, as.numeric(ans_broyden$fvec), sol_end)
  }
  
  rec_env$stage <- "broyden"
  
  # Decide whether to switch to Newton
  do_switch <- isTRUE(switch_to_newton)
  if (do_switch) {
    fmax   <- max(abs(ans_broyden$fvec %||% fn(ans_broyden$x)))
    iters  <- ans_broyden$iter %||% 0L
    do_switch <- is.finite(fmax) && (fmax <= switch_ftol) && (iters >= switch_iter_min)
  }
  
  # -------------------- Phase 2: pure Newton (optional) --------------------
  if (do_switch) {
    rec_env$stage <- "newton"
    z_start <- ans_broyden$x
    jf <- function(z) jac_fd(z)  # supply Jacobian explicitly
    ans_newton <- tryCatch(
      nleqslv::nleqslv(
        x = z_start, fn = fn, jac = jf, method = "Newton",
        global = "none", xscalm = "auto", control = newton_control
      ),
      error = function(e) e
    )
    
    # If Newton succeeds and improves max|f|, take it; otherwise keep Broyden result
    if (!inherits(ans_newton, "error")) {
      f_b <- max(abs(ans_broyden$fvec %||% fn(ans_broyden$x)))
      f_n <- max(abs(ans_newton$fvec  %||% fn(ans_newton$x)))
      final_ans <- if (is.finite(f_n) && (f_n <= f_b)) ans_newton else ans_broyden
      chosen <- if (identical(final_ans, ans_newton)) "newton" else "broyden"
    } else {
      final_ans <- ans_broyden
      chosen <- "broyden"
      warning(sprintf("Newton phase skipped (error: %s). Keeping Broyden solution.", conditionMessage(ans_newton)))
    }
  } else {
    final_ans <- ans_broyden
    chosen <- "broyden"
  }
  
  # -------------------- Build final suboptimal H and return --------------------
  xs_hat <- x_from_z(final_ans$x)
  sol_subopt <- build_suboptimal_H(
    p, xs_hat[["xL"]], xs_hat[["xk"]], xs_hat[["xl"]], xs_hat[["xU"]]
  )
  
  hist_df <- NULL
  if (isTRUE(record_iterates)) {
    if (length(rec_env$hist)) {
      hist_list <- lapply(seq_along(rec_env$hist), function(k) {
        h <- rec_env$hist[[k]]
        data.frame(
          step = k,
          stage = as.character(h$stage),
          xL = as.numeric(h$xL),
          xk = as.numeric(h$xk),
          xl = as.numeric(h$xl),
          xU = as.numeric(h$xU),
          gamma = as.numeric(h$gamma),
          fmax  = as.numeric(h$fmax),
          Hp_xL  = as.numeric(h$Hp_xL),
          dHp_xk = as.numeric(h$dHp_xk),
          dHp_xl = as.numeric(h$dHp_xl),
          Hp_xU  = as.numeric(h$Hp_xU),
          stringsAsFactors = FALSE
        )
      })
      
      hist_df <- do.call(rbind.data.frame, hist_list)
      rownames(hist_df) <- NULL
    } else {
      # keep as a 0-row DF WITH columns so downstream checks are meaningful
      hist_df <- data.frame(
        step=integer(), stage=character(),
        xL=numeric(), xk=numeric(), xl=numeric(), xU=numeric(),
        gamma=numeric(), fmax=numeric(),
        Hp_xL=numeric(), dHp_xk=numeric(), dHp_xl=numeric(), Hp_xU=numeric(),
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    H      = sol_subopt$H,
    Hp     = sol_subopt$Hp,
    gamma  = sol_subopt$gamma,
    pieces = sol_subopt$pieces,
    params = sol_subopt$params,
    x      = sol_subopt$x,
    nleqslv_primary = ans_broyden,
    nleqslv_refined = if (exists("ans_newton")) ans_newton else NULL,
    chosen_solver   = chosen,
    iter_history = hist_df
  )
}

# =================== Solver-convergence animation (NEW) ===================

# Precompute fixed axes (so frames don't jump)
solver_anim_limits <- function(history, params,
                               n_grid = 700,
                               pad_x_frac = 0.10,
                               top_blank = 0.10, bottom_blank = 0.05) {
  stopifnot(nrow(history) >= 1)
  
  # numeric safety
  history$xL <- as.numeric(history$xL)
  history$xU <- as.numeric(history$xU)
  history$xl <- as.numeric(history$xl)
  
  xL_min <- min(history$xL, na.rm=TRUE)
  x_max  <- max(pmax(history$xU, history$xl), na.rm=TRUE)
  
  if (!is.finite(xL_min) || !is.finite(x_max)) {
    stop("Non-finite x-range in iter_history. Check that iter_history has finite xL/xU/xl values.")
  }
  
  W <- x_max - xL_min
  if (!is.finite(W) || W <= 0) W <- 1  # fallback
  
  x_from <- xL_min - pad_x_frac * W
  x_to   <- x_max  + pad_x_frac * W
  
  xg <- seq(x_from, x_to, length.out = n_grid)
  
  # H range across *all* iterates
  ymin <- +Inf; ymax <- -Inf
  for (k in seq_len(nrow(history))) {
    sol_k <- build_suboptimal_H(params,
                                history$xL[k], history$xk[k], history$xl[k], history$xU[k])
    yk <- sol_k$H(xg)
    ymin <- min(ymin, yk, na.rm=TRUE)
    ymax <- max(ymax, yk, na.rm=TRUE)
  }
  
  frac <- 1 - top_blank - bottom_blank
  ylim_H <- c(ymin, ymax)
  ylim_H <- ylim_H + diff(ylim_H)/frac * c(-bottom_blank, top_blank)
  
  # gamma limits
  gy <- history$gamma
  ylim_g <- range(gy, finite=TRUE)
  if (diff(ylim_g) == 0) ylim_g <- ylim_g + c(-1,1) * 0.05 * max(1, abs(ylim_g[1]))
  ylim_g <- ylim_g + 0.06 * diff(ylim_g) * c(-1,1)
  
  list(xg=xg, x_from=x_from, x_to=x_to, ylim_H=ylim_H, ylim_g=ylim_g)
}

plot_frame_solver_convergence <- function(history, params, i,
                                          lims,
                                          out_dir="frames/solver", base="solver",
                                          digits=4,
                                          width_in=10.5, height_in=7.0, dpi=200,
                                          show=FALSE, save=TRUE,
                                          show_axes = TRUE,
                                          bty="n", box_lwd=1) {
  n <- nrow(history)
  stopifnot(i >= 1, i <= n)
  
  fname <- sprintf("%s_%0*d", base, digits, i)
  
  plotfun <- function() {
    op <- par(no.readonly=TRUE); on.exit(par(op), add=TRUE)
    
    # Layout:
    #   row 1: gamma (left) + text (right)
    #   row 2: H     (left) + legend (right)
    layout(matrix(c(1,2,
                    3,4), nrow=2, byrow=TRUE),
           widths=c(0.78, 0.22), heights=c(0.42, 0.58))
    
    # ---------------- Panel 1: gamma vs step ----------------
    par(mar=c(3.2, 3.4, 0.6, 0.3))
    k  <- seq_len(n)
    gy <- history$gamma
    gy_mask <- .mask_after(gy, i)
    
    plot(k, gy, type="n",
         xlab=if (isTRUE(show_axes)) "iteration" else "",
         ylab=if (isTRUE(show_axes)) expression(gamma) else "",
         xaxt=if (isTRUE(show_axes)) "s" else "n",
         yaxt=if (isTRUE(show_axes)) "s" else "n",
         ylim=lims$ylim_g, bty="n")
    
    lines(k, gy_mask, lwd=2)
    points(i, gy[i], pch=16, cex=1.5)
    
    if (!identical(bty, "n")) box(bty=bty, lwd=box_lwd)
    
    # ---------------- Panel 2: text column (C1 jumps) ----------------
    par(mar=c(3.2, 0.6, 0.4, 0.6))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    
    x_shift    <- 0.0
    y_shift    <- -0.1
    line_space <- 2.5
    cex_txt    <- 1.5
    
    x0 <- 0 + x_shift
    y0 <- 1 + y_shift
    dy <- 0.06 * line_space
    
    lines_expr <- list(
      bquote(step~.(i)~"/"~.(n)),
      # quote(""),
      bquote(H^minute*(underline(x)~"+")==" "*.(formatC(history$Hp_xL[i], format="e", digits=3))),
      bquote(Delta*H^minute*(x^kappa)==" "*.(formatC(history$dHp_xk[i], format="e", digits=3))),
      bquote(Delta*H^minute*(x^lambda)==" "*.(formatC(history$dHp_xl[i], format="e", digits=3))),
      bquote(H^minute*(bar(x)~"-")==" "*.(formatC(history$Hp_xU[i], format="e", digits=3)))
    )
    
    y <- y0
    for (e in lines_expr) {
      if (is.character(e) && e == "") {
        y <- y - dy
      } else {
        text(x0, y, labels=e, adj=c(0,1), xpd=NA, cex=cex_txt)
        y <- y - dy
      }
    }
    
    # ---------------- Panel 3: current H(x) ----------------
    sol_i <- build_suboptimal_H(params,
                                history$xL[i], history$xk[i], history$xl[i], history$xU[i])
    
    xs <- sol_i$x
    xg <- lims$xg
    yg <- sol_i$H(xg)
    
    par(mar=c(3.2, 3.4, 0.6, 0.3))
    plot(xg, yg, type="l", lwd=2,
         xlab=if (isTRUE(show_axes)) "x" else "",
         ylab=if (isTRUE(show_axes)) expression(H(x)) else "",
         xaxt=if (isTRUE(show_axes)) "s" else "n",
         yaxt=if (isTRUE(show_axes)) "s" else "n",
         xlim=c(lims$x_from, lims$x_to), ylim=lims$ylim_H, bty="n", cex.axis=1.5)
    
    # thresholds (consistent colors with your H plots if you like)
    cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22")
    ltys <- c(1,1,1,1)
    abline(v=xs$xL, col=cols[1], lty=ltys[1], lwd=2)
    abline(v=xs$xk, col=cols[2], lty=ltys[2], lwd=2)
    abline(v=xs$xl, col=cols[3], lty=ltys[3], lwd=2)
    abline(v=xs$xU, col=cols[4], lty=ltys[4], lwd=2)
    
    abline(h=-params$u, lty=2, col="gray50", lwd=2)
    abline(h=0,       lty=3, col="gray80", lwd=2)
    abline(h=params$l, lty=4, col="gray50", lwd=2)
    
    if (!identical(bty, "n")) box(bty=bty, lwd=box_lwd)
    
    
    # ---------------- Panel 4: legend (right of H, below text) ----------------
    par(mar=c(0.8, 0.6, 0.6, 0.6))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    
    # line/marker styles used in panels
    cols_thr <- c("#B22222","#1E90FF","#8A2BE2","#228B22")
    ltys_thr <- c(1,1,1,1)
    
    legend("center", bty="n", xpd=NA, cex=1.5, y.intersp=1.05,
           legend = c(expression(H(x)),
                      expression(underline(x)),
                      expression(x^kappa),
                      expression(x^lambda),
                      expression(bar(x)),
                      expression(-c[U]),
                      expression(c[D]),
                      expression(gamma)),
           lty    = c(1, ltys_thr[1], ltys_thr[2], ltys_thr[3], ltys_thr[4], 2, 4, 1),
           lwd    = c(2, 2, 2, 2, 2, 2, 2, 2),
           col    = c("black", cols_thr[1], cols_thr[2], cols_thr[3], cols_thr[4],
                      "gray50", "gray50", "black"),
           pch    = c(NA, NA, NA, NA, NA, NA, NA, 16),
           pt.cex = c(NA, NA, NA, NA, NA, NA, NA, 1.2))
    
    
  }
  
  render_and_save_png(fname, plotfun, dir=out_dir,
                      width_in=width_in, height_in=height_in, dpi=dpi,
                      show=show, save=save, match_current=FALSE)
}

save_animation_frames_solver_convergence <- function(sol_opt, params,
                                                     out_dir="frames/solver",
                                                     prefix="solver",
                                                     every=1L, digits=4,
                                                     width_in=10.5, height_in=7.0, dpi=200,
                                                     show=FALSE, save=TRUE,
                                                     show_axes=TRUE,
                                                     bty="n", box_lwd=1) {
  history <- sol_opt$iter_history
  
  # If something coerced it to a matrix, fix it
  if (is.matrix(history)) history <- as.data.frame(history, stringsAsFactors = FALSE)
  
  if (is.null(history)) stop("sol_opt$iter_history is NULL. Run solve_optimal_barriers(..., record_iterates=TRUE).")
  if (!is.data.frame(history)) stop(sprintf("iter_history must be a data.frame; got class: %s", paste(class(history), collapse=", ")))
  
  req <- c("xL","xk","xl","xU","gamma","fmax","stage","Hp_xL","dHp_xk","dHp_xl","Hp_xU")
  miss <- setdiff(req, names(history))
  if (length(miss)) stop("iter_history is missing columns: ", paste(miss, collapse=", "))
  
  # force numeric (except stage)
  num_cols <- setdiff(req, "stage")
  history[num_cols] <- lapply(history[num_cols], function(v) suppressWarnings(as.numeric(v)))
  
  # keep only finite barrier rows
  keep <- is.finite(history$xL) & is.finite(history$xU) & is.finite(history$xl) & is.finite(history$xk)
  history <- history[keep, , drop=FALSE]
  if (nrow(history) == 0) stop("iter_history has no finite (xL,xk,xl,xU). Recorder didn't capture valid iterates.")
  
  lims <- solver_anim_limits(history, params)
  
  n <- nrow(history)
  for (i in seq.int(1L, n, by=as.integer(every))) {
    plot_frame_solver_convergence(history, params, i, lims,
                                  out_dir=out_dir, base=prefix, digits=digits,
                                  width_in=width_in, height_in=height_in, dpi=dpi,
                                  show=show, save=save,
                                  show_axes=show_axes, bty=bty, box_lwd=box_lwd)
    if (isTRUE(show)) Sys.sleep(0.5)  # NEW
    }
  invisible(lims)
}

# ===================== 3-left + 1-right "canvas" frame =====================
# Layout:
#   Left column:  D_t (top), X̄_t (middle), U_t (bottom)
#   Right column: J_t (finite-horizon ergodic cost), same total height as left stack

# --- draw a step control panel in the CURRENT figure region ---
.draw_step_panel <- function(t, y, i, ylim, col, ylab_expr,
                             show_axis_labels,
                             bty, box_lwd,
                             show_x_axes = FALSE, show_y_axes = TRUE,
                             mar = c(0.4, 3.05, 0.4, 0.3)) {
  yp <- .mask_after(y, i)
  par(mar = mar)
  
  # xaxt0 <- if (isTRUE(show_x_axes)) "s" else "n"
  # yaxt0 <- if (isTRUE(show_y_axes)) "s" else "n"
  xaxt0 <- "n"
  yaxt0 <- "n"
  
  xlab0 <- if (isTRUE(show_axis_labels) && isTRUE(show_x_axes)) "time" else ""
  ylab0 <- if (isTRUE(show_axis_labels) && isTRUE(show_y_axes)) ylab_expr else ""
  
  # --- match x-width to reflected panel (reserve same right strip) ---
  xlim0 <- range(t)
  pad_x <- 0.12 * diff(xlim0)
  xlim  <- c(xlim0[1], xlim0[2] + pad_x)
  
  plot(t, yp, type="s", xlab=xlab0, ylab=ylab0, col=col,
       ylim=ylim, xlim=xlim, xaxt=xaxt0, yaxt=yaxt0, lwd=2, bty="n")
  
  if (isTRUE(show_x_axes)) {
    tmax <- max(t, finite=TRUE)
    axis(1, at=c(0, tmax), labels=c("0", "T"), cex.axis=1.5)
  }
  if (isTRUE(show_y_axes)) {
    at <- c(0, ylim[2])
    axis(2, at=at, labels=c("0", ""), cex.axis=1.5)   # max unlabeled
  }
  
  
  if (!identical(bty, "n")) box(bty=bty, lwd=box_lwd)
}

# --- draw reflected state panel in the CURRENT figure region ---
.draw_reflected_panel_frame <- function(sol, sim, i,
                                        show_ambiguity, show_barriers,
                                        show_axis_labels,
                                        show_x_axes, show_y_axes,  # <-- NEW
                                        bty, box_lwd,
                                        top_blank=0.075, bottom_blank=0.05,
                                        draw_legend=FALSE,
                                        mar = c(0.4, 3.05, 0.4, 0.3)) {
  xs <- sol$x # can be removed
  t  <- sim$time
  Xp <- .mask_after(sim$X, i)
  
  # NEW: time-dependent reflecting boundaries if available
  xL_vec <- if (!is.null(sim$xL)) as.numeric(sim$xL) else rep(as.numeric(xs$xL), length(t))
  xU_vec <- if (!is.null(sim$xU)) as.numeric(sim$xU) else rep(as.numeric(xs$xU), length(t))
  xLmin  <- min(xL_vec)
  xUmax  <- max(xU_vec)
  
  
  # y-limits
  Hcorr <- xUmax - xLmin
  frac  <- 1 - top_blank - bottom_blank
  Htot  <- Hcorr/frac
  ylim  <- c(xLmin - Htot*bottom_blank, xUmax + Htot*top_blank)
  
  par(mar = mar)
  
  xaxt0 <- "n"
  yaxt0 <- "n"
  
  # In stacked-left layout, usually no x-label here (bottom panel will carry it)
  xlab0 <- ""
  ylab0 <- if (isTRUE(show_axis_labels)) "" else ""
  
  # --- reserve a thin strip on the right for the cost scale (NEW) ---
  xlim0 <- range(t)
  pad_x <- 0.12 * diff(xlim0)
  xlim  <- c(xlim0[1], xlim0[2] + pad_x)
  
  plot(t, Xp, type="n", xlab=xlab0, ylab=ylab0,
       ylim=ylim, xlim=xlim, xaxt=xaxt0, yaxt=yaxt0, bty="n")   # <-- add xlim
  
  # --- quadratic-cost shading inside the control band + colorbar (REPLACE your gray rect) ---
  tmax <- xlim0[2]
  
  ybreaks <- seq(xLmin, xUmax, length.out = 101) # resolution of gradient
  ymid    <- (ybreaks[-1] + ybreaks[-length(ybreaks)]) / 2
  
  cmin <- if (xs$xL <= 0 && xs$xU >= 0) 0 else min(xs$xL^2, xs$xU^2)
  cmax <- max(xs$xL^2, xs$xU^2)
  z    <- (ymid^2 - cmin) / (cmax - cmin + 1e-12)               # normalize to [0,1]
  
  fill <- adjustcolor(gray(1 - 0.75*z), alpha.f = 0.45)         # greyscale gradient
  
  # NEW: gradient shading constrained to time-dependent band
  cmin <- if (xLmin <= 0 && xUmax >= 0) 0 else min(xLmin^2, xUmax^2)
  cmax <- max(xLmin^2, xUmax^2)
  z    <- (ymid^2 - cmin) / (cmax - cmin + 1e-12)
  
  fill <- adjustcolor(gray(1 - 0.75*z), alpha.f = 0.45)   # (keeps your colorbar fill)
  
  pal  <- adjustcolor(gray(1 - 0.75*seq(0,1,length.out=256)), alpha.f = 0.45)
  Z    <- matrix(rep(z, each=length(t)), nrow=length(t), ncol=length(ymid))
  Z[outer(xL_vec, ymid, `>`)] <- NA   # below xL(t)
  Z[outer(xU_vec, ymid, `<`)] <- NA   # above xU(t)
  
  image(t, ymid, Z, col=pal, add=TRUE, useRaster=TRUE)
  
  # --- add barriers and ambiguity thresholds
  cols4 <- c("#CE6600","#A3CF5B","#8A2BE2","#009BCE") # xL,xk,xl,xU
  ltys  <- c(2,2,2,2)
  
  if (isTRUE(show_barriers)) {
    lines(t, xL_vec, col=cols4[1], lty=ltys[1], lwd=3)
    lines(t, xU_vec, col=cols4[4], lty=ltys[4], lwd=3)
  }
  if (isTRUE(show_ambiguity)) {
    lines(t, rep(xs$xk, length(t)), col=cols4[2], lty=ltys[2], lwd=2)
    lines(t, rep(xs$xl, length(t)), col=cols4[3], lty=ltys[3], lwd=2)
  }
  
  # colorbar in the reserved strip
  x0 <- tmax + 0.35*pad_x
  x1 <- tmax + 0.55*pad_x
  rect(xleft=x0, xright=x1,
       ybottom=ybreaks[-length(ybreaks)], ytop=ybreaks[-1],
       col=fill, border=NA)
  rect(x0, xLmin, x1, xUmax, border="gray30", lwd=0.8)
  points((x0 + x1)/2, Xp[i], pch=16, col="black", cex=1.5, xpd=NA)  # NEW: marker on the scale
  
  # ticks labelled in cost units (c(x)=x^2)
  ct <- pretty(c(cmin, cmax), n=4)
  ct <- ct[ct >= cmin & ct <= cmax]
  sgn <- if (xUmax >= 0) 1 else -1
  yt  <- sgn*sqrt(ct)
  ok  <- yt >= xLmin & yt <= xUmax
  ct  <- ct[ok]; yt <- yt[ok]
  
  segments(x1, yt, x1 + 0.10*pad_x, yt, col="gray30", xpd=NA)
  text(x1 + 0.12*pad_x, yt, labels=formatC(ct, format="g", digits=3),
       adj=c(0,0.5), cex=1.5, xpd=NA)
  text((x0+x1)/2, xUmax, labels=expression(c(x)==x^2), pos=3, cex=1.5, xpd=NA)
  
  lines(t, Xp, lwd=1.2, col="black")
  points(t[i], Xp[i], pch=16, col="black", cex=1)  # NEW: highlight current state
  
  if (isTRUE(show_x_axes)) {
    tmax <- max(t, finite=TRUE)
    axis(1, at=c(0, tmax), labels=c("0", "T"), cex.axis=1.5)
  }
  
  if (isTRUE(show_y_axes)) {
    at   <- c(ylim[1], 0, ylim[2])
    labs <- c("", "0", "")
    keep <- at >= ylim[1] & at <= ylim[2]
    axis(2, at=at[keep], labels=labs[keep], cex.axis=1.5)
  }
  # usually only y-axis here on the left stack
  
  if (isTRUE(draw_legend)) {
    leg <- list(legend=c(expression(bar(X)[t])), lty=1, col="black", lwd=1.2)
    if (isTRUE(show_barriers)) {
      leg$legend <- c(leg$legend, expression(underline(x), bar(x)))
      leg$lty    <- c(leg$lty, 1, 1)
      leg$col    <- c(leg$col, cols4[1], cols4[4])
      leg$lwd    <- c(leg$lwd, 2, 2)
    }
    if (isTRUE(show_ambiguity)) {
      leg$legend <- c(leg$legend, expression(x^kappa, x^lambda))
      leg$lty    <- c(leg$lty, 2, 2)
      leg$col    <- c(leg$col, cols4[2], cols4[3])
      leg$lwd    <- c(leg$lwd, 2, 2)
    }
    legend("topleft", bty="n", horiz=TRUE, x.intersp=0.5, seg.len=2,
           legend=leg$legend, lty=leg$lty, col=leg$col, lwd=leg$lwd)
  }
  
  if (!identical(bty, "n")) box(bty=bty, lwd=box_lwd)
}

# --- draw cost panel in the CURRENT figure region ---
.draw_cost_panel <- function(cost_df, i, gamma, show_gamma,
                             show_axis_labels,
                             show_x_axes, show_y_axes,  # <-- NEW
                             bty, box_lwd,
                             mar = c(2.6, 3.05, 0.6, 0.8)) {
  t  <- cost_df$time
  Jp <- .mask_after(cost_df$avg, i)
  
  ylim <- range(cost_df$avg, finite=TRUE)
  if (!all(is.finite(ylim))) ylim <- c(0, 1)
  if (isTRUE(show_gamma) && is.finite(gamma)) ylim <- range(ylim, gamma)
  
  par(mar = mar)
  
  # xaxt0 <- if (isTRUE(show_x_axes)) "s" else "n"
  # yaxt0 <- if (isTRUE(show_y_axes)) "s" else "n"
  
  xaxt0 <- "n"
  yaxt0 <- "n"
  
  xlab0 <- if (isTRUE(show_axis_labels)) "time" else ""
  ylab0 <- if (isTRUE(show_axis_labels)) expression(J[t]) else ""
  
  plot(t, Jp, type="l", xlab=xlab0, ylab=ylab0, col="#CF2120",
       ylim=ylim, xaxt=xaxt0, yaxt=yaxt0, lwd=2, bty="n")
  
  if (isTRUE(show_x_axes)) {
    tmax <- max(t, finite=TRUE)
    axis(1, at=c(0, tmax), labels=c("0", "T"), cex.axis=1.5)
  }
  if (isTRUE(show_y_axes)) {
    at <- c(0, ylim[2])
    axis(2, at=at, labels=c("0", ""), cex.axis=1.5)   # max unlabeled
  }
  
  if (isTRUE(show_gamma) && is.finite(gamma)) abline(h=gamma, lty=3)
  if (!identical(bty, "n")) box(bty=bty, lwd=box_lwd)
}

.draw_legend_panel <- function(show_ambiguity = TRUE,
                               show_barriers  = TRUE,
                               show_gamma = TRUE,
                               bty="n", box_lwd=1,
                               mar=c(0.4, 0.4, 0.4, 0.4)) {
  par(mar = mar)
  plot.new()
  
  leg <- c("State process",
           "Cumulative downwards interventions",
           "Cumulative upwards interventions",
           "Running-average cost")
  col <- c("black", "#009BCE", "#CE6600", "#CF2120")
  lty <- c(1, 1, 1, 1)
  lwd <- c(2, 2, 2, 2)
  
  # barriers (xL(t), xU(t))
  if (isTRUE(show_barriers)) {
    leg <- c(leg, "Push-down barrier", "Push-up barrier")
    col <- c(col, "#CE6600", "#009BCE")   
    lty <- c(lty, 2, 2)
    lwd <- c(lwd, 2, 2)
  }
  
  # ambiguity thresholds (x^kappa, x^lambda)
  if (isTRUE(show_ambiguity)) {
    leg <- c(leg, "Drift ambiguity threshold", "Intensity ambiguity threshold")
    col <- c(col, "#A3CF5B", "#8A2BE2")   
    lty <- c(lty, 2, 2)
    lwd <- c(lwd, 2, 2)
  }
  
  # ambiguity thresholds (x^kappa, x^lambda)
  if (isTRUE(show_gamma)) {
    leg <- c(leg, "Ergodic value")
    col <- c(col, "black")   
    lty <- c(lty, 3)
    lwd <- c(lwd, 1)
  }
  
  legend("left", bty="n",
         legend = leg,
         col    = col,
         lty    = lty,
         lwd    = lwd,
         x.intersp = 0.6, y.intersp = 1.2,
         cex = par("cex") * 3,   # <-- NEW (tune 1.2–1.6)
         text.width = strwidth("Cumulative downwards interventions"))
  
  if (!identical(bty, "n")) box(bty=bty, lwd=box_lwd)
}

# --- single-frame canvas saver (one PNG per i) ---
plot_frame_canvas_3left_1right <- function(
    sol, params, sim, i = length(sim$time),
    cost_df = NULL,
    out_dir="frames/canvas", base="canvas", digits=5,
    width_in=10.5, height_in=5.5, dpi=200,
    # toggles
    show_ambiguity = TRUE, show_barriers = TRUE, show_gamma = TRUE,
    show_axis_labels = TRUE, show_x_axes = TRUE, show_y_axes = TRUE,
    bty = "o", box_lwd = 1,
    # layout tuning
    widths  = c(0.58, 0.42),
    heights = c(0.7, 1.7, 0.7),
    # per-panel margins
    mar_D  = c(0.2, 3.05, 0.4, 0.3),
    mar_X  = c(0.2, 3.05, 0.4, 0.3),
    mar_U  = c(2.6, 3.05, 0.4, 0.3),
    mar_J  = c(2.6, 3.05, 0.6, 0.8),
    draw_legend_reflected = FALSE,
    top_blank=0.075, bottom_blank=0.05,
    show=FALSE, save=TRUE
) {
  if (is.null(cost_df)) cost_df <- finite_horizon_ergodic_cost(sim, params)
  fname <- sprintf("%s_%0*d", base, digits, i)
  
  # fixed y-lims for controls across frames
  ylim_D <- range(sim$L, finite=TRUE); if (!all(is.finite(ylim_D))) ylim_D <- c(0, 1)
  ylim_U <- range(sim$U, finite=TRUE); if (!all(is.finite(ylim_U))) ylim_U <- c(0, 1)
  
  plotfun <- function() {
    op <- par(no.readonly=TRUE); on.exit(par(op), add=TRUE)
    
    # 3×2 layout:
    # row1: D (left) | LEGEND (right)
    # row2: X (left) | COST   (right, spans rows 2-3)
    # row3: U (left) | COST
    layout(matrix(c(1, 4,
                    2, 5,
                    3, 5), nrow=3, byrow=TRUE),
           widths = widths, heights = heights)
    
    # Panel 1 (top-left): D_t (your sim$L)
    .draw_step_panel(sim$time, sim$L, i, ylim=ylim_D, col="#009BCE",
                     ylab_expr=expression(D[t]),
                     show_axis_labels=show_axis_labels, 
                     show_x_axes=show_x_axes, show_y_axes=show_y_axes,
                     bty=bty, box_lwd=box_lwd,
                     mar=mar_D)
    
    
    # Panel 2 (middle-left): X̄_t
    .draw_reflected_panel_frame(sol, sim, i,
                                show_ambiguity=show_ambiguity,
                                show_barriers=show_barriers,
                                show_axis_labels=show_axis_labels,
                                show_x_axes=show_x_axes, show_y_axes=show_y_axes,  # <-- changed
                                bty=bty, box_lwd=box_lwd,
                                top_blank=top_blank, bottom_blank=bottom_blank,
                                draw_legend=draw_legend_reflected,
                                mar=mar_X)
    
    # Panel 3 (bottom-left): U_t
    .draw_step_panel(sim$time, sim$U, i, ylim=ylim_U, col="#CE6600",
                     ylab_expr=expression(U[t]),
                     show_axis_labels=show_axis_labels, 
                     show_x_axes=TRUE, show_y_axes=show_y_axes,
                     bty=bty, box_lwd=box_lwd,
                     mar=mar_U)
    
    # Panel 4 (top-right): legend
    .draw_legend_panel(bty="n",show_ambiguity=show_ambiguity, show_barriers=show_barriers,
                       show_gamma=show_gamma)
    
    # Panel 5 (right, rows 2-3): J_t
    .draw_cost_panel(cost_df, i, gamma=sol$gamma %||% NA_real_, show_gamma=show_gamma,
                     show_axis_labels=show_axis_labels, 
                     show_x_axes=TRUE, show_y_axes=show_y_axes,  # <-- changed
                     bty=bty, box_lwd=box_lwd,
                     mar=mar_J)
  }
  
  render_and_save_png(fname, plotfun, dir=out_dir,
                      width_in=width_in, height_in=height_in, dpi=dpi,
                      show=show, save=save, match_current=FALSE)
}

# --- batch exporter for the canvas frames ---
save_animation_frames_canvas <- function(
    sol, params, sim,
    out_dir="frames/canvas", prefix="canvas",
    every=1L, digits=5,
    width_in=10.5, height_in=5.5, dpi=200,
    show_ambiguity=TRUE, show_barriers=TRUE, show_gamma=TRUE,
    show_axis_labels=TRUE, 
    show_x_axes=TRUE, show_y_axes=TRUE, 
    bty="o", box_lwd=1,
    widths=c(0.58,0.42), heights=c(0.7,1.7,0.7),
    draw_legend_reflected=FALSE,
    top_blank=0.075, bottom_blank=0.05,
    show=FALSE, save=TRUE
) {
  cost_df <- finite_horizon_ergodic_cost(sim, params)
  n <- length(sim$time)
  
  for (i in seq.int(2L, n, by=as.integer(every))) {
    plot_frame_canvas_3left_1right(
      sol, params, sim, i,
      cost_df = cost_df,
      out_dir = out_dir,
      base    = prefix,
      digits  = digits,
      width_in=width_in, height_in=height_in, dpi=dpi,
      show_ambiguity=show_ambiguity, show_barriers=show_barriers, show_gamma=show_gamma,
      show_axis_labels=show_axis_labels, 
      show_x_axes=show_x_axes, show_y_axes=show_y_axes, 
      bty=bty, box_lwd=box_lwd,
      widths=widths, heights=heights,
      draw_legend_reflected=draw_legend_reflected,
      top_blank=top_blank, bottom_blank=bottom_blank,
      show=show, save=save
    )
    if (isTRUE(show)) Sys.sleep(0.5)  # NEW
  }
  invisible(list(cost=cost_df))
}


optimize_png_frames_magick <- function(in_dir, out_dir,
                                       scale = 0.7,   # 0.7 means 70% linear size
                                       colors = 256)  # try 128 if you want smaller
{
  if (!requireNamespace("magick", quietly = TRUE))
    stop("Install magick: install.packages('magick')")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  files <- list.files(in_dir, pattern="\\.png$", full.names=TRUE)
  for (f in files) {
    img <- magick::image_read(f)
    
    if (!is.null(scale) && is.finite(scale) && scale > 0 && scale != 1) {
      img <- magick::image_scale(img, paste0(round(scale*100), "%"))
    }
    
    # quantize to palette (keeps alpha reasonably well)
    img <- magick::image_quantize(img, max = colors)
    
    magick::image_write(img, path = file.path(out_dir, basename(f)),
                        format = "png")
  }
}



# ================================ Example =================================
# Parameters
p <- make_params(b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1, l = 1)
xL0 <- -0.5; xk0 <- -0.1; xl0 <- 0.5; xU0 <- 1

# (Assumes you already built sol_opt and sim as in your script.)
# Optimal H (outer solver)
sol_opt <- solve_optimal_barriers(p, xL0, xk0, xl0, xU0, verbose=TRUE, record_iterates = FALSE,
                                  switch_to_newton = FALSE, tol_regime_switch = 1e-3)

# 1
# # sim     <- simulate_reflected_jd(T=10, dt = 0.5, params=p, thresholds=sol_opt$x, seed=123)
sim <- simulate_reflected_jd(
  T=75, dt=0.25, params=p,
  thresholds=list(
    xL=function(t) sol_opt$x$xL + 0.05*sin(t/5),
    xU=function(t) sol_opt$x$xU + 0.1*cos(t/5),
    xk=sol_opt$x$xk,
    xl=sol_opt$x$xl
  ),
  seed=666
)

  # After you have sol_opt and sim:
save_animation_frames_canvas(
  sol=sol_opt, params=p, sim=sim,
  out_dir="frames/erg_sing_control",
  prefix="erg_sing_control",
  every=1, digits=5,
  width_in=12, height_in=9, dpi=80,
  show_ambiguity=FALSE, show_gamma=FALSE, show_barriers=TRUE,
  show_axis_labels=FALSE,
  show_x_axes=FALSE, show_y_axes=TRUE,
  bty="n",
  draw_legend_reflected=FALSE, show = FALSE, save = TRUE,
  top_blank = 0.02, bottom_blank = 0.02, heights = c(0.95, 1.5, 1.1)
)

# optimize_png_frames_magick("frames/erg_sing_control", "frames/erg_sing_control_opt", scale=0.5, colors=2048)

# 2
sim     <- simulate_reflected_jd(T=75, dt = 0.25, params=p, 
                                 thresholds=sol_opt$x, seed=666)

# After you have sol_opt and sim:
save_animation_frames_canvas(
  sol=sol_opt, params=p, sim=sim,
  out_dir="frames/robust_erg_sing_control",
  prefix="robust_erg_sing_control",
  every=1, digits=5,
  width_in=12, height_in=9, dpi=80,
  show_ambiguity=TRUE, show_gamma=TRUE, show_barriers=TRUE,
  show_axis_labels=FALSE,
  show_x_axes=FALSE, show_y_axes=TRUE,
  bty="n",
  draw_legend_reflected=FALSE, show = FALSE, save = TRUE,
  top_blank = 0.02, bottom_blank = 0.02, heights = c(0.95, 1.5, 1.1)
)

# --- Solver convergence
sol_opt <- solve_optimal_barriers(
  p, xL0, xk0, xl0, xU0,
  verbose=TRUE,
  switch_to_newton = FALSE,
  tol_regime_switch = 1e-3,
  record_iterates=TRUE, record_xtol=0
)

save_animation_frames_solver_convergence(
  sol_opt, p,
  out_dir="frames/solver_convergence",
  prefix="solver_convergence",
  every=1, digits=4,
  width_in=12, height_in=9, dpi=150,
  show=TRUE, save=TRUE,
  show_axes=TRUE, bty="n"
)
