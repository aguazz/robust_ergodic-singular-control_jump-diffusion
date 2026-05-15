# Plot wrappers for Shiny outputs.
#
# The app reuses the repo's base-graphics plotting functions where possible and
# adds a few lightweight diagnostics that are app-specific.

app_expand_range <- function(x, frac = 0.06, fallback = c(-1, 1)) {
  rng <- range(x, finite = TRUE)
  if (!all(is.finite(rng))) rng <- fallback
  if (diff(rng) == 0) {
    bump <- max(1, abs(rng[1]))
    rng <- rng + c(-1, 1) * 0.05 * bump
  }
  rng + c(-1, 1) * diff(rng) * frac
}

app_solution_x_grid <- function(sol, n = 800, pad_x_frac = 0.12) {
  xs <- unlist(sol$x)[c("xL", "xk", "xl", "xU")]
  x_max <- max(xs[["xU"]], xs[["xl"]])
  width <- x_max - xs[["xL"]]
  if (!is.finite(width) || width <= 0) width <- 1
  seq(xs[["xL"]] - pad_x_frac * width,
      x_max + pad_x_frac * width,
      length.out = n)
}

app_safe_plot_par <- function(expr, mfrow = NULL, mar = NULL) {
  saved_names <- c("mfrow", "mar", "mgp", "tcl", "xaxs", "yaxs",
                   "cex", "cex.axis", "cex.lab", "cex.main")
  op <- par(saved_names)
  on.exit(par(op), add = TRUE)
  if (!is.null(mfrow)) par(mfrow = mfrow)
  if (!is.null(mar)) par(mar = mar)
  force(expr)
}

app_threshold_style <- function() {
  style <- .threshold_style()
  style$ltys[c("xL", "xU")] <- 2
  style
}

app_draw_threshold_lines <- function(xs, ylim, include_labels = FALSE) {
  style <- app_threshold_style()
  abline(v = xs, col = style$cols, lty = style$ltys, lwd = 1.8)
  if (isTRUE(include_labels)) {
    labs <- expression(underline(x), x^kappa, x^lambda, bar(x))
    text(
      x = xs,
      y = ylim[2],
      labels = labs,
      col = style$cols,
      pos = 3,
      xpd = NA,
      cex = 0.85
    )
  }
  invisible(NULL)
}

app_plot_H <- function(sol, params) {
  xg <- app_solution_x_grid(sol)
  yg <- sol$H(xg)
  xs <- unlist(sol$x)[c("xL", "xk", "xl", "xU")]
  style <- app_threshold_style()
  ylim <- app_expand_range(c(yg, -params$u, params$l), frac = 0.08)

  app_safe_plot_par({
    par(mar = c(3.0, 3.2, 0.8, 0.8), mgp = c(2.1, 0.65, 0),
        tcl = -0.25, xaxs = "i")
    plot(
      xg, yg,
      type = "l",
      lwd = 1.7,
      col = "black",
      xlab = "x",
      ylab = expression(H(x)),
      ylim = ylim
    )
    app_draw_threshold_lines(xs, ylim)
    abline(h = -params$u, lty = 3, col = "#777777")
    abline(h = 0, lty = 3, col = "#BBBBBB")
    abline(h = params$l, lty = 4, col = "#777777")
    legend(
      "topleft",
      legend = expression(H(x), underline(x), x^kappa, x^lambda, bar(x), -c[U], c[D]),
      lty = c(1, style$ltys, 3, 4),
      col = c("black", style$cols, "#777777", "#777777"),
      lwd = c(1.7, rep(1.8, 4), 1.2, 1.2),
      bty = "n",
      cex = 0.95,
      ncol = 1
    )
    grid(col = "gray90")
    box()
  })
}

app_plot_H_prime <- function(sol, params) {
  xg <- app_solution_x_grid(sol)
  yg <- sol$Hp(xg)
  xs <- unlist(sol$x)[c("xL", "xk", "xl", "xU")]
  style <- app_threshold_style()
  for (b in xs) {
    j <- which.min(abs(xg - b))
    for (k in c(j - 1L, j, j + 1L)) {
      if (k >= 1L && k <= length(yg)) yg[k] <- NA_real_
    }
  }
  ylim <- app_expand_range(c(yg, -params$u, params$l), frac = 0.08)

  app_safe_plot_par({
    par(mar = c(3.0, 3.2, 0.8, 0.8), mgp = c(2.1, 0.65, 0),
        tcl = -0.25, xaxs = "i")
    plot(
      xg, yg,
      type = "l",
      lwd = 1.7,
      col = "black",
      xlab = "x",
      ylab = expression(H^minute * (x)),
      ylim = ylim
    )
    app_draw_threshold_lines(xs, ylim)
    abline(h = -params$u, lty = 3, col = "#777777")
    abline(h = 0, lty = 3, col = "#BBBBBB")
    abline(h = params$l, lty = 4, col = "#777777")
    legend(
      "topleft",
      legend = expression(H^minute * (x), underline(x), x^kappa, x^lambda, bar(x), -c[U], c[D]),
      lty = c(1, style$ltys, 3, 4),
      col = c("black", style$cols, "#777777", "#777777"),
      lwd = c(1.7, rep(1.8, 4), 1.2, 1.2),
      bty = "n",
      cex = 0.95,
      ncol = 1
    )
    grid(col = "gray90")
    box()
  })
}

app_plot_reflected_with_controls <- function(sol, params, sim) {
  t <- sim$time
  xlim <- range(t, finite = TRUE)
  if (!all(is.finite(xlim)) || diff(xlim) <= 0) {
    xlim <- c(0, 1)
  }
  thresholds <- .thresholds_for_sim_plot(sol, sim)
  xL <- .as_plot_path(thresholds$xL, t, "xL")
  xU <- .as_plot_path(thresholds$xU, t, "xU")
  xk <- .as_plot_path(thresholds$xk, t, "xk")
  xl <- .as_plot_path(thresholds$xl, t, "xl")
  state_values <- c(sim$X, xL, xk, xl, xU)
  state_data_range <- range(state_values, finite = TRUE)
  if (!all(is.finite(state_data_range)) || diff(state_data_range) <= 0) {
    state_data_range <- app_expand_range(state_values, frac = 0.08)
  }
  state_top_blank <- 0.075
  state_bottom_blank <- 0.05
  state_data_share <- 1 - state_top_blank - state_bottom_blank
  if (!is.finite(state_data_share) || state_data_share <= 0) {
    stop("State-panel blank shares must sum to less than 1.", call. = FALSE)
  }
  state_total_height <- diff(state_data_range) / state_data_share
  y_state <- c(
    state_data_range[1] - state_bottom_blank * state_total_height,
    state_data_range[2] + state_top_blank * state_total_height
  )
  y_D <- app_expand_range(sim$L, frac = 0.05, fallback = c(0, 1))
  y_U <- app_expand_range(sim$U, frac = 0.05, fallback = c(0, 1))
  style <- app_threshold_style()

  app_safe_plot_par({
    layout(matrix(seq_len(3), ncol = 1), heights = c(0.78, 1.55, 0.78))
    on.exit(layout(1), add = TRUE)
    par(cex = 1, xaxs = "i")

    par(mar = c(0.8, 3.2, 0.8, 0.8), mgp = c(2.1, 0.65, 0),
        tcl = -0.25, cex.lab = 1.05, cex.axis = 1.0, xaxs = "i")
    plot(t, sim$L, type = "s", lwd = 1.8, col = style$cols[["xU"]],
         xlab = "", ylab = expression(D[t]), xaxt = "n", xlim = xlim, ylim = y_D)
    grid(col = "gray90")
    box()

    par(mar = c(0.8, 3.2, 0.2, 0.8), mgp = c(2.1, 0.65, 0),
        tcl = -0.25, cex.lab = 1.05, cex.axis = 1.0, xaxs = "i")
    plot(t, sim$X, type = "l", lwd = 1.2, col = "black",
         xlab = "", ylab = expression(X[t]), xaxt = "n", xlim = xlim, ylim = y_state)
    polygon(c(t, rev(t)), c(xL, rev(xU)),
            col = adjustcolor("gray85", 0.5), border = NA)
    lines(t, sim$X, lwd = 1.2, col = "black")
    lines(t, xL, col = style$cols[["xL"]], lty = style$ltys[["xL"]], lwd = 1.8)
    lines(t, xk, col = style$cols[["xk"]], lty = style$ltys[["xk"]], lwd = 1.8)
    lines(t, xl, col = style$cols[["xl"]], lty = style$ltys[["xl"]], lwd = 1.8)
    lines(t, xU, col = style$cols[["xU"]], lty = style$ltys[["xU"]], lwd = 1.8)
    legend(
      "topleft",
      legend = expression(X[t], underline(x), x^kappa, x^lambda, bar(x)),
      col = c("black", style$cols),
      lty = c(1, style$ltys),
      lwd = c(1.2, rep(1.8, 4)),
      bty = "n",
      cex = 0.95,
      ncol = 5
    )
    grid(col = "gray90")
    box()

    par(mar = c(3.0, 3.2, 0.2, 0.8), mgp = c(2.1, 0.65, 0),
        tcl = -0.25, cex.lab = 1.05, cex.axis = 1.0, xaxs = "i")
    plot(t, sim$U, type = "s", lwd = 1.8, col = style$cols[["xL"]],
         xlab = "time", ylab = expression(U[t]), xlim = xlim, ylim = y_U)
    grid(col = "gray90")
    box()
  })
}

app_plot_cost <- function(cost_df, gamma = NA_real_) {
  app_safe_plot_par({
    par(mar = c(3.0, 3.3, 0.8, 0.8), mgp = c(2.1, 0.65, 0), tcl = -0.25)

    ylim <- app_expand_range(c(cost_df$avg, gamma))
    plot(
      cost_df$time,
      cost_df$avg,
      type = "l",
      lwd = 2,
      col = "#F28E2B",
      xlab = "time",
      ylab = expression(J[t]),
      ylim = ylim
    )
    if (is.finite(gamma)) {
      abline(h = gamma, lty = 3, lwd = 1.5, col = "gray35")
      legend(
        "topright",
        legend = c(expression(J[t]), expression(gamma)),
        col = c("#F28E2B", "gray35"),
        lty = c(1, 3),
        lwd = c(2, 1.5),
        bty = "n"
      )
    }
    grid(col = "gray88")
    box()
  })
}

app_plot_solver_history <- function(sol) {
  hist <- sol$iter_history
  if (is.null(hist) || !is.data.frame(hist) || !nrow(hist)) {
    plot.new()
    text(0.5, 0.5, "Enable recorded iterates and solve again.")
    return(invisible(NULL))
  }

  app_safe_plot_par({
    style <- app_threshold_style()
    par(mfrow = c(1, 2))

    par(mar = c(3.1, 3.4, 0.8, 0.8), mgp = c(2.1, 0.65, 0), tcl = -0.25)
    plot(
      hist$step,
      hist$fmax,
      type = "b",
      log = "y",
      pch = 16,
      lwd = 1.5,
      xlab = "iteration",
      ylab = "max residual"
    )
    grid(col = "gray88")
    box()

    par(mar = c(3.1, 3.4, 0.8, 0.8), mgp = c(2.1, 0.65, 0), tcl = -0.25)
    matplot(
      hist$step,
      as.matrix(hist[, c("xL", "xk", "xl", "xU")]),
      type = "l",
      lty = style$ltys,
      lwd = 2,
      col = style$cols,
      xlab = "iteration",
      ylab = "threshold"
    )
    legend(
      "right",
      legend = expression(underline(x), x^kappa, x^lambda, bar(x)),
      lty = style$ltys,
      lwd = 2,
      col = style$cols,
      bty = "n"
    )
    grid(col = "gray88")
    box()
  })
}

app_plot_sweep_thresholds <- function(sweep_obj, include_gamma = TRUE) {
  plot_sweep(
    sweep_obj,
    show = TRUE,
    save = FALSE,
    title = FALSE,
    plot_gamma = isTRUE(include_gamma),
    gamma_layout = "stacked",
    tick_cex = 1,
    axis_title_cex = 1,
    margins = c(3.0, 3.2, 0.9, 0.8),
    axis_mgp = c(2.1, 0.65, 0),
    legend_cex = 0.9
  )
}
