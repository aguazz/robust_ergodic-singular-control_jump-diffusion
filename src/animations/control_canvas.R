# Control/state/cost animation canvas frames.
#
# Layout:
#   Left column:  D_t (top), X_t (middle), U_t (bottom).
#   Right column: J_t (finite-horizon running-average cost), spanning the lower
#   two rows. The top-right panel is a compact legend for Beamer slides.
#
# Loading this file defines functions only.

.draw_step_panel <- function(t, y, i, ylim, col, ylab_expr,
                             show_axis_labels,
                             bty, box_lwd,
                             show_x_axes = FALSE, show_y_axes = TRUE,
                             mar = c(0.4, 3.05, 0.4, 0.3),
                             tick_cex = 1.5,
                             axis_title_cex = NULL,
                             axis_mgp = NULL,
                             tcl = NULL,
                             base_cex = NULL,
                             mex = NULL) {
  yp <- .mask_after(y, i)
  .animation_apply_par(
    margins = mar,
    axis_mgp = axis_mgp,
    tcl = tcl,
    tick_cex = tick_cex,
    axis_title_cex = axis_title_cex,
    base_cex = base_cex,
    mex = mex
  )
  
  xlab0 <- if (isTRUE(show_axis_labels) && isTRUE(show_x_axes)) "time" else ""
  ylab0 <- if (isTRUE(show_axis_labels) && isTRUE(show_y_axes)) ylab_expr else ""
  
  # Match the reflected-state panel width by reserving the same right strip.
  xlim0 <- range(t)
  pad_x <- 0.12 * diff(xlim0)
  xlim <- c(xlim0[1], xlim0[2] + pad_x)
  
  plot(t, yp, type = "s", xlab = xlab0, ylab = ylab0, col = col,
       ylim = ylim, xlim = xlim, xaxt = "n", yaxt = "n", lwd = 2, bty = "n")
  
  if (isTRUE(show_x_axes)) {
    .animation_axis_T(1, t, tick_cex = tick_cex)
  }
  if (isTRUE(show_y_axes)) {
    axis(2, at = c(0, ylim[2]), labels = c("0", ""), cex.axis = tick_cex)
  }
  if (!identical(bty, "n")) box(bty = bty, lwd = box_lwd)
}

.draw_reflected_panel_frame <- function(sol, sim, i,
                                        show_ambiguity, show_barriers,
                                        show_axis_labels,
                                        show_x_axes, show_y_axes,
                                        bty, box_lwd,
                                        top_blank = 0.075, bottom_blank = 0.05,
                                        draw_legend = FALSE,
                                        mar = c(0.4, 3.05, 0.4, 0.3),
                                        tick_cex = 1.5,
                                        axis_title_cex = NULL,
                                        cost_scale_cex = 1.5,
                                        marker_cex = 1.5,
                                        state_point_cex = 1,
                                        legend_cex = NULL,
                                        axis_mgp = NULL,
                                        tcl = NULL,
                                        base_cex = NULL,
                                        mex = NULL) {
  t <- sim$time
  Xp <- .mask_after(sim$X, i)
  
  thresholds <- .thresholds_for_sim_plot(sol, sim)
  xL_vec <- .as_plot_path(thresholds$xL, t, "xL")
  xU_vec <- .as_plot_path(thresholds$xU, t, "xU")
  xk_vec <- .as_plot_path(thresholds$xk, t, "xk")
  xl_vec <- .as_plot_path(thresholds$xl, t, "xl")
  xLmin <- min(xL_vec, finite = TRUE)
  xUmax <- max(xU_vec, finite = TRUE)
  
  Hcorr <- xUmax - xLmin
  frac <- 1 - top_blank - bottom_blank
  Htot <- Hcorr / frac
  ylim <- c(xLmin - Htot * bottom_blank, xUmax + Htot * top_blank)
  
  .animation_apply_par(
    margins = mar,
    axis_mgp = axis_mgp,
    tcl = tcl,
    tick_cex = tick_cex,
    axis_title_cex = axis_title_cex,
    base_cex = base_cex,
    mex = mex
  )
  
  xlim0 <- range(t)
  pad_x <- 0.12 * diff(xlim0)
  xlim <- c(xlim0[1], xlim0[2] + pad_x)
  
  plot(t, Xp, type = "n", xlab = "", ylab = "",
       ylim = ylim, xlim = xlim, xaxt = "n", yaxt = "n", bty = "n")
  
  # Quadratic running-cost shading, clipped to the time-varying reflecting band.
  ybreaks <- seq(xLmin, xUmax, length.out = 101)
  ymid <- (ybreaks[-1] + ybreaks[-length(ybreaks)]) / 2
  cmin <- if (xLmin <= 0 && xUmax >= 0) 0 else min(xLmin^2, xUmax^2)
  cmax <- max(xLmin^2, xUmax^2)
  z <- (ymid^2 - cmin) / (cmax - cmin + 1e-12)
  fill <- adjustcolor(gray(1 - 0.75 * z), alpha.f = 0.45)
  pal <- adjustcolor(gray(1 - 0.75 * seq(0, 1, length.out = 256)), alpha.f = 0.45)
  Z <- matrix(rep(z, each = length(t)), nrow = length(t), ncol = length(ymid))
  Z[outer(xL_vec, ymid, `>`)] <- NA
  Z[outer(xU_vec, ymid, `<`)] <- NA
  image(t, ymid, Z, col = pal, add = TRUE, useRaster = TRUE)
  
  colors <- .animation_colors()
  cols4 <- c(
    colors[["lower_barrier"]],
    colors[["drift_threshold"]],
    colors[["intensity_threshold"]],
    colors[["upper_barrier"]]
  )
  ltys <- c(2, 2, 2, 2)
  
  if (isTRUE(show_barriers)) {
    lines(t, xL_vec, col = cols4[1], lty = ltys[1], lwd = 3)
    lines(t, xU_vec, col = cols4[4], lty = ltys[4], lwd = 3)
  }
  if (isTRUE(show_ambiguity)) {
    lines(t, xk_vec, col = cols4[2], lty = ltys[2], lwd = 2)
    lines(t, xl_vec, col = cols4[3], lty = ltys[3], lwd = 2)
  }
  
  # Running-cost color scale in the reserved strip.
  tmax <- xlim0[2]
  x0 <- tmax + 0.35 * pad_x
  x1 <- tmax + 0.55 * pad_x
  rect(xleft = x0, xright = x1,
       ybottom = ybreaks[-length(ybreaks)], ytop = ybreaks[-1],
       col = fill, border = NA)
  rect(x0, xLmin, x1, xUmax, border = "gray30", lwd = 0.8)
  points((x0 + x1) / 2, Xp[i], pch = 16, col = "black", cex = marker_cex, xpd = NA)
  
  ct <- pretty(c(cmin, cmax), n = 4)
  ct <- ct[ct >= cmin & ct <= cmax]
  sgn <- if (xUmax >= 0) 1 else -1
  yt <- sgn * sqrt(ct)
  ok <- yt >= xLmin & yt <= xUmax
  ct <- ct[ok]
  yt <- yt[ok]
  
  segments(x1, yt, x1 + 0.10 * pad_x, yt, col = "gray30", xpd = NA)
  text(x1 + 0.12 * pad_x, yt, labels = formatC(ct, format = "g", digits = 3),
       adj = c(0, 0.5), cex = cost_scale_cex, xpd = NA)
  text((x0 + x1) / 2, xUmax, labels = expression(c(x) == x^2),
       pos = 3, cex = cost_scale_cex, xpd = NA)
  
  lines(t, Xp, lwd = 1.2, col = "black")
  points(t[i], Xp[i], pch = 16, col = "black", cex = state_point_cex)
  
  if (isTRUE(show_x_axes)) {
    .animation_axis_T(1, t, tick_cex = tick_cex)
  }
  if (isTRUE(show_y_axes)) {
    at <- c(ylim[1], 0, ylim[2])
    labs <- c("", "0", "")
    keep <- at >= ylim[1] & at <= ylim[2]
    axis(2, at = at[keep], labels = labs[keep], cex.axis = tick_cex)
  }
  
  if (isTRUE(draw_legend)) {
    leg <- list(legend = c(expression(bar(X)[t])), lty = 1, col = "black", lwd = 1.2)
    if (isTRUE(show_barriers)) {
      leg$legend <- c(leg$legend, expression(underline(x), bar(x)))
      leg$lty <- c(leg$lty, 1, 1)
      leg$col <- c(leg$col, cols4[1], cols4[4])
      leg$lwd <- c(leg$lwd, 2, 2)
    }
    if (isTRUE(show_ambiguity)) {
      leg$legend <- c(leg$legend, expression(x^kappa, x^lambda))
      leg$lty <- c(leg$lty, 2, 2)
      leg$col <- c(leg$col, cols4[2], cols4[3])
      leg$lwd <- c(leg$lwd, 2, 2)
    }
    legend_args <- list(
      x = "topleft",
      bty = "n",
      horiz = TRUE,
      x.intersp = 0.5,
      seg.len = 2,
      legend = leg$legend,
      lty = leg$lty,
      col = leg$col,
      lwd = leg$lwd
    )
    if (!is.null(legend_cex)) legend_args$cex <- legend_cex
    do.call(legend, legend_args)
  }
  
  if (!identical(bty, "n")) box(bty = bty, lwd = box_lwd)
}

.draw_cost_panel <- function(cost_df, i, gamma, show_gamma,
                             show_axis_labels,
                             show_x_axes, show_y_axes,
                             bty, box_lwd,
                             mar = c(2.6, 3.05, 0.6, 0.8),
                             tick_cex = 1.5,
                             axis_title_cex = NULL,
                             axis_mgp = NULL,
                             tcl = NULL,
                             base_cex = NULL,
                             mex = NULL) {
  t <- cost_df$time
  Jp <- .mask_after(cost_df$avg, i)
  
  ylim <- range(cost_df$avg, finite = TRUE)
  if (!all(is.finite(ylim))) ylim <- c(0, 1)
  if (isTRUE(show_gamma) && is.finite(gamma)) ylim <- range(ylim, gamma)
  
  .animation_apply_par(
    margins = mar,
    axis_mgp = axis_mgp,
    tcl = tcl,
    tick_cex = tick_cex,
    axis_title_cex = axis_title_cex,
    base_cex = base_cex,
    mex = mex
  )
  
  xlab0 <- if (isTRUE(show_axis_labels)) "time" else ""
  ylab0 <- if (isTRUE(show_axis_labels)) expression(J[t]) else ""
  
  colors <- .animation_colors()
  plot(t, Jp, type = "l", xlab = xlab0, ylab = ylab0, col = colors[["cost"]],
       ylim = ylim, xaxt = "n", yaxt = "n", lwd = 2, bty = "n")
  
  if (isTRUE(show_x_axes)) {
    .animation_axis_T(1, t, tick_cex = tick_cex)
  }
  if (isTRUE(show_y_axes)) {
    axis(2, at = c(0, ylim[2]), labels = c("0", ""), cex.axis = tick_cex)
  }
  
  if (isTRUE(show_gamma) && is.finite(gamma)) abline(h = gamma, lty = 3)
  if (!identical(bty, "n")) box(bty = bty, lwd = box_lwd)
}

.draw_legend_panel <- function(show_ambiguity = TRUE,
                               show_barriers  = TRUE,
                               show_gamma = TRUE,
                               bty = "n", box_lwd = 1,
                               mar = c(0.4, 0.4, 0.4, 0.4),
                               legend_cex = NULL,
                               legend_min_cex = 1,
                               legend_ncol = NULL,
                               legend_y_intersp = 1.05,
                               base_cex = NULL,
                               mex = NULL) {
  .animation_apply_par(margins = mar, base_cex = base_cex, mex = mex)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  
  colors <- .animation_colors()
  base_legend <- data.frame(
    label = c(
      "State process",
      "Cumulative downwards interventions",
      "Cumulative upwards interventions",
      "Running-average cost"
    ),
    col = c(
      colors[["state"]],
      colors[["downward_control"]],
      colors[["upward_control"]],
      colors[["cost"]]
    ),
    lty = c(1, 1, 1, 1),
    lwd = c(2, 2, 2, 2),
    stringsAsFactors = FALSE
  )
  barrier_legend <- data.frame(
    label = c("Push-down barrier", "Push-up barrier"),
    col = c(colors[["lower_barrier"]], colors[["upper_barrier"]]),
    lty = c(2, 2),
    lwd = c(2, 2),
    stringsAsFactors = FALSE
  )
  ambiguity_legend <- data.frame(
    label = c("Drift ambiguity threshold", "Intensity ambiguity threshold"),
    col = c(colors[["drift_threshold"]], colors[["intensity_threshold"]]),
    lty = c(2, 2),
    lwd = c(2, 2),
    stringsAsFactors = FALSE
  )
  gamma_legend <- data.frame(
    label = "Ergodic value",
    col = "black",
    lty = 3,
    lwd = 1,
    stringsAsFactors = FALSE
  )
  
  legend_df <- base_legend
  if (isTRUE(show_barriers)) {
    legend_df <- rbind(legend_df, barrier_legend)
  }
  if (isTRUE(show_ambiguity)) {
    legend_df <- rbind(legend_df, ambiguity_legend)
  }
  if (isTRUE(show_gamma)) {
    legend_df <- rbind(legend_df, gamma_legend)
  }
  full_legend_df <- rbind(base_legend, barrier_legend, ambiguity_legend, gamma_legend)
  
  if (is.null(legend_ncol)) {
    legend_ncol <- 1L
  }
  legend_ncol <- as.integer(legend_ncol)[1]
  if (!is.finite(legend_ncol) || legend_ncol < 1L) {
    stop("legend_ncol must be a positive integer.")
  }
  
  cex <- if (is.null(legend_cex)) 1.5 else as.numeric(legend_cex)[1]
  if (!is.finite(cex) || cex <= 0) {
    stop("legend_cex must be NULL or a positive number.")
  }
  legend_min_cex <- as.numeric(legend_min_cex)[1]
  if (!is.finite(legend_min_cex) || legend_min_cex <= 0) {
    stop("legend_min_cex must be a positive number.")
  }
  
  longest_label <- "Cumulative downwards interventions"
  legend_args <- list(
    x = "left",
    bty = "n",
    x.intersp = 0.6,
    y.intersp = legend_y_intersp,
    seg.len = 2,
    ncol = legend_ncol
  )
  
  fit_args <- c(
    legend_args,
    list(
      legend = full_legend_df$label,
      col = full_legend_df$col,
      lty = full_legend_df$lty,
      lwd = full_legend_df$lwd
    )
  )
  
  draw_args <- c(
    legend_args,
    list(
      legend = legend_df$label,
      col = legend_df$col,
      lty = legend_df$lty,
      lwd = legend_df$lwd
    )
  )
  
  # Size the legend against the full robust entry set. That keeps the visual
  # scale consistent between the robust and non-robust canvas animations.
  repeat {
    fit <- do.call(
      legend,
      c(
        fit_args,
        list(
          cex = cex,
          text.width = strwidth(longest_label, cex = cex),
          plot = FALSE
        )
      )
    )
    fits_width <- is.null(fit$rect$w) || fit$rect$w <= 0.96
    fits_height <- is.null(fit$rect$h) || fit$rect$h <= 0.92
    if ((fits_width && fits_height) || cex <= legend_min_cex) break
    cex <- max(legend_min_cex, cex * 0.95)
  }
  
  do.call(
    legend,
    c(
      draw_args,
      list(
        cex = cex,
        text.width = strwidth(longest_label, cex = cex)
      )
    )
  )
  
  if (!identical(bty, "n")) box(bty = bty, lwd = box_lwd)
}

plot_frame_canvas_3left_1right <- function(
    sol, params, sim, i = length(sim$time),
    cost_df = NULL,
    out_dir = "frames/canvas", base = "canvas", digits = 5,
    width_in = 10.5, height_in = 5.5, dpi = 200,
    show_ambiguity = TRUE, show_barriers = TRUE, show_gamma = TRUE,
    show_axis_labels = TRUE, show_x_axes = TRUE, show_y_axes = TRUE,
    bottom_x_axis = TRUE, cost_x_axis = TRUE,
    bty = "o", box_lwd = 1,
    widths = c(0.58, 0.42),
    heights = c(0.7, 1.7, 0.7),
    mar_D = c(0.2, 3.05, 0.4, 0.3),
    mar_X = c(0.2, 3.05, 0.4, 0.3),
    mar_U = c(2.6, 3.05, 0.4, 0.3),
    mar_J = c(2.6, 3.05, 0.6, 0.8),
    draw_legend_reflected = FALSE,
    top_blank = 0.075, bottom_blank = 0.05,
    tick_cex = 1.5,
    axis_title_cex = NULL,
    cost_scale_cex = 1.5,
    legend_cex = NULL,
    legend_min_cex = 1,
    legend_ncol = NULL,
    marker_cex = 1.5,
    state_point_cex = 1,
    axis_mgp = NULL,
    tcl = NULL,
    base_cex = NULL,
    mex = NULL,
    show = FALSE, save = TRUE
) {
  stopifnot(i >= 1L, i <= length(sim$time))
  if (is.null(cost_df)) cost_df <- finite_horizon_ergodic_cost(sim, params)
  fname <- .frame_name(base, digits, i)
  
  ylim_D <- .animation_control_ylim(sim$L)
  ylim_U <- .animation_control_ylim(sim$U)
  colors <- .animation_colors()
  
  plotfun <- function() {
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    layout(
      matrix(c(1, 4,
               2, 5,
               3, 5), nrow = 3, byrow = TRUE),
      widths = widths,
      heights = heights
    )
    
    .draw_step_panel(
      sim$time, sim$L, i,
      ylim = ylim_D,
      col = colors[["downward_control"]],
      ylab_expr = expression(D[t]),
      show_axis_labels = show_axis_labels,
      show_x_axes = show_x_axes,
      show_y_axes = show_y_axes,
      bty = bty,
      box_lwd = box_lwd,
      mar = mar_D,
      tick_cex = tick_cex,
      axis_title_cex = axis_title_cex,
      axis_mgp = axis_mgp,
      tcl = tcl,
      base_cex = base_cex,
      mex = mex
    )
    
    .draw_reflected_panel_frame(
      sol, sim, i,
      show_ambiguity = show_ambiguity,
      show_barriers = show_barriers,
      show_axis_labels = show_axis_labels,
      show_x_axes = show_x_axes,
      show_y_axes = show_y_axes,
      bty = bty,
      box_lwd = box_lwd,
      top_blank = top_blank,
      bottom_blank = bottom_blank,
      draw_legend = draw_legend_reflected,
      mar = mar_X,
      tick_cex = tick_cex,
      axis_title_cex = axis_title_cex,
      cost_scale_cex = cost_scale_cex,
      marker_cex = marker_cex,
      state_point_cex = state_point_cex,
      legend_cex = NULL,
      axis_mgp = axis_mgp,
      tcl = tcl,
      base_cex = base_cex,
      mex = mex
    )
    
    .draw_step_panel(
      sim$time, sim$U, i,
      ylim = ylim_U,
      col = colors[["upward_control"]],
      ylab_expr = expression(U[t]),
      show_axis_labels = show_axis_labels,
      show_x_axes = bottom_x_axis,
      show_y_axes = show_y_axes,
      bty = bty,
      box_lwd = box_lwd,
      mar = mar_U,
      tick_cex = tick_cex,
      axis_title_cex = axis_title_cex,
      axis_mgp = axis_mgp,
      tcl = tcl,
      base_cex = base_cex,
      mex = mex
    )
    
    .draw_legend_panel(
      bty = "n",
      show_ambiguity = show_ambiguity,
      show_barriers = show_barriers,
      show_gamma = show_gamma,
      legend_cex = legend_cex,
      legend_min_cex = legend_min_cex,
      legend_ncol = legend_ncol,
      base_cex = base_cex,
      mex = mex
    )
    
    .draw_cost_panel(
      cost_df, i,
      gamma = sol$gamma %||% NA_real_,
      show_gamma = show_gamma,
      show_axis_labels = show_axis_labels,
      show_x_axes = cost_x_axis,
      show_y_axes = show_y_axes,
      bty = bty,
      box_lwd = box_lwd,
      mar = mar_J,
      tick_cex = tick_cex,
      axis_title_cex = axis_title_cex,
      axis_mgp = axis_mgp,
      tcl = tcl,
      base_cex = base_cex,
      mex = mex
    )
  }
  
  render_and_save_png(fname, plotfun, dir = out_dir,
                      width_in = width_in, height_in = height_in, dpi = dpi,
                      show = show, save = save, match_current = FALSE)
}

save_animation_frames_canvas <- function(
    sol, params, sim,
    out_dir = "frames/canvas", prefix = "canvas",
    every = 1L, digits = 5,
    width_in = 10.5, height_in = 5.5, dpi = 200,
    show_ambiguity = TRUE, show_barriers = TRUE, show_gamma = TRUE,
    show_axis_labels = TRUE,
    show_x_axes = TRUE, show_y_axes = TRUE,
    bottom_x_axis = TRUE, cost_x_axis = TRUE,
    bty = "o", box_lwd = 1,
    widths = c(0.58, 0.42), heights = c(0.7, 1.7, 0.7),
    draw_legend_reflected = FALSE,
    top_blank = 0.075, bottom_blank = 0.05,
    tick_cex = 1.5,
    axis_title_cex = NULL,
    cost_scale_cex = 1.5,
    legend_cex = NULL,
    legend_min_cex = 1,
    legend_ncol = NULL,
    marker_cex = 1.5,
    state_point_cex = 1,
    axis_mgp = NULL,
    tcl = NULL,
    base_cex = NULL,
    mex = NULL,
    frame_pause = 0.5,
    show = FALSE, save = TRUE
) {
  cost_df <- finite_horizon_ergodic_cost(sim, params)
  n <- length(sim$time)
  
  for (i in seq.int(2L, n, by = as.integer(every))) {
    plot_frame_canvas_3left_1right(
      sol, params, sim, i,
      cost_df = cost_df,
      out_dir = out_dir,
      base = prefix,
      digits = digits,
      width_in = width_in,
      height_in = height_in,
      dpi = dpi,
      show_ambiguity = show_ambiguity,
      show_barriers = show_barriers,
      show_gamma = show_gamma,
      show_axis_labels = show_axis_labels,
      show_x_axes = show_x_axes,
      show_y_axes = show_y_axes,
      bottom_x_axis = bottom_x_axis,
      cost_x_axis = cost_x_axis,
      bty = bty,
      box_lwd = box_lwd,
      widths = widths,
      heights = heights,
      draw_legend_reflected = draw_legend_reflected,
      top_blank = top_blank,
      bottom_blank = bottom_blank,
      tick_cex = tick_cex,
      axis_title_cex = axis_title_cex,
      cost_scale_cex = cost_scale_cex,
      legend_cex = legend_cex,
      legend_min_cex = legend_min_cex,
      legend_ncol = legend_ncol,
      marker_cex = marker_cex,
      state_point_cex = state_point_cex,
      axis_mgp = axis_mgp,
      tcl = tcl,
      base_cex = base_cex,
      mex = mex,
      show = show,
      save = save
    )
    if (isTRUE(show) && is.finite(frame_pause) && frame_pause > 0) {
      Sys.sleep(frame_pause)
    }
  }
  invisible(list(cost = cost_df))
}
