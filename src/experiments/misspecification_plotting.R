# Misspecification-grid plotting helpers.
#
# Notation:
#   delta, eps
#     Horizontal grid coordinates in the misspecification experiments.
#   zmat
#     Matrix of values for the selected surface, indexed by delta rows and eps
#     columns.
#   RMC
#     Robust misspecification cost, plotted by default in percentage points.
#
# Loading this file defines functions only.

.validate_margin_spec <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 4L || any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("'%s' must be a numeric vector of length 4 with non-negative entries.", name))
  }
  x
}

.validate_fig_spec <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 4L || any(!is.finite(x))) {
    stop(sprintf("'%s' must be a numeric vector of length 4.", name))
  }
  if (any(x < 0 | x > 1)) {
    stop(sprintf("'%s' entries must lie between 0 and 1.", name))
  }
  if (x[1] >= x[2] || x[3] >= x[4]) {
    stop(sprintf("'%s' must satisfy left < right and bottom < top.", name))
  }
  x
}

.validate_range_spec <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 2L || any(!is.finite(x))) {
    stop(sprintf("'%s' must be a numeric vector of length 2.", name))
  }
  x <- sort(x)
  if (x[1] == x[2]) {
    stop(sprintf("'%s' must have distinct endpoints.", name))
  }
  x
}

.misspecification_surface_plot_meta <- function(surface) {
  list(
    z_axis_title = switch(
      surface,
      relative_cost = expression(plain(RMC) / 100 == (gamma[NR*","*W] - gamma[R]) / gamma[R]),
      relative_cost_pct = expression(plain(RMC)),
      gamma_difference = expression(gamma[NR*","*W] - gamma[R]),
      gamma_with_ambiguity = expression(gamma[R]),
      gamma_zero_ambiguity_policy = expression(gamma[NR*","*W])
    ),
    color_scale_title = switch(
      surface,
      relative_cost = expression(plain(RMC) / 100),
      relative_cost_pct = expression(plain(RMC)),
      gamma_difference = expression(gamma[NR*","*W] - gamma[R]),
      gamma_with_ambiguity = expression(gamma[R]),
      gamma_zero_ambiguity_policy = expression(gamma[NR*","*W])
    ),
    main = switch(
      surface,
      relative_cost = "RMC / 100",
      relative_cost_pct = "RMC",
      gamma_difference = "Excess ergodic cost of the non-robust policy",
      gamma_with_ambiguity = "Robust ergodic cost",
      gamma_zero_ambiguity_policy = "Non-robust policy under the worst-case model"
    )
  )
}

.validate_unit_interval_scalar <- function(x, name) {
  x <- as.numeric(x)[1]
  if (!is.finite(x) || x < 0 || x > 1) {
    stop(sprintf("'%s' must be a single number between 0 and 1.", name))
  }
  x
}

.validate_nonnegative_scalar <- function(x, name) {
  x <- as.numeric(x)[1]
  if (!is.finite(x) || x < 0) {
    stop(sprintf("'%s' must be a single non-negative number.", name))
  }
  x
}

.validate_xy_nudge <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 2L || any(!is.finite(x))) {
    stop(sprintf("'%s' must be a numeric vector of length 2.", name))
  }
  x
}

.coerce_plot_label <- function(label) {
  if (is.null(label) || !length(label)) return("")
  if (is.language(label) && !is.expression(label)) {
    return(as.expression(label))
  }
  label
}

.draw_parametrizer_legend_key <- function(
    label,
    corner = c("topleft", "bottomleft"),
    inset = c(0.04, 0.06),
    line_length = 0.18,
    cex = 1,
    col = "#333333",
    lty = 1,
    lwd = 1.5
) {
  label <- .coerce_plot_label(label)
  corner <- match.arg(corner)
  inset <- as.numeric(inset)
  if (length(inset) != 2L || any(!is.finite(inset)) ||
      any(inset < 0 | inset > 1)) {
    stop("'inset' must be a numeric vector of length 2 with entries in [0, 1].")
  }
  line_length <- .validate_unit_interval_scalar(line_length, "line_length")
  cex <- .validate_nonnegative_scalar(cex, "cex")
  
  usr <- par("usr")
  xspan <- diff(usr[1:2])
  yspan <- diff(usr[3:4])
  if (!is.finite(xspan) || !is.finite(yspan) || xspan <= 0 || yspan <= 0) {
    return(invisible(NULL))
  }
  
  x0 <- usr[1] + inset[1] * xspan
  x1 <- min(usr[2], x0 + line_length * xspan)
  y0 <- if (identical(corner, "topleft")) {
    usr[4] - inset[2] * yspan
  } else {
    usr[3] + inset[2] * yspan
  }
  
  total_width <- x1 - x0
  if (!is.finite(total_width) || total_width <= 0) {
    return(invisible(NULL))
  }
  
  gap_width <- max(
    1.25 * strwidth(label, cex = cex, units = "user"),
    0.025 * xspan
  )
  gap_width <- min(gap_width, 0.8 * total_width)
  seg_width <- 0.5 * (total_width - gap_width)
  
  if (seg_width > 0) {
    segments(x0, y0, x0 + seg_width, y0, col = col, lty = lty, lwd = lwd)
    segments(x1 - seg_width, y0, x1, y0, col = col, lty = lty, lwd = lwd)
  }
  
  text(x0 + total_width/2, y0, labels = label, cex = cex, col = col)
  
  invisible(NULL)
}

.draw_curve_with_gap <- function(x, y, col, lty, lwd, x_gap_left, x_gap_right) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2L) {
    return(invisible(NULL))
  }
  
  gap_invalid <- !is.finite(x_gap_left) || !is.finite(x_gap_right) ||
    x_gap_left >= x_gap_right ||
    x_gap_right <= min(x) ||
    x_gap_left >= max(x)
  if (gap_invalid) {
    lines(x, y, col = col, lty = lty, lwd = lwd)
    return(invisible(NULL))
  }
  
  x_gap_left <- max(x_gap_left, min(x))
  x_gap_right <- min(x_gap_right, max(x))
  if (x_gap_left >= x_gap_right) {
    lines(x, y, col = col, lty = lty, lwd = lwd)
    return(invisible(NULL))
  }
  
  y_gap_left <- approx(x, y, xout = x_gap_left, ties = mean, rule = 1)$y
  y_gap_right <- approx(x, y, xout = x_gap_right, ties = mean, rule = 1)$y
  
  x_left <- x[x < x_gap_left]
  y_left <- y[x < x_gap_left]
  if (is.finite(y_gap_left) && x_gap_left > min(x)) {
    x_left <- c(x_left, x_gap_left)
    y_left <- c(y_left, y_gap_left)
  }
  
  x_right <- x[x > x_gap_right]
  y_right <- y[x > x_gap_right]
  if (is.finite(y_gap_right) && x_gap_right < max(x)) {
    x_right <- c(x_gap_right, x_right)
    y_right <- c(y_gap_right, y_right)
  }
  
  if (length(x_left) >= 2L) {
    lines(x_left, y_left, col = col, lty = lty, lwd = lwd)
  }
  if (length(x_right) >= 2L) {
    lines(x_right, y_right, col = col, lty = lty, lwd = lwd)
  }
  
  invisible(NULL)
}

.curve_label_info <- function(x, y, target_x, gap_width) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2L) return(NULL)
  if (target_x < min(x) || target_x > max(x)) return(NULL)
  
  y_target <- approx(x, y, xout = target_x, ties = mean, rule = 1)$y
  if (!is.finite(y_target)) return(NULL)
  
  left_x <- max(min(x), target_x - gap_width / 2)
  right_x <- min(max(x), target_x + gap_width / 2)
  if (left_x >= right_x) return(NULL)
  
  list(x = target_x, y = y_target, left_x = left_x, right_x = right_x)
}

.persp_axis_edge_candidates <- function(axis, xlim, ylim, zlim) {
  axis <- match.arg(axis, c("x", "y", "z"))
  x0 <- xlim[1]; x1 <- xlim[2]
  y0 <- ylim[1]; y1 <- ylim[2]
  z0 <- zlim[1]; z1 <- zlim[2]
  
  switch(
    axis,
    x = list(
      rbind(c(x0, y0, z0), c(x1, y0, z0)),
      rbind(c(x0, y1, z0), c(x1, y1, z0)),
      rbind(c(x0, y0, z1), c(x1, y0, z1)),
      rbind(c(x0, y1, z1), c(x1, y1, z1))
    ),
    y = list(
      rbind(c(x0, y0, z0), c(x0, y1, z0)),
      rbind(c(x1, y0, z0), c(x1, y1, z0)),
      rbind(c(x0, y0, z1), c(x0, y1, z1)),
      rbind(c(x1, y0, z1), c(x1, y1, z1))
    ),
    z = list(
      rbind(c(x0, y0, z0), c(x0, y0, z1)),
      rbind(c(x1, y0, z0), c(x1, y0, z1)),
      rbind(c(x0, y1, z0), c(x0, y1, z1)),
      rbind(c(x1, y1, z0), c(x1, y1, z1))
    )
  )
}

.resolve_persp_axis_edge <- function(axis, pmat, xlim, ylim, zlim) {
  candidates <- .persp_axis_edge_candidates(axis, xlim, ylim, zlim)
  center_proj <- trans3d(mean(xlim), mean(ylim), mean(zlim), pmat = pmat)
  center_xy <- c(center_proj$x, center_proj$y)
  
  mids <- do.call(
    rbind,
    lapply(candidates, function(edge) {
      proj <- trans3d(edge[, 1], edge[, 2], edge[, 3], pmat = pmat)
      c(mean(proj$x), mean(proj$y))
    })
  )
  
  pick <- switch(
    axis,
    x = order(mids[, 2], -abs(mids[, 1] - center_xy[1]))[1],
    y = order(mids[, 2], -abs(mids[, 1] - center_xy[1]))[1],
    z = order(-abs(mids[, 1] - center_xy[1]), mids[, 2])[1]
  )
  
  candidates[[pick]]
}

.persp_axis_geometry <- function(axis, pmat, xlim, ylim, zlim) {
  edge <- .resolve_persp_axis_edge(axis, pmat = pmat, xlim = xlim, ylim = ylim, zlim = zlim)
  proj <- trans3d(edge[, 1], edge[, 2], edge[, 3], pmat = pmat)
  p0 <- c(proj$x[1], proj$y[1])
  p1 <- c(proj$x[2], proj$y[2])
  v <- p1 - p0
  len <- sqrt(sum(v^2))
  if (!is.finite(len) || len <= 0) {
    return(NULL)
  }
  
  dir <- v / len
  normal <- c(-dir[2], dir[1])
  center_proj <- trans3d(mean(xlim), mean(ylim), mean(zlim), pmat = pmat)
  center_xy <- c(center_proj$x, center_proj$y)
  mid <- 0.5 * (p0 + p1)
  if (sum(normal * (mid - center_xy)) < 0) {
    normal <- -normal
  }
  
  angle <- atan2(dir[2], dir[1]) * 180 / pi
  if (angle > 90) angle <- angle - 180
  if (angle < -90) angle <- angle + 180
  
  list(
    edge = edge,
    p0 = p0,
    p1 = p1,
    v = v,
    dir = dir,
    normal = normal,
    angle = angle
  )
}

.validate_positive_scalar_or_null <- function(x, name) {
  if (is.null(x)) return(NULL)
  x <- as.numeric(x)[1]
  if (!is.finite(x) || x <= 0) {
    stop(sprintf("'%s' must be NULL or a single positive number.", name))
  }
  x
}

.map_axis_values_to_persp <- function(values, from, stretch = 1, scale = TRUE) {
  dims <- dim(values)
  dim_names <- dimnames(values)
  values <- as.numeric(values)
  from <- as.numeric(from)
  span <- diff(from)
  if (!is.finite(span) || span <= 0) {
    out <- rep(0.5 * stretch, length(values))
  } else if (isTRUE(scale)) {
    out <- (values - from[1]) / span * stretch
  } else {
    out <- (values - from[1]) * stretch
  }
  
  if (!is.null(dims)) {
    dim(out) <- dims
    dimnames(out) <- dim_names
  }
  out
}

.pretty_ticks_within <- function(lim, n = 5) {
  lim <- sort(as.numeric(lim))
  at <- pretty(lim, n = n)
  tol <- sqrt(.Machine$double.eps) * max(1, diff(lim))
  at <- at[is.finite(at) & at >= (lim[1] - tol) & at <= (lim[2] + tol)]
  if (!length(at)) {
    at <- lim
  }
  sort(unique(at))
}

.format_projected_axis_labels <- function(values) {
  values <- as.numeric(values)
  tol <- sqrt(.Machine$double.eps) * max(1, max(abs(values), na.rm = TRUE))
  values[abs(values) < tol] <- 0
  format(signif(values, 8), trim = TRUE, scientific = FALSE)
}

.draw_manual_persp_axis_ticks <- function(
    axis,
    pmat,
    xlim,
    ylim,
    zlim,
    axis_limits,
    at = NULL,
    n = 5,
    cex = 1,
    col = par("col.axis"),
    font = par("font.axis"),
    lwd = 1,
    tcl = par("tcl"),
    ticktype = c("detailed", "simple"),
    axis_mgp = c(2.2, 0.7, 0),
    draw_axis_line = TRUE
) {
  ticktype <- match.arg(ticktype)
  geom <- .persp_axis_geometry(axis, pmat = pmat, xlim = xlim, ylim = ylim, zlim = zlim)
  if (is.null(geom)) {
    return(invisible(NULL))
  }
  
  axis_limits <- sort(as.numeric(axis_limits))
  plot_limits <- switch(axis, x = xlim, y = ylim, z = zlim)
  if (is.null(at)) {
    at <- .pretty_ticks_within(axis_limits, n = n)
  } else {
    at <- sort(unique(as.numeric(at)))
    at <- at[is.finite(at)]
  }
  tol <- sqrt(.Machine$double.eps) * max(1, diff(axis_limits))
  at <- at[at >= (axis_limits[1] - tol) & at <= (axis_limits[2] + tol)]
  if (!length(at)) {
    return(invisible(NULL))
  }
  
  at_plot <- .map_axis_values_to_persp(
    values = at,
    from = axis_limits,
    stretch = diff(plot_limits),
    scale = TRUE
  ) + plot_limits[1]
  
  coords <- switch(
    axis,
    x = cbind(at_plot, rep(geom$edge[1, 2], length(at_plot)), rep(geom$edge[1, 3], length(at_plot))),
    y = cbind(rep(geom$edge[1, 1], length(at_plot)), at_plot, rep(geom$edge[1, 3], length(at_plot))),
    z = cbind(rep(geom$edge[1, 1], length(at_plot)), rep(geom$edge[1, 2], length(at_plot)), at_plot)
  )
  proj_ticks <- trans3d(coords[, 1], coords[, 2], coords[, 3], pmat = pmat)
  tick_xy <- cbind(proj_ticks$x, proj_ticks$y)
  
  char_size <- max(
    strwidth("0", units = "user", cex = cex),
    strheight("0", units = "user", cex = cex)
  )
  tick_scale <- abs(as.numeric(tcl)[1])
  if (!is.finite(tick_scale)) tick_scale <- 0.25
  tick_length <- if (tick_scale == 0) 0 else {
    (0.35 + if (identical(ticktype, "detailed")) 0.12 else 0) *
      char_size * (tick_scale / 0.25)
  }
  label_gap <- tick_length + (0.35 + max(axis_mgp[2], 0)) * char_size
  
  if (isTRUE(draw_axis_line)) {
    segments(geom$p0[1], geom$p0[2], geom$p1[1], geom$p1[2], col = col, lwd = lwd, xpd = NA)
  }
  if (tick_length > 0) {
    tick_end <- tick_xy + matrix(geom$normal, nrow = nrow(tick_xy), ncol = 2, byrow = TRUE) * tick_length
    segments(tick_xy[, 1], tick_xy[, 2], tick_end[, 1], tick_end[, 2], col = col, lwd = lwd, xpd = NA)
  }
  
  label_xy <- tick_xy + matrix(geom$normal, nrow = nrow(tick_xy), ncol = 2, byrow = TRUE) * label_gap
  text(
    label_xy[, 1],
    label_xy[, 2],
    labels = .format_projected_axis_labels(at),
    cex = cex,
    col = col,
    font = font,
    srt = geom$angle,
    adj = c(0.5, 0.5),
    xpd = NA
  )
  
  invisible(at)
}

.draw_manual_persp_axis_title <- function(
    label,
    axis,
    pmat,
    xlim,
    ylim,
    zlim,
    at = 0.5,
    pad = 0.08,
    nudge = c(0, 0),
    cex = 1,
    col = par("col.lab"),
    font = par("font.lab")
) {
  label <- .coerce_plot_label(label)
  if (is.character(label) && (!length(label) || !nzchar(label[1]))) {
    return(invisible(NULL))
  }
  
  geom <- .persp_axis_geometry(axis, pmat = pmat, xlim = xlim, ylim = ylim, zlim = zlim)
  if (is.null(geom)) {
    return(invisible(NULL))
  }
  
  usr <- par("usr")
  pad_units <- pad * max(diff(usr[1:2]), diff(usr[3:4]))
  nudge_units <- c(
    nudge[1] * diff(usr[1:2]),
    nudge[2] * diff(usr[3:4])
  )
  pos <- geom$p0 + at * geom$v + pad_units * geom$normal + nudge_units
  
  text(
    pos[1],
    pos[2],
    labels = label,
    srt = geom$angle,
    adj = c(0.5, 0.5),
    xpd = NA,
    cex = cex,
    col = col,
    font = font
  )
  
  invisible(pos)
}

.resolve_color_scale_fig <- function(
    color_scale_position,
    color_scale_fig,
    color_scale_outer_fig
) {
  if (identical(color_scale_position, "figure")) {
    return(color_scale_fig)
  }
  
  omd <- par("omd")
  plt <- par("plt")
  if (!all(is.finite(omd)) || !all(is.finite(plt))) {
    stop("Could not determine the current plotting regions for the color scale.")
  }
  
  x_bounds <- c(omd[2], 1)
  if ((x_bounds[2] - x_bounds[1]) <= 0) {
    stop("Set outer_margins[4] > 0 to use color_scale_position = 'right_outer'.")
  }
  
  y_bounds <- c(
    omd[3] + (omd[4] - omd[3]) * plt[3],
    omd[3] + (omd[4] - omd[3]) * plt[4]
  )
  if ((y_bounds[2] - y_bounds[1]) <= 0) {
    stop("The plot region is too small to place the color scale in the right outer margin.")
  }
  
  c(
    x_bounds[1] + (x_bounds[2] - x_bounds[1]) * color_scale_outer_fig[1],
    x_bounds[1] + (x_bounds[2] - x_bounds[1]) * color_scale_outer_fig[2],
    y_bounds[1] + (y_bounds[2] - y_bounds[1]) * color_scale_outer_fig[3],
    y_bounds[1] + (y_bounds[2] - y_bounds[1]) * color_scale_outer_fig[4]
  )
}

.draw_vertical_color_scale_panel <- function(
    zlim,
    color_palette,
    color_scale_title,
    color_scale_axis_mgp,
    color_scale_nticks,
    color_scale_title_cex,
    color_scale_tick_cex,
    color_scale_border
) {
  op_scale <- par(
    mgp = color_scale_axis_mgp,
    xpd = NA,
    cex.axis = color_scale_tick_cex
  )
  on.exit(par(op_scale), add = TRUE)
  
  plot.new()
  plot.window(xlim = c(0, 1), ylim = zlim, xaxs = "i", yaxs = "i")
  
  y_breaks <- seq(zlim[1], zlim[2], length.out = length(color_palette) + 1L)
  for (k in seq_along(color_palette)) {
    rect(0, y_breaks[k], 1, y_breaks[k + 1L], col = color_palette[k], border = NA)
  }
  box(col = color_scale_border)
  tick_at <- pretty(zlim, n = color_scale_nticks)
  tick_tol <- sqrt(.Machine$double.eps) * max(1, diff(zlim))
  tick_at <- tick_at[tick_at >= (zlim[1] - tick_tol) & tick_at <= (zlim[2] + tick_tol)]
  if (!length(tick_at)) {
    tick_at <- zlim
  }
  axis(4, at = tick_at, las = 1)
  mtext(color_scale_title, side = 3, line = 0.15, cex = color_scale_title_cex)
  
  invisible(NULL)
}

.draw_vertical_color_scale <- function(
    zlim,
    color_palette,
    color_scale_title,
    legend_fig,
    color_scale_margins,
    color_scale_axis_mgp,
    color_scale_nticks,
    color_scale_title_cex,
    color_scale_tick_cex,
    color_scale_border
) {
  op_legend <- par(
    fig = legend_fig,
    mar = color_scale_margins,
    new = TRUE
  )
  on.exit(par(op_legend), add = TRUE)
  
  .draw_vertical_color_scale_panel(
    zlim = zlim,
    color_palette = color_palette,
    color_scale_title = color_scale_title,
    color_scale_axis_mgp = color_scale_axis_mgp,
    color_scale_nticks = color_scale_nticks,
    color_scale_title_cex = color_scale_title_cex,
    color_scale_tick_cex = color_scale_tick_cex,
    color_scale_border = color_scale_border
  )
}

plot_misspecification_cost_3d <- function(
    grid_obj,                                                   # output of misspecification_cost_grid()
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"), # which matrix in grid_obj to plot as the surface height
    out_dir = file.path("figures", "misspecification"),         # directory where the figure is saved
    out_name = NULL,                                            # base filename for saved outputs
    show = interactive(),                                       # if TRUE, draw on the current device
    save = FALSE,                                               # if TRUE, save PNG and PDF copies
    width_in = 6.67,                                            # saved figure width in inches
    height_in = 4.67,                                           # saved figure height in inches
    dpi = 300,                                                  # PNG resolution used on save
    match_current = FALSE,                                      # if TRUE, save using the current device size instead of width_in/height_in
    title = TRUE,                                               # if TRUE, use the default plot title
    main = NULL,                                                # custom main title; overrides the default when supplied
    show_x_axis_title = TRUE,                                   # if TRUE, show the x-axis title
    show_y_axis_title = TRUE,                                   # if TRUE, show the y-axis title
    show_z_axis_title = TRUE,                                   # if TRUE, show the z-axis title
    x_axis_title = expression(delta),                           # x-axis title shown along the projected x-axis
    y_axis_title = expression(epsilon),                         # y-axis title shown along the projected y-axis
    z_axis_title = NULL,                                        # z-axis title; if NULL, use the default for the chosen surface
    axis_title_cex = 1,                                         # size of the axis titles
    x_axis_title_at = 0.5,                                      # position of the x-axis title along its projected axis (0 = start, 1 = end)
    y_axis_title_at = 0.5,                                      # position of the y-axis title along its projected axis (0 = start, 1 = end)
    z_axis_title_at = 0.5,                                      # position of the z-axis title along its projected axis (0 = start, 1 = end)
    x_axis_title_pad = 0.08,                                    # outward distance of the x-axis title from its projected axis
    y_axis_title_pad = 0.08,                                    # outward distance of the y-axis title from its projected axis
    z_axis_title_pad = 0.10,                                    # outward distance of the z-axis title from its projected axis
    tick_cex = 1.2,                                             # size of the axis tick labels
    title_cex = 1,                                              # size of the main title
    margins = c(3.2, 3.2, 1.8, 0.8),                           # margins around the 3D panel, in par(mar) order: bottom, left, top, right
    axis_mgp = c(2.2, 0.7, 0),                                 # spacing control used by the custom projected tick labels
    tcl = -0.25,                                                # tick length
    theta = 35,                                                 # horizontal viewing angle
    phi = 25,                                                   # vertical viewing angle
    r = sqrt(3),                                                # distance of the eyepoint from the center of the 3D box
    d = 2,                                                      # strength of the perspective effect
    scale = TRUE,                                               # if TRUE, normalize x, y, and z before applying the internal stretch; if FALSE, preserve their original relative scales
    auto_stretch = TRUE,                                        # if TRUE, compute default internal x/y stretch factors from the actual panel aspect ratio
    x_stretch = NULL,                                           # manual stretch applied to the internal delta direction; NULL uses the automatic value when available
    y_stretch = NULL,                                           # manual stretch applied to the internal epsilon direction; NULL uses the automatic value when available
    expand = 0.7,                                               # vertical expansion factor of the surface
    shade = 0.35,                                               # amount of lighting/shading on the surface
    ltheta = -135,                                              # horizontal angle of the light source
    lphi = 0,                                                   # vertical angle of the light source
    use_z_colormap = TRUE,                                      # if TRUE, color the surface according to z-values
    color_min = "#1E90FF",                                      # color used for the minimum z-values
    color_max = "#B22222",                                      # color used for the maximum z-values
    color_palette = NULL,                                       # custom vector of colors; overrides color_min/color_max when supplied
    color_steps = 64,                                           # number of colors in the generated palette
    add_color_scale = TRUE,                                     # if TRUE, draw a vertical color scale in a narrow panel to the right
    color_scale_width = 0.14,                                   # fraction of the total figure width reserved for the right-hand color-scale panel
    color_scale_title = NULL,                                   # title shown above the color scale
    color_scale_nticks = 5,                                     # approximate number of tick marks on the color scale
    color_scale_title_cex = 0.9,                                # size of the color-scale title
    color_scale_tick_cex = 0.85,                                # size of the color-scale tick labels
    color_scale_border = "#444444",                             # border color of the color scale
    col = "#7FA7A2",                                            # fallback constant surface color when use_z_colormap = FALSE
    border = NA,                                                # border color of the surface facets
    box = TRUE,                                                 # if TRUE, draw the 3D bounding box
    axes = TRUE,                                                # if TRUE, draw custom projected axes and tick labels using the true, un-stretched values
    nticks = 5,                                                 # approximate number of tick marks per custom projected axis
    ticktype = c("detailed", "simple"),                         # style used for the custom projected tick marks
    zlim = NULL,                                                # custom z-range; if NULL, compute it from the plotted surface
    restore_par = TRUE                                          # if FALSE, do not restore par() on exit; useful inside a shared layout
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  ticktype <- match.arg(ticktype)
  surface_meta <- .misspecification_surface_plot_meta(surface)
  
  margins <- .validate_margin_spec(margins, "margins")
  x_axis_title_at <- .validate_unit_interval_scalar(x_axis_title_at, "x_axis_title_at")
  y_axis_title_at <- .validate_unit_interval_scalar(y_axis_title_at, "y_axis_title_at")
  z_axis_title_at <- .validate_unit_interval_scalar(z_axis_title_at, "z_axis_title_at")
  x_axis_title_pad <- .validate_nonnegative_scalar(x_axis_title_pad, "x_axis_title_pad")
  y_axis_title_pad <- .validate_nonnegative_scalar(y_axis_title_pad, "y_axis_title_pad")
  z_axis_title_pad <- .validate_nonnegative_scalar(z_axis_title_pad, "z_axis_title_pad")
  x_stretch <- .validate_positive_scalar_or_null(x_stretch, "x_stretch")
  y_stretch <- .validate_positive_scalar_or_null(y_stretch, "y_stretch")
  color_scale_width <- .validate_unit_interval_scalar(color_scale_width, "color_scale_width")
  if (color_scale_width <= 0 || color_scale_width >= 1) {
    stop("'color_scale_width' must be strictly between 0 and 1.")
  }
  color_steps <- as.integer(color_steps)[1]
  if (!is.finite(color_steps) || color_steps < 2L) {
    stop("'color_steps' must be an integer >= 2.")
  }
  
  if (is.null(grid_obj$delta_values) || is.null(grid_obj$eps_values)) {
    stop("grid_obj must be the output of misspecification_cost_grid().")
  }
  
  delta_values <- as.numeric(grid_obj$delta_values)
  eps_values <- as.numeric(grid_obj$eps_values)
  zmat <- grid_obj[[surface]]
  
  zmat <- .ensure_misspecification_surface_has_finite_values(grid_obj, surface)
  if (length(delta_values) < 2L || length(eps_values) < 2L) {
    stop("A 3D surface plot needs at least two delta values and two eps values.")
  }
  
  delta_ord <- order(delta_values)
  eps_ord <- order(eps_values)
  delta_values <- delta_values[delta_ord]
  eps_values <- eps_values[eps_ord]
  zmat <- zmat[delta_ord, eps_ord, drop = FALSE]
  zmat[!is.finite(zmat)] <- NA_real_
  
  if (is.null(zlim)) {
    zlim <- range(zmat, finite = TRUE)
  }
  if (diff(zlim) == 0) {
    bump <- max(1, abs(zlim[1]))
    zlim <- zlim + c(-1, 1) * 0.05 * bump
  }
  
  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(c(color_min, color_max))(color_steps)
  } else {
    color_palette <- as.character(color_palette)
  }
  if (!length(color_palette)) {
    stop("color_palette must contain at least one color.")
  }
  
  draw_color_scale <- isTRUE(use_z_colormap) && isTRUE(add_color_scale)
  
  if (is.null(z_axis_title)) {
    z_axis_title <- surface_meta$z_axis_title
  }
  if (is.null(color_scale_title)) {
    color_scale_title <- surface_meta$color_scale_title
  }
  if (is.null(main)) {
    main <- if (isTRUE(title)) surface_meta$main else ""
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_3d")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  if (isTRUE(use_z_colormap)) {
    n_delta <- nrow(zmat)
    n_eps <- ncol(zmat)
    z11 <- zmat[-n_delta, -n_eps, drop = FALSE]
    z21 <- zmat[-1,      -n_eps, drop = FALSE]
    z12 <- zmat[-n_delta, -1,     drop = FALSE]
    z22 <- zmat[-1,      -1,      drop = FALSE]
    facet_ok <- is.finite(z11) & is.finite(z21) & is.finite(z12) & is.finite(z22)
    facet_z <- matrix(NA_real_, nrow = n_delta - 1L, ncol = n_eps - 1L)
    facet_z[facet_ok] <- (z11[facet_ok] + z21[facet_ok] + z12[facet_ok] + z22[facet_ok]) / 4
    
    breaks <- seq(zlim[1], zlim[2], length.out = length(color_palette) + 1L)
    facet_cols <- matrix(grDevices::adjustcolor("white", alpha.f = 0),
                         nrow = n_delta - 1L, ncol = n_eps - 1L)
    bins <- cut(facet_z[facet_ok], breaks = breaks, include.lowest = TRUE, labels = FALSE)
    facet_cols[facet_ok] <- color_palette[bins]
  } else {
    facet_cols <- col
  }
  
  plotfun <- function() {
    if (isTRUE(restore_par)) {
      op <- par(c("mar", "mgp", "tcl", "cex.axis", "cex.lab", "cex.main"))
      on.exit(par(op), add = TRUE)
    }
    
    if (draw_color_scale) {
      layout(matrix(c(1, 2), nrow = 1), widths = c(1 - color_scale_width, color_scale_width))
      on.exit(layout(1), add = TRUE)
    }
    
    par(
      mar = margins,
      mgp = axis_mgp,
      tcl = tcl,
      cex.axis = tick_cex,
      cex.lab = axis_title_cex,
      cex.main = title_cex
    )
    
    panel_pin <- par("pin")
    panel_aspect <- if (all(is.finite(panel_pin)) && panel_pin[2] > 0) {
      panel_pin[1] / panel_pin[2]
    } else {
      1
    }
    auto_x_stretch <- max(1, panel_aspect)
    auto_y_stretch <- max(1, 1 / panel_aspect)
    x_stretch_use <- if (isTRUE(auto_stretch) && is.null(x_stretch)) auto_x_stretch else (x_stretch %||% 1)
    y_stretch_use <- if (isTRUE(auto_stretch) && is.null(y_stretch)) auto_y_stretch else (y_stretch %||% 1)
    
    delta_plot <- .map_axis_values_to_persp(
      values = delta_values,
      from = range(delta_values),
      stretch = x_stretch_use,
      scale = scale
    )
    eps_plot <- .map_axis_values_to_persp(
      values = eps_values,
      from = range(eps_values),
      stretch = y_stretch_use,
      scale = scale
    )
    z_plot <- .map_axis_values_to_persp(
      values = zmat,
      from = zlim,
      stretch = 1,
      scale = scale
    )
    zlim_plot <- .map_axis_values_to_persp(
      values = zlim,
      from = zlim,
      stretch = 1,
      scale = scale
    )
    
    pmat <- persp(
      x = delta_plot,
      y = eps_plot,
      z = z_plot,
      xlab = "",
      ylab = "",
      zlab = "",
      main = main,
      theta = theta,
      phi = phi,
      r = r,
      d = d,
      scale = FALSE,
      expand = expand,
      shade = shade,
      ltheta = ltheta,
      lphi = lphi,
      col = facet_cols,
      border = border,
      box = box,
      axes = FALSE,
      nticks = nticks,
      ticktype = ticktype,
      zlim = zlim_plot
    )
    
    if (isTRUE(axes)) {
      .draw_manual_persp_axis_ticks(
        axis = "x",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        axis_limits = range(delta_values),
        n = nticks,
        cex = tick_cex,
        lwd = 1,
        tcl = tcl,
        ticktype = ticktype,
        axis_mgp = axis_mgp,
        draw_axis_line = !isTRUE(box)
      )
      .draw_manual_persp_axis_ticks(
        axis = "y",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        axis_limits = range(eps_values),
        n = nticks,
        cex = tick_cex,
        lwd = 1,
        tcl = tcl,
        ticktype = ticktype,
        axis_mgp = axis_mgp,
        draw_axis_line = !isTRUE(box)
      )
      .draw_manual_persp_axis_ticks(
        axis = "z",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        axis_limits = zlim,
        n = nticks,
        cex = tick_cex,
        lwd = 1,
        tcl = tcl,
        ticktype = ticktype,
        axis_mgp = axis_mgp,
        draw_axis_line = !isTRUE(box)
      )
    }
    
    if (isTRUE(show_x_axis_title)) {
      .draw_manual_persp_axis_title(
        label = x_axis_title,
        axis = "x",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        at = x_axis_title_at,
        pad = x_axis_title_pad,
        cex = axis_title_cex
      )
    }
    if (isTRUE(show_y_axis_title)) {
      .draw_manual_persp_axis_title(
        label = y_axis_title,
        axis = "y",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        at = y_axis_title_at,
        pad = y_axis_title_pad,
        cex = axis_title_cex
      )
    }
    if (isTRUE(show_z_axis_title)) {
      .draw_manual_persp_axis_title(
        label = z_axis_title,
        axis = "z",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        at = z_axis_title_at,
        pad = z_axis_title_pad,
        cex = axis_title_cex
      )
    }
    
    if (draw_color_scale) {
      par(
        mar = c(margins[1], 0.35, margins[3], 1.45),
        tcl = tcl,
        cex.lab = axis_title_cex,
        cex.main = title_cex
      )
      .draw_vertical_color_scale_panel(
        zlim = zlim,
        color_palette = color_palette,
        color_scale_title = color_scale_title,
        color_scale_axis_mgp = c(1.2, 0.25, 0),
        color_scale_nticks = color_scale_nticks,
        color_scale_title_cex = color_scale_title_cex,
        color_scale_tick_cex = color_scale_tick_cex,
        color_scale_border = color_scale_border
      )
    }
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_misspecification_cost_contour <- function(
    grid_obj,                                                   # output of misspecification_cost_grid()
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"), # which matrix in grid_obj to plot
    out_dir = file.path("figures", "misspecification"),         # directory where the figure is saved
    out_name = NULL,                                            # base filename for saved outputs
    show = interactive(),                                       # if TRUE, draw on the current device
    save = FALSE,                                               # if TRUE, save PNG and PDF copies
    width_in = 6.67,                                            # saved figure width in inches
    height_in = 4.67,                                           # saved figure height in inches
    dpi = 300,                                                  # PNG resolution used on save
    match_current = FALSE,                                      # if TRUE, save using the current device size instead
    title = TRUE,                                               # if TRUE, use the default plot title
    main = NULL,                                                # custom main title; overrides the default when supplied
    show_x_axis_title = TRUE,                                   # if TRUE, show the x-axis title
    show_y_axis_title = TRUE,                                   # if TRUE, show the y-axis title
    x_axis_title = expression(delta),                           # x-axis title
    y_axis_title = expression(epsilon),                         # y-axis title
    axis_title_cex = 1,                                         # size of the axis titles
    tick_cex = 1.2,                                             # size of the axis tick labels
    title_cex = 1,                                              # size of the main title
    base_cex = 1,                                               # global cex scaling inside the panel
    mex = 1,                                                    # margin expansion factor
    margins = c(3.2, 3.2, 1.8, 0.8),                           # par("mar")
    outer_margins = c(0, 0, 0, 0),                              # par("oma")
    axis_mgp = c(2.2, 0.7, 0),                                 # par("mgp")
    tcl = -0.25,                                                # tick length
    axes = TRUE,                                                # if TRUE, draw axes
    box = TRUE,                                                 # if TRUE, draw a box around the panel
    x_axis_at = NULL,                                           # custom x-axis tick positions
    y_axis_at = NULL,                                           # custom y-axis tick positions
    x_axis_labels = TRUE,                                       # x-axis labels (logical or custom vector)
    y_axis_labels = TRUE,                                       # y-axis labels (logical or custom vector)
    axis_las = 1,                                               # orientation of axis tick labels
    xlim = NULL,                                                # custom x-range
    ylim = NULL,                                                # custom y-range
    zlim = NULL,                                                # custom z-range for colors/contours
    asp = NA,                                                   # aspect ratio; NA leaves it unconstrained
    xaxs = "i",                                                 # x-axis style used by plot.window()
    yaxs = "i",                                                 # y-axis style used by plot.window()
    add_grid = FALSE,                                           # if TRUE, add a 2D reference grid
    grid_col = "#DDDDDD",                                       # grid color
    grid_lty = 3,                                               # grid line type
    grid_lwd = 0.8,                                             # grid line width
    use_z_colormap = FALSE,                                     # if TRUE, draw a filled color background; default leaves the panel white
    color_min = "#1E90FF",                                      # minimum-color anchor
    color_max = "#B22222",                                      # maximum-color anchor
    color_palette = NULL,                                       # custom vector of colors; overrides color_min/color_max
    color_steps = 64,                                           # number of colors in the generated palette
    useRaster = FALSE,                                          # passed to image()
    interpolate = FALSE,                                        # passed to image()
    add_contours = TRUE,                                        # if TRUE, overlay contour lines
    contour_levels = NULL,                                      # explicit contour levels; if NULL, use pretty()
    contour_nlevels = 10,                                       # target number of contour levels when contour_levels is NULL
    contour_drawlabels = TRUE,                                  # if TRUE, label contour lines
    contour_labcex = 0.8,                                       # contour-label size
    contour_method = c("flattest", "simple", "edge"),          # label-placement method for contour()
    contour_lwd = 1.3,                                          # contour-line width
    contour_lty = 1,                                            # contour-line type
    contour_col = "#333333",                                    # contour-line color
    contour_vfont = c("sans serif", "bold"),                    # contour-label vector font
    add_points = FALSE,                                         # if TRUE, show the evaluated (delta, eps) grid points
    point_pch = 16,                                             # point character for the grid points
    point_cex = 0.55,                                           # point size for the grid points
    point_col = "#111111",                                      # point color for the grid points
    add_parametrizer_legend = TRUE,                             # if TRUE and add_contours = TRUE, draw a small broken-line key for the contour level variable
    parametrizer_legend_label = NULL,                           # label shown in the contour-level key; NULL -> surface-specific default
    parametrizer_legend_inset = c(0.04, 0.08),                 # inset of the contour-level key from the lower-left corner
    parametrizer_legend_line_length = 0.18,                     # line length of the contour-level key as a fraction of panel width
    parametrizer_legend_cex = NULL,                             # text size of the contour-level key; NULL -> tick_cex
    parametrizer_legend_col = NULL,                             # color of the contour-level key; NULL -> contour_col
    parametrizer_legend_lty = NULL,                             # line type of the contour-level key; NULL -> contour_lty
    parametrizer_legend_lwd = NULL,                             # line width of the contour-level key; NULL -> contour_lwd
    add_color_scale = FALSE,                                    # if TRUE and use_z_colormap = TRUE, add the color-scale legend
    color_scale_position = c("figure", "right_outer"),          # legend placement mode
    color_scale_title = NULL,                                   # title shown above the color-scale legend
    color_scale_fig = c(0.86, 0.94, 0.26, 0.82),                # legend location in device coordinates when color_scale_position = "figure"
    color_scale_outer_fig = c(0.18, 0.82, 0.06, 0.94),          # legend location within the right outer-margin strip when color_scale_position = "right_outer"
    color_scale_outer_margin = 4.5,                             # minimum right outer margin reserved automatically for "right_outer"
    color_scale_margins = c(0.4, 0.1, 0.8, 1.8),                # margins used inside the legend panel
    color_scale_axis_mgp = c(1.2, 0.25, 0),                    # placement of ticks/title inside the legend panel
    color_scale_nticks = 5,                                     # number of tick marks on the color-scale legend
    color_scale_title_cex = 0.9,                                # size of the color-scale title
    color_scale_tick_cex = 0.85,                                # size of the color-scale tick labels
    color_scale_border = "#444444",                             # border color of the color-scale legend
    restore_par = TRUE                                          # if FALSE, do not restore par() on exit; useful inside a shared layout
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  contour_method <- match.arg(contour_method)
  color_scale_position <- match.arg(color_scale_position)
  surface_meta <- .misspecification_surface_plot_meta(surface)
  
  margins <- .validate_margin_spec(margins, "margins")
  outer_margins <- .validate_margin_spec(outer_margins, "outer_margins")
  color_scale_outer_margin <- as.numeric(color_scale_outer_margin)[1]
  if (!is.finite(color_scale_outer_margin) || color_scale_outer_margin < 0) {
    stop("'color_scale_outer_margin' must be a non-negative number.")
  }
  if (isTRUE(use_z_colormap) && isTRUE(add_color_scale) &&
      identical(color_scale_position, "right_outer")) {
    outer_margins[4] <- max(outer_margins[4], color_scale_outer_margin)
  }
  if (identical(color_scale_position, "figure")) {
    color_scale_fig <- .validate_fig_spec(color_scale_fig, "color_scale_fig")
  } else {
    color_scale_outer_fig <- .validate_fig_spec(color_scale_outer_fig, "color_scale_outer_fig")
  }
  
  color_steps <- as.integer(color_steps)[1]
  if (!is.finite(color_steps) || color_steps < 2L) {
    stop("'color_steps' must be an integer >= 2.")
  }
  contour_nlevels <- as.integer(contour_nlevels)[1]
  if (!is.finite(contour_nlevels) || contour_nlevels < 1L) {
    stop("'contour_nlevels' must be an integer >= 1.")
  }
  
  if (is.null(grid_obj$delta_values) || is.null(grid_obj$eps_values)) {
    stop("grid_obj must be the output of misspecification_cost_grid().")
  }
  
  delta_values <- as.numeric(grid_obj$delta_values)
  eps_values <- as.numeric(grid_obj$eps_values)
  zmat <- grid_obj[[surface]]
  
  zmat <- .ensure_misspecification_surface_has_finite_values(grid_obj, surface)
  if (length(delta_values) < 2L || length(eps_values) < 2L) {
    stop("A contour plot needs at least two delta values and two eps values.")
  }
  if (!isTRUE(use_z_colormap) && !isTRUE(add_contours) && !isTRUE(add_points)) {
    stop("At least one of use_z_colormap, add_contours, or add_points must be TRUE.")
  }
  
  delta_ord <- order(delta_values)
  eps_ord <- order(eps_values)
  delta_values <- delta_values[delta_ord]
  eps_values <- eps_values[eps_ord]
  zmat <- zmat[delta_ord, eps_ord, drop = FALSE]
  zmat[!is.finite(zmat)] <- NA_real_
  
  if (is.null(xlim)) {
    xlim <- range(delta_values, finite = TRUE)
  } else {
    xlim <- .validate_range_spec(xlim, "xlim")
  }
  if (is.null(ylim)) {
    ylim <- range(eps_values, finite = TRUE)
  } else {
    ylim <- .validate_range_spec(ylim, "ylim")
  }
  
  if (is.null(zlim)) {
    zlim <- range(zmat, finite = TRUE)
  } else {
    zlim <- .validate_range_spec(zlim, "zlim")
  }
  if (diff(zlim) == 0) {
    bump <- max(1, abs(zlim[1]))
    zlim <- zlim + c(-1, 1) * 0.05 * bump
  }
  
  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(c(color_min, color_max))(color_steps)
  } else {
    color_palette <- as.character(color_palette)
  }
  if (!length(color_palette)) {
    stop("color_palette must contain at least one color.")
  }
  
  if (is.null(contour_levels)) {
    contour_levels <- pretty(zlim, n = contour_nlevels)
  }
  contour_levels <- sort(unique(as.numeric(contour_levels)))
  contour_levels <- contour_levels[is.finite(contour_levels)]
  contour_levels <- contour_levels[
    contour_levels >= zlim[1] & contour_levels <= zlim[2]
  ]
  if (isTRUE(add_contours) && !length(contour_levels)) {
    stop("No contour levels fall inside 'zlim'.")
  }
  
  if (is.null(color_scale_title)) {
    color_scale_title <- surface_meta$color_scale_title
  }
  if (is.null(parametrizer_legend_label)) {
    parametrizer_legend_label <- surface_meta$color_scale_title
  }
  if (is.null(parametrizer_legend_cex)) {
    parametrizer_legend_cex <- tick_cex
  }
  if (is.null(parametrizer_legend_col)) {
    parametrizer_legend_col <- contour_col
  }
  if (is.null(parametrizer_legend_lty)) {
    parametrizer_legend_lty <- contour_lty
  }
  if (is.null(parametrizer_legend_lwd)) {
    parametrizer_legend_lwd <- contour_lwd
  }
  if (is.null(main)) {
    main <- if (isTRUE(title)) surface_meta$main else ""
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_contour")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  plotfun <- function() {
    if (isTRUE(restore_par)) {
      op <- par(c("oma", "mar", "mgp", "tcl", "cex", "mex",
                  "cex.axis", "cex.lab", "cex.main"))
      on.exit(par(op), add = TRUE)
    }
    
    par_args <- list(
      mar = margins,
      mgp = axis_mgp,
      tcl = tcl,
      cex = base_cex,
      mex = mex,
      cex.axis = tick_cex,
      cex.lab = axis_title_cex,
      cex.main = title_cex
    )
    if (isTRUE(restore_par)) {
      par_args$oma <- outer_margins
    }
    do.call(par, par_args)
    
    plot.new()
    window_args <- list(xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
    if (!is.na(asp)) window_args$asp <- asp
    do.call(plot.window, window_args)
    
    if (isTRUE(use_z_colormap)) {
      image(
        x = delta_values,
        y = eps_values,
        z = zmat,
        zlim = zlim,
        col = color_palette,
        add = TRUE,
        xaxs = "i",
        yaxs = "i",
        useRaster = useRaster,
        interpolate = interpolate
      )
    }
    
    if (isTRUE(add_grid)) {
      grid(col = grid_col, lty = grid_lty, lwd = grid_lwd)
    }
    
    if (isTRUE(add_contours)) {
      contour(
        x = delta_values,
        y = eps_values,
        z = zmat,
        levels = contour_levels,
        add = TRUE,
        drawlabels = contour_drawlabels,
        method = contour_method,
        labcex = contour_labcex,
        lwd = contour_lwd,
        lty = contour_lty,
        col = contour_col,
        vfont = contour_vfont,
        axes = FALSE,
        frame.plot = FALSE
      )
    }
    
    if (isTRUE(add_points)) {
      pts <- expand.grid(delta = delta_values, eps = eps_values)
      points(pts$delta, pts$eps, pch = point_pch, cex = point_cex, col = point_col)
    }
    
    if (isTRUE(axes)) {
      axis(1, at = x_axis_at, labels = x_axis_labels, las = axis_las)
      axis(2, at = y_axis_at, labels = y_axis_labels, las = axis_las)
    }
    if (isTRUE(box)) {
      box()
    }
    
    title(
      main = main,
      xlab = if (show_x_axis_title) x_axis_title else "",
      ylab = if (show_y_axis_title) y_axis_title else ""
    )
    
    if (isTRUE(add_parametrizer_legend) && isTRUE(add_contours)) {
      .draw_parametrizer_legend_key(
        label = parametrizer_legend_label,
        corner = "bottomleft",
        inset = parametrizer_legend_inset,
        line_length = parametrizer_legend_line_length,
        cex = parametrizer_legend_cex,
        col = parametrizer_legend_col,
        lty = parametrizer_legend_lty,
        lwd = parametrizer_legend_lwd
      )
    }
    
    if (isTRUE(use_z_colormap) && isTRUE(add_color_scale)) {
      legend_fig <- .resolve_color_scale_fig(
        color_scale_position = color_scale_position,
        color_scale_fig = color_scale_fig,
        color_scale_outer_fig = color_scale_outer_fig
      )
      .draw_vertical_color_scale(
        zlim = zlim,
        color_palette = color_palette,
        color_scale_title = color_scale_title,
        legend_fig = legend_fig,
        color_scale_margins = color_scale_margins,
        color_scale_axis_mgp = color_scale_axis_mgp,
        color_scale_nticks = color_scale_nticks,
        color_scale_title_cex = color_scale_title_cex,
        color_scale_tick_cex = color_scale_tick_cex,
        color_scale_border = color_scale_border
      )
    }
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_misspecification_cost_vs_delta <- function(
    grid_obj,                                                   # output of misspecification_cost_grid()
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"), # which matrix in grid_obj to plot on the y-axis
    eps_values = NULL,                                          # optional subset of epsilon values to display
    out_dir = file.path("figures", "misspecification"),         # directory where the figure is saved
    out_name = NULL,                                            # base filename for saved outputs
    show = interactive(),                                       # if TRUE, draw on the current device
    save = FALSE,                                               # if TRUE, save PNG and PDF copies
    width_in = 6.67,                                            # saved figure width in inches
    height_in = 4.67,                                           # saved figure height in inches
    dpi = 300,                                                  # PNG resolution used on save
    match_current = FALSE,                                      # if TRUE, save using the current device size instead
    title = TRUE,                                               # if TRUE, use the default plot title
    main = NULL,                                                # custom main title; overrides the default when supplied
    show_x_axis_title = TRUE,                                   # if TRUE, show the x-axis title
    show_y_axis_title = TRUE,                                   # if TRUE, show the y-axis title
    x_axis_title = expression(delta),                           # x-axis title
    y_axis_title = NULL,                                        # y-axis title; if NULL, use the default for the chosen surface
    axis_title_cex = 1,                                         # size of the axis titles
    tick_cex = 1.2,                                             # size of the axis tick labels
    title_cex = 1,                                              # size of the main title
    base_cex = 1,                                               # global cex scaling inside the panel
    mex = 1,                                                    # margin expansion factor
    margins = c(3.2, 3.6, 1.8, 0.8),                           # par("mar")
    outer_margins = c(0, 0, 0, 0),                              # par("oma")
    axis_mgp = c(2.2, 0.7, 0),                                 # par("mgp")
    tcl = -0.25,                                                # tick length
    axes = TRUE,                                                # if TRUE, draw axes
    box = TRUE,                                                 # if TRUE, draw a box around the panel
    x_axis_at = NULL,                                           # custom x-axis tick positions
    y_axis_at = NULL,                                           # custom y-axis tick positions
    x_axis_labels = TRUE,                                       # x-axis labels (logical or custom vector)
    y_axis_labels = TRUE,                                       # y-axis labels (logical or custom vector)
    axis_las = 1,                                               # orientation of axis tick labels
    xlim = NULL,                                                # custom x-range
    ylim = NULL,                                                # custom y-range
    xaxs = "i",                                                 # x-axis style used by plot.window()
    yaxs = "i",                                                 # y-axis style used by plot.window()
    add_grid = FALSE,                                           # if TRUE, add a 2D reference grid
    grid_col = "#DDDDDD",                                       # grid color
    grid_lty = 3,                                               # grid line type
    grid_lwd = 0.8,                                             # grid line width
    type = c("l", "p", "b", "o"),                              # line/point style used for the epsilon curves
    cols = NULL,                                                # line colors; generated automatically when NULL
    ltys = 1,                                                   # line types for the epsilon curves
    lwds = 2,                                                   # line widths for the epsilon curves
    add_points = FALSE,                                         # if TRUE, add points even when type = "l"
    point_pch = 16,                                             # point character for the epsilon curves
    point_cex = 0.7,                                            # point size for the epsilon curves
    point_cols = NULL,                                          # point colors; defaults to cols
    label_eps_on_curve = TRUE,                                  # if TRUE, annotate each curve with epsilon directly on the line
    label_x = 0.9,                                              # target delta value where epsilon labels are placed
    label_gap = NULL,                                           # horizontal gap left in the line around the label; NULL -> automatic
    label_cex = 0.85,                                           # size of the inline epsilon labels
    label_col = NULL,                                           # label colors; defaults to the curve colors
    add_parametrizer_legend = TRUE,                             # if TRUE, draw a small broken-line key for the parametrizing variable
    add_legend = FALSE,                                         # if TRUE, show the epsilon legend
    legend_title = expression(epsilon),                         # legend title
    legend_labels = NULL,                                       # custom legend labels; if NULL, use the epsilon values
    legend_position = "topleft",                                # legend position keyword or c(x, y)
    legend_inset = 0.01,                                        # legend inset when using a keyword position
    legend_ncol = 1,                                            # number of legend columns
    legend_cex = 0.9,                                           # legend text size
    legend_bty = "n",                                           # legend box type
    legend_seg_len = 2.5,                                       # legend line-segment length
    restore_par = TRUE                                          # if FALSE, do not restore par() on exit; useful inside a shared layout
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  type <- match.arg(type)
  surface_meta <- .misspecification_surface_plot_meta(surface)
  
  margins <- .validate_margin_spec(margins, "margins")
  outer_margins <- .validate_margin_spec(outer_margins, "outer_margins")
  
  if (is.null(grid_obj$delta_values) || is.null(grid_obj$eps_values)) {
    stop("grid_obj must be the output of misspecification_cost_grid().")
  }
  
  delta_values <- as.numeric(grid_obj$delta_values)
  eps_grid <- as.numeric(grid_obj$eps_values)
  zmat <- grid_obj[[surface]]
  
  zmat <- .ensure_misspecification_surface_has_finite_values(grid_obj, surface)
  if (!length(delta_values) || !length(eps_grid)) {
    stop("The selected grid must contain at least one delta value and one eps value.")
  }
  
  delta_ord <- order(delta_values)
  eps_ord <- order(eps_grid)
  delta_values <- delta_values[delta_ord]
  eps_grid <- eps_grid[eps_ord]
  zmat <- zmat[delta_ord, eps_ord, drop = FALSE]
  zmat[!is.finite(zmat)] <- NA_real_
  
  if (is.null(eps_values)) {
    eps_plot <- eps_grid
  } else {
    eps_values <- as.numeric(eps_values)
    if (!length(eps_values) || any(!is.finite(eps_values)) ||
        any(eps_values < 0 | eps_values > 1)) {
      stop("'eps_values' must be a non-empty numeric vector with entries in [0, 1].")
    }
    
    tol <- sqrt(.Machine$double.eps)
    idx <- vapply(
      eps_values,
      function(v) {
        hits <- which(abs(eps_grid - v) <= tol * max(1, abs(v)))
        if (length(hits)) hits[1] else NA_integer_
      },
      integer(1)
    )
    if (anyNA(idx)) {
      missing_eps <- format(signif(eps_values[is.na(idx)], 8), trim = TRUE, scientific = FALSE)
      stop(sprintf("The following eps_values were not found in grid_obj: %s",
                   paste(missing_eps, collapse = ", ")))
    }
    
    eps_plot <- eps_grid[idx]
    zmat <- zmat[, idx, drop = FALSE]
  }
  
  keep_curve <- apply(zmat, 2, function(y) any(is.finite(y)))
  if (!all(keep_curve)) {
    warning(sprintf("Dropping %d epsilon curve(s) with no finite values.", sum(!keep_curve)))
    eps_plot <- eps_plot[keep_curve]
    zmat <- zmat[, keep_curve, drop = FALSE]
  }
  if (!length(eps_plot) || !any(is.finite(zmat))) {
    stop("The selected surface does not contain any finite values to plot.")
  }
  
  if (is.null(xlim)) {
    xlim <- range(delta_values, finite = TRUE)
  } else {
    xlim <- .validate_range_spec(xlim, "xlim")
  }
  if (is.null(ylim)) {
    ylim <- range(zmat, finite = TRUE)
  } else {
    ylim <- .validate_range_spec(ylim, "ylim")
  }
  if (diff(ylim) == 0) {
    bump <- max(1, abs(ylim[1]))
    ylim <- ylim + c(-1, 1) * 0.05 * bump
  }
  
  label_x <- as.numeric(label_x)[1]
  if (!is.finite(label_x)) {
    stop("'label_x' must be a single finite number.")
  }
  if (!is.null(label_gap)) {
    label_gap <- as.numeric(label_gap)[1]
    if (!is.finite(label_gap) || label_gap < 0) {
      stop("'label_gap' must be NULL or a single non-negative number.")
    }
  }
  
  n_curves <- ncol(zmat)
  if (is.null(cols)) {
    cols <- tryCatch(
      grDevices::hcl.colors(n_curves, palette = "Dark 3"),
      error = function(e) grDevices::rainbow(n_curves, end = 0.85)
    )
  }
  cols <- rep_len(cols, n_curves)
  ltys <- rep_len(ltys, n_curves)
  lwds <- rep_len(lwds, n_curves)
  point_pch <- rep_len(point_pch, n_curves)
  point_cex <- rep_len(point_cex, n_curves)
  point_cols <- rep_len(point_cols %||% cols, n_curves)
  label_cols <- rep_len(label_col %||% cols, n_curves)
  
  if (is.null(legend_labels)) {
    legend_labels <- format(signif(eps_plot, 8), trim = TRUE, scientific = FALSE)
  } else if (length(legend_labels) != n_curves) {
    stop("'legend_labels' must have the same length as the number of plotted epsilon curves.")
  }
  
  if (is.null(y_axis_title)) {
    y_axis_title <- surface_meta$color_scale_title
  }
  if (is.null(main)) {
    main <- if (isTRUE(title)) paste(surface_meta$main, "vs delta") else ""
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_vs_delta")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  show_lines <- type %in% c("l", "b", "o")
  show_points <- isTRUE(add_points) || type %in% c("p", "b", "o")
  
  plotfun <- function() {
    if (isTRUE(restore_par)) {
      op <- par(c("oma", "mar", "mgp", "tcl", "cex", "mex",
                  "cex.axis", "cex.lab", "cex.main"))
      on.exit(par(op), add = TRUE)
    }
    
    par_args <- list(
      mar = margins,
      mgp = axis_mgp,
      tcl = tcl,
      cex = base_cex,
      mex = mex,
      cex.axis = tick_cex,
      cex.lab = axis_title_cex,
      cex.main = title_cex
    )
    if (isTRUE(restore_par)) {
      par_args$oma <- outer_margins
    }
    do.call(par, par_args)
    
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
    
    if (isTRUE(add_grid)) {
      grid(col = grid_col, lty = grid_lty, lwd = grid_lwd)
    }
    
    for (j in seq_len(n_curves)) {
      yj <- zmat[, j]
      okj <- is.finite(yj)
      if (!any(okj)) next
      
      xj <- delta_values[okj]
      yj <- yj[okj]
      label_info <- NULL
      if (isTRUE(label_eps_on_curve)) {
        gap_width_j <- label_gap
        if (is.null(gap_width_j)) {
          gap_width_j <- max(
            1.15 * strwidth(legend_labels[j], cex = label_cex, units = "user"),
            0.015 * diff(xlim)
          )
        }
        label_info <- .curve_label_info(
          x = xj,
          y = yj,
          target_x = label_x,
          gap_width = gap_width_j
        )
      }
      
      if (show_lines) {
        if (is.null(label_info)) {
          lines(xj, yj, col = cols[j], lty = ltys[j], lwd = lwds[j])
        } else {
          .draw_curve_with_gap(
            x = xj,
            y = yj,
            col = cols[j],
            lty = ltys[j],
            lwd = lwds[j],
            x_gap_left = label_info$left_x,
            x_gap_right = label_info$right_x
          )
        }
      }
      if (show_points) {
        if (!is.null(label_info)) {
          keep_pts <- (xj < label_info$left_x) | (xj > label_info$right_x)
          x_pts <- xj[keep_pts]
          y_pts <- yj[keep_pts]
        } else {
          x_pts <- xj
          y_pts <- yj
        }
        if (length(x_pts)) {
          points(x_pts, y_pts, pch = point_pch[j], cex = point_cex[j], col = point_cols[j])
        }
      }
      if (!is.null(label_info)) {
        text(
          x = label_info$x,
          y = label_info$y,
          labels = legend_labels[j],
          cex = label_cex,
          col = label_cols[j],
          xpd = NA
        )
      }
    }
    
    if (isTRUE(axes)) {
      axis(1, at = x_axis_at, labels = x_axis_labels, las = axis_las)
      axis(2, at = y_axis_at, labels = y_axis_labels, las = axis_las)
    }
    if (isTRUE(box)) {
      box()
    }
    
    title(
      main = main,
      xlab = if (show_x_axis_title) x_axis_title else "",
      ylab = if (show_y_axis_title) y_axis_title else ""
    )
    
    if (isTRUE(add_parametrizer_legend)) {
      .draw_parametrizer_legend_key(
        label = expression(epsilon),
        cex = tick_cex
      )
    }
    
    if (isTRUE(add_legend)) {
      legend_args <- list(
        legend = legend_labels,
        title = legend_title,
        col = cols,
        ncol = legend_ncol,
        cex = legend_cex,
        bty = legend_bty,
        seg.len = legend_seg_len,
        merge = TRUE
      )
      
      if (show_lines) {
        legend_args$lty <- ltys
        legend_args$lwd <- lwds
      }
      if (show_points) {
        legend_args$pch <- point_pch
        legend_args$pt.cex <- point_cex
      }
      
      if (is.character(legend_position) && length(legend_position) == 1L) {
        legend_args$x <- legend_position
        legend_args$inset <- legend_inset
      } else if (is.numeric(legend_position) && length(legend_position) == 2L &&
                 all(is.finite(legend_position))) {
        legend_args$x <- legend_position[1]
        legend_args$y <- legend_position[2]
      } else {
        stop("'legend_position' must be a single legend keyword or a numeric vector c(x, y).")
      }
      
      do.call(legend, legend_args)
    }
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_misspecification_cost_vs_epsilon <- function(
    grid_obj,                                                   # output of misspecification_cost_grid()
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"), # which matrix in grid_obj to plot on the y-axis
    delta_values = NULL,                                        # optional subset of delta values to display
    out_dir = file.path("figures", "misspecification"),         # directory where the figure is saved
    out_name = NULL,                                            # base filename for saved outputs
    show = interactive(),                                       # if TRUE, draw on the current device
    save = FALSE,                                               # if TRUE, save PNG and PDF copies
    width_in = 6.67,                                            # saved figure width in inches
    height_in = 4.67,                                           # saved figure height in inches
    dpi = 300,                                                  # PNG resolution used on save
    match_current = FALSE,                                      # if TRUE, save using the current device size instead
    title = TRUE,                                               # if TRUE, use the default plot title
    main = NULL,                                                # custom main title; overrides the default when supplied
    show_x_axis_title = TRUE,                                   # if TRUE, show the x-axis title
    show_y_axis_title = TRUE,                                   # if TRUE, show the y-axis title
    x_axis_title = expression(epsilon),                         # x-axis title
    y_axis_title = NULL,                                        # y-axis title; if NULL, use the default for the chosen surface
    axis_title_cex = 1,                                         # size of the axis titles
    tick_cex = 1.2,                                             # size of the axis tick labels
    title_cex = 1,                                              # size of the main title
    base_cex = 1,                                               # global cex scaling inside the panel
    mex = 1,                                                    # margin expansion factor
    margins = c(3.2, 3.6, 1.8, 0.8),                           # par("mar")
    outer_margins = c(0, 0, 0, 0),                              # par("oma")
    axis_mgp = c(2.2, 0.7, 0),                                 # par("mgp")
    tcl = -0.25,                                                # tick length
    axes = TRUE,                                                # if TRUE, draw axes
    box = TRUE,                                                 # if TRUE, draw a box around the panel
    x_axis_at = NULL,                                           # custom x-axis tick positions
    y_axis_at = NULL,                                           # custom y-axis tick positions
    x_axis_labels = TRUE,                                       # x-axis labels (logical or custom vector)
    y_axis_labels = TRUE,                                       # y-axis labels (logical or custom vector)
    axis_las = 1,                                               # orientation of axis tick labels
    xlim = NULL,                                                # custom x-range
    ylim = NULL,                                                # custom y-range
    xaxs = "i",                                                 # x-axis style used by plot.window()
    yaxs = "i",                                                 # y-axis style used by plot.window()
    add_grid = FALSE,                                           # if TRUE, add a 2D reference grid
    grid_col = "#DDDDDD",                                       # grid color
    grid_lty = 3,                                               # grid line type
    grid_lwd = 0.8,                                             # grid line width
    type = c("l", "p", "b", "o"),                              # line/point style used for the delta curves
    cols = NULL,                                                # line colors; generated automatically when NULL
    ltys = 1,                                                   # line types for the delta curves
    lwds = 2,                                                   # line widths for the delta curves
    add_points = FALSE,                                         # if TRUE, add points even when type = "l"
    point_pch = 16,                                             # point character for the delta curves
    point_cex = 0.7,                                            # point size for the delta curves
    point_cols = NULL,                                          # point colors; defaults to cols
    label_delta_on_curve = TRUE,                                # if TRUE, annotate each curve with delta directly on the line
    label_x = 0.9,                                              # target epsilon value where delta labels are placed
    label_gap = NULL,                                           # horizontal gap left in the line around the label; NULL -> automatic
    label_cex = 0.85,                                           # size of the inline delta labels
    label_col = NULL,                                           # label colors; defaults to the curve colors
    add_parametrizer_legend = TRUE,                             # if TRUE, draw a small broken-line key for the parametrizing variable
    add_legend = FALSE,                                         # if TRUE, show the delta legend
    legend_title = expression(delta),                           # legend title
    legend_labels = NULL,                                       # custom legend labels; if NULL, use the delta values
    legend_position = "topleft",                                # legend position keyword or c(x, y)
    legend_inset = 0.01,                                        # legend inset when using a keyword position
    legend_ncol = 1,                                            # number of legend columns
    legend_cex = 0.9,                                           # legend text size
    legend_bty = "n",                                           # legend box type
    legend_seg_len = 2.5,                                       # legend line-segment length
    restore_par = TRUE                                          # if FALSE, do not restore par() on exit; useful inside a shared layout
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  type <- match.arg(type)
  surface_meta <- .misspecification_surface_plot_meta(surface)
  
  margins <- .validate_margin_spec(margins, "margins")
  outer_margins <- .validate_margin_spec(outer_margins, "outer_margins")
  
  if (is.null(grid_obj$delta_values) || is.null(grid_obj$eps_values)) {
    stop("grid_obj must be the output of misspecification_cost_grid().")
  }
  
  delta_grid <- as.numeric(grid_obj$delta_values)
  eps_values <- as.numeric(grid_obj$eps_values)
  zmat <- grid_obj[[surface]]
  
  zmat <- .ensure_misspecification_surface_has_finite_values(grid_obj, surface)
  if (!length(delta_grid) || !length(eps_values)) {
    stop("The selected grid must contain at least one delta value and one eps value.")
  }
  
  delta_ord <- order(delta_grid)
  eps_ord <- order(eps_values)
  delta_grid <- delta_grid[delta_ord]
  eps_values <- eps_values[eps_ord]
  zmat <- zmat[delta_ord, eps_ord, drop = FALSE]
  zmat[!is.finite(zmat)] <- NA_real_
  
  if (is.null(delta_values)) {
    delta_plot <- delta_grid
  } else {
    delta_values <- as.numeric(delta_values)
    if (!length(delta_values) || any(!is.finite(delta_values)) || any(delta_values < 0)) {
      stop("'delta_values' must be a non-empty numeric vector with non-negative entries.")
    }
    
    tol <- sqrt(.Machine$double.eps)
    idx <- vapply(
      delta_values,
      function(v) {
        hits <- which(abs(delta_grid - v) <= tol * max(1, abs(v)))
        if (length(hits)) hits[1] else NA_integer_
      },
      integer(1)
    )
    if (anyNA(idx)) {
      missing_delta <- format(signif(delta_values[is.na(idx)], 8), trim = TRUE, scientific = FALSE)
      stop(sprintf("The following delta_values were not found in grid_obj: %s",
                   paste(missing_delta, collapse = ", ")))
    }
    
    delta_plot <- delta_grid[idx]
    zmat <- zmat[idx, , drop = FALSE]
  }
  
  keep_curve <- apply(zmat, 1, function(y) any(is.finite(y)))
  if (!all(keep_curve)) {
    warning(sprintf("Dropping %d delta curve(s) with no finite values.", sum(!keep_curve)))
    delta_plot <- delta_plot[keep_curve]
    zmat <- zmat[keep_curve, , drop = FALSE]
  }
  if (!length(delta_plot) || !any(is.finite(zmat))) {
    stop("The selected surface does not contain any finite values to plot.")
  }
  
  if (is.null(xlim)) {
    xlim <- range(eps_values, finite = TRUE)
  } else {
    xlim <- .validate_range_spec(xlim, "xlim")
  }
  if (is.null(ylim)) {
    ylim <- range(zmat, finite = TRUE)
  } else {
    ylim <- .validate_range_spec(ylim, "ylim")
  }
  if (diff(ylim) == 0) {
    bump <- max(1, abs(ylim[1]))
    ylim <- ylim + c(-1, 1) * 0.05 * bump
  }
  
  label_x <- as.numeric(label_x)[1]
  if (!is.finite(label_x)) {
    stop("'label_x' must be a single finite number.")
  }
  if (!is.null(label_gap)) {
    label_gap <- as.numeric(label_gap)[1]
    if (!is.finite(label_gap) || label_gap < 0) {
      stop("'label_gap' must be NULL or a single non-negative number.")
    }
  }
  
  n_curves <- nrow(zmat)
  if (is.null(cols)) {
    cols <- tryCatch(
      grDevices::hcl.colors(n_curves, palette = "Dark 3"),
      error = function(e) grDevices::rainbow(n_curves, end = 0.85)
    )
  }
  cols <- rep_len(cols, n_curves)
  ltys <- rep_len(ltys, n_curves)
  lwds <- rep_len(lwds, n_curves)
  point_pch <- rep_len(point_pch, n_curves)
  point_cex <- rep_len(point_cex, n_curves)
  point_cols <- rep_len(point_cols %||% cols, n_curves)
  label_cols <- rep_len(label_col %||% cols, n_curves)
  
  if (is.null(legend_labels)) {
    legend_labels <- format(signif(delta_plot, 8), trim = TRUE, scientific = FALSE)
  } else if (length(legend_labels) != n_curves) {
    stop("'legend_labels' must have the same length as the number of plotted delta curves.")
  }
  
  if (is.null(y_axis_title)) {
    y_axis_title <- surface_meta$color_scale_title
  }
  if (is.null(main)) {
    main <- if (isTRUE(title)) paste(surface_meta$main, "vs epsilon") else ""
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_vs_epsilon")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  show_lines <- type %in% c("l", "b", "o")
  show_points <- isTRUE(add_points) || type %in% c("p", "b", "o")
  
  plotfun <- function() {
    if (isTRUE(restore_par)) {
      op <- par(c("oma", "mar", "mgp", "tcl", "cex", "mex",
                  "cex.axis", "cex.lab", "cex.main"))
      on.exit(par(op), add = TRUE)
    }
    
    par_args <- list(
      mar = margins,
      mgp = axis_mgp,
      tcl = tcl,
      cex = base_cex,
      mex = mex,
      cex.axis = tick_cex,
      cex.lab = axis_title_cex,
      cex.main = title_cex
    )
    if (isTRUE(restore_par)) {
      par_args$oma <- outer_margins
    }
    do.call(par, par_args)
    
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
    
    if (isTRUE(add_grid)) {
      grid(col = grid_col, lty = grid_lty, lwd = grid_lwd)
    }
    
    for (j in seq_len(n_curves)) {
      yj <- zmat[j, ]
      okj <- is.finite(yj)
      if (!any(okj)) next
      
      xj <- eps_values[okj]
      yj <- yj[okj]
      label_info <- NULL
      if (isTRUE(label_delta_on_curve)) {
        gap_width_j <- label_gap
        if (is.null(gap_width_j)) {
          gap_width_j <- max(
            1.15 * strwidth(legend_labels[j], cex = label_cex, units = "user"),
            0.015 * diff(xlim)
          )
        }
        label_info <- .curve_label_info(
          x = xj,
          y = yj,
          target_x = label_x,
          gap_width = gap_width_j
        )
      }
      
      if (show_lines) {
        if (is.null(label_info)) {
          lines(xj, yj, col = cols[j], lty = ltys[j], lwd = lwds[j])
        } else {
          .draw_curve_with_gap(
            x = xj,
            y = yj,
            col = cols[j],
            lty = ltys[j],
            lwd = lwds[j],
            x_gap_left = label_info$left_x,
            x_gap_right = label_info$right_x
          )
        }
      }
      if (show_points) {
        if (!is.null(label_info)) {
          keep_pts <- (xj < label_info$left_x) | (xj > label_info$right_x)
          x_pts <- xj[keep_pts]
          y_pts <- yj[keep_pts]
        } else {
          x_pts <- xj
          y_pts <- yj
        }
        if (length(x_pts)) {
          points(x_pts, y_pts, pch = point_pch[j], cex = point_cex[j], col = point_cols[j])
        }
      }
      if (!is.null(label_info)) {
        text(
          x = label_info$x,
          y = label_info$y,
          labels = legend_labels[j],
          cex = label_cex,
          col = label_cols[j],
          xpd = NA
        )
      }
    }
    
    if (isTRUE(axes)) {
      axis(1, at = x_axis_at, labels = x_axis_labels, las = axis_las)
      axis(2, at = y_axis_at, labels = y_axis_labels, las = axis_las)
    }
    if (isTRUE(box)) {
      box()
    }
    
    title(
      main = main,
      xlab = if (show_x_axis_title) x_axis_title else "",
      ylab = if (show_y_axis_title) y_axis_title else ""
    )
    
    if (isTRUE(add_parametrizer_legend)) {
      .draw_parametrizer_legend_key(
        label = expression(delta),
        cex = tick_cex
      )
    }
    
    if (isTRUE(add_legend)) {
      legend_args <- list(
        legend = legend_labels,
        title = legend_title,
        col = cols,
        ncol = legend_ncol,
        cex = legend_cex,
        bty = legend_bty,
        seg.len = legend_seg_len,
        merge = TRUE
      )
      
      if (show_lines) {
        legend_args$lty <- ltys
        legend_args$lwd <- lwds
      }
      if (show_points) {
        legend_args$pch <- point_pch
        legend_args$pt.cex <- point_cex
      }
      
      if (is.character(legend_position) && length(legend_position) == 1L) {
        legend_args$x <- legend_position
        legend_args$inset <- legend_inset
      } else if (is.numeric(legend_position) && length(legend_position) == 2L &&
                 all(is.finite(legend_position))) {
        legend_args$x <- legend_position[1]
        legend_args$y <- legend_position[2]
      } else {
        stop("'legend_position' must be a single legend keyword or a numeric vector c(x, y).")
      }
      
      do.call(legend, legend_args)
    }
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_misspecification_cost_canvas <- function(
    grid_obj,
    grid_obj_vs_delta = grid_obj,
    grid_obj_vs_epsilon = grid_obj,
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"),
    out_dir = file.path("figures", "misspecification"),
    out_name = NULL,
    show = interactive(),
    save = FALSE,
    width_in = 12,
    height_in = 8.5,
    dpi = 300,
    match_current = FALSE,
    panel_layout = c("1x3", "3x1"),
    column_widths = NULL,
    row_heights = NULL,
    contour_args = list(),
    vs_delta_args = list(),
    vs_epsilon_args = list()
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  panel_layout <- match.arg(panel_layout)
  
  validate_positive_vector <- function(x, name, n) {
    x <- as.numeric(x)
    if (length(x) != n || any(!is.finite(x)) || any(x <= 0)) {
      stop(sprintf("'%s' must be a numeric vector of length %d with positive entries.", name, n))
    }
    x
  }
  validate_list <- function(x, name) {
    if (!is.list(x)) {
      stop(sprintf("'%s' must be a list.", name))
    }
    x
  }
  
  if (identical(panel_layout, "1x3")) {
    if (is.null(column_widths)) column_widths <- c(1, 1, 1)
    if (is.null(row_heights)) row_heights <- 1
    column_widths <- validate_positive_vector(column_widths, "column_widths", 3L)
    row_heights <- validate_positive_vector(row_heights, "row_heights", 1L)
  } else {
    if (is.null(column_widths)) column_widths <- 1
    if (is.null(row_heights)) row_heights <- c(1, 1, 1)
    column_widths <- validate_positive_vector(column_widths, "column_widths", 1L)
    row_heights <- validate_positive_vector(row_heights, "row_heights", 3L)
  }
  contour_args <- validate_list(contour_args, "contour_args")
  vs_delta_args <- validate_list(vs_delta_args, "vs_delta_args")
  vs_epsilon_args <- validate_list(vs_epsilon_args, "vs_epsilon_args")
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_canvas")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  build_panel_args <- function(defaults, extras, grid) {
    args <- utils::modifyList(defaults, extras)
    args$grid_obj <- grid
    args$surface <- surface
    args$out_dir <- out_dir
    args$show <- TRUE
    args$save <- FALSE
    args$title <- FALSE
    args$restore_par <- FALSE
    args
  }
  
  contour_call <- build_panel_args(
    defaults = list(add_color_scale = FALSE),
    extras = contour_args,
    grid = grid_obj
  )
  contour_call$add_color_scale <- FALSE
  
  vs_delta_call <- build_panel_args(
    defaults = list(),
    extras = vs_delta_args,
    grid = grid_obj_vs_delta
  )
  vs_epsilon_call <- build_panel_args(
    defaults = list(),
    extras = vs_epsilon_args,
    grid = grid_obj_vs_epsilon
  )
  
  plotfun <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    layout_mat <- if (identical(panel_layout, "1x3")) {
      matrix(1:3, nrow = 1, byrow = TRUE)
    } else {
      matrix(1:3, ncol = 1, byrow = TRUE)
    }
    layout(layout_mat, widths = column_widths, heights = row_heights)
    on.exit(layout(1), add = TRUE)
    
    do.call(plot_misspecification_cost_vs_delta, vs_delta_call)
    do.call(plot_misspecification_cost_contour, contour_call)
    do.call(plot_misspecification_cost_vs_epsilon, vs_epsilon_call)
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_misspecification_cost_example_set <- function(
    tag = NULL,
    b_value = NULL,
    grids = NULL,
    r = 1,
    sigma = 1,
    mu = 1,
    u = 1,
    l = 1,
    surface = c("relative_cost_pct", "relative_cost"),
    table_surface = NULL,
    vs_delta_label_x = 0.9,
    canvas_panel_layout = "1x3",
    canvas_column_widths = NULL,
    canvas_row_heights = NULL,
    canvas_width_in = NULL,
    canvas_height_in = NULL
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  if (is.null(table_surface)) {
    table_surface <- surface
  }
  table_surface <- .misspecification_reporting_surface(
    match.arg(table_surface, c("relative_cost_pct", "relative_cost"))
  )
  
  if (!is.null(grids)) {
    grid_params <- .misspecification_grid_model_params(grids$main_grid)
    if (!is.null(grid_params)) {
      b_value <- grid_params$b
      r <- grid_params$r
      sigma <- grid_params$sigma
      mu <- grid_params$mu
      u <- grid_params$u
      l <- grid_params$l
    }
  }
  
  tag <- .resolve_misspecification_example_tag(tag = tag, b_value = b_value, r = r)
  
  if (is.null(grids)) {
    grids <- load_misspecification_cost_example_data(
      tag = tag,
      b_value = b_value,
      r = r, sigma = sigma, mu = mu, u = u, l = l
    )
  }
  
  spec <- .misspecification_cost_example_spec(
    tag = tag,
    b_value = b_value,
    r = r, sigma = sigma, mu = mu, u = u, l = l
  )
  base_out_dir <- spec$out_dir
  base_name <- spec$base_name
  persp_r <- 0.5
  
  misspec_grid <- grids$main_grid
  misspec_grid_cost_vs_delta <- grids$grid_vs_delta
  misspec_grid_cost_vs_epsilon <- grids$grid_vs_epsilon
  misspec_grid_text <- grids$text_grid
  
  if (is.null(canvas_width_in)) {
    canvas_width_in <- if (identical(canvas_panel_layout, "3x1")) 4.67 else 12
  }
  if (is.null(canvas_height_in)) {
    canvas_height_in <- if (identical(canvas_panel_layout, "3x1")) 12 else 4.67
  }
  if (is.null(canvas_column_widths)) {
    canvas_column_widths <- if (identical(canvas_panel_layout, "3x1")) 1 else c(1, 1, 1)
  }
  if (is.null(canvas_row_heights)) {
    canvas_row_heights <- if (identical(canvas_panel_layout, "3x1")) c(1, 1, 1) else 1
  }
  
  plot_misspecification_cost_3d(
    misspec_grid,
    surface = surface,
    show = TRUE,
    save = TRUE,
    add_color_scale = FALSE,
    height_in = 4.67, width_in = 6.67,
    margins = c(0.9, 2.1, 0.0, 0.0),
    x_axis_title_pad = 0.12,
    y_axis_title_pad = 0.12,
    z_axis_title_pad = 0.1,
    x_axis_title_at = 0.5,
    y_axis_title_at = 0.5,
    z_axis_title_at = 0.5,
    r = persp_r,
    expand = 0.65,
    title = FALSE,
    tick_cex = 1,
    axis_title_cex = 1.4
  )
  plot_misspecification_cost_contour(
    misspec_grid,
    surface = surface,
    show = TRUE,
    save = TRUE,
    title = FALSE,
    show_x_axis_title = TRUE,
    show_y_axis_title = TRUE,
    axis_title_cex = 1.4,
    tick_cex = 1,
    margins = c(3.5, 3.8, 1.5, 1),
    axis_mgp = c(2.6, 0.8, 0),
    outer_margins = c(0, 0, 0, 0),
    add_grid = TRUE,
    contour_nlevels = 12,
    use_z_colormap = FALSE,
    add_color_scale = FALSE,
    contour_lwd = 1.5, contour_labcex = 1
  )
  plot_misspecification_cost_vs_delta(
    misspec_grid_cost_vs_delta,
    surface = surface,
    show = TRUE,
    save = TRUE,
    title = FALSE,
    show_x_axis_title = TRUE,
    show_y_axis_title = TRUE,
    axis_title_cex = 1.4,
    tick_cex = 1,
    margins = c(3.8, 3.8, 1.5, 1),
    axis_mgp = c(2.6, 0.8, 0),
    add_grid = TRUE,
    label_eps_on_curve = TRUE,
    label_x = vs_delta_label_x,
    add_legend = FALSE, cols = "black"
  )
  plot_misspecification_cost_vs_epsilon(
    misspec_grid_cost_vs_epsilon,
    surface = surface,
    show = TRUE,
    save = TRUE,
    title = FALSE,
    show_x_axis_title = TRUE,
    show_y_axis_title = TRUE,
    axis_title_cex = 1.4,
    tick_cex = 1,
    margins = c(3.8, 3.8, 1.5, 1),
    axis_mgp = c(2.6, 0.8, 0),
    add_grid = TRUE,
    label_delta_on_curve = TRUE,
    label_x = 0.9,
    add_legend = FALSE, cols = "black"
  )
  plot_misspecification_cost_canvas(
    misspec_grid,
    grid_obj_vs_delta = misspec_grid_cost_vs_delta,
    grid_obj_vs_epsilon = misspec_grid_cost_vs_epsilon,
    surface = surface,
    out_dir = base_out_dir,
    out_name = paste0(base_name, "_", surface, "_canvas"),
    show = TRUE,
    save = TRUE,
    width_in = canvas_width_in,
    height_in = canvas_height_in,
    panel_layout = canvas_panel_layout,
    column_widths = canvas_column_widths,
    row_heights = canvas_row_heights,
    contour_args = list(
      show_x_axis_title = TRUE,
      show_y_axis_title = TRUE,
      axis_title_cex = 1.1,
      tick_cex = 0.9,
      margins = c(3.2, 3.4, 0.8, 0.5),
      axis_mgp = c(2.2, 0.7, 0),
      outer_margins = c(0, 0, 0, 0),
      add_grid = TRUE,
      contour_nlevels = 12,
      contour_lwd = 1.3,
      contour_labcex = 0.9
    ),
    vs_delta_args = list(
      show_x_axis_title = TRUE,
      show_y_axis_title = TRUE,
      axis_title_cex = 1.1,
      tick_cex = 0.9,
      margins = c(3.2, 3.4, 0.8, 0.5),
      axis_mgp = c(2.2, 0.7, 0),
      outer_margins = c(0, 0, 0, 0),
      add_grid = TRUE,
      label_eps_on_curve = TRUE,
      label_x = vs_delta_label_x,
      add_legend = FALSE,
      cols = "black"
    ),
    vs_epsilon_args = list(
      show_x_axis_title = TRUE,
      show_y_axis_title = TRUE,
      axis_title_cex = 1.1,
      tick_cex = 0.9,
      margins = c(3.2, 3.4, 0.8, 0.5),
      axis_mgp = c(2.2, 0.7, 0),
      outer_margins = c(0, 0, 0, 0),
      add_grid = TRUE,
      label_delta_on_curve = TRUE,
      label_x = 0.9,
      add_legend = FALSE,
      cols = "black"
    )
  )
  
  tex_obj <- export_misspecification_cost_latex(
    misspec_grid_text,
    surface = table_surface,
    out_dir = base_out_dir,
    out_name = paste0("rmc_table_", tag),
    save = TRUE
  )
  
  invisible(list(
    main_grid = misspec_grid,
    grid_vs_delta = misspec_grid_cost_vs_delta,
    grid_vs_epsilon = misspec_grid_cost_vs_epsilon,
    text_grid = misspec_grid_text,
    tex = tex_obj
  ))
}

