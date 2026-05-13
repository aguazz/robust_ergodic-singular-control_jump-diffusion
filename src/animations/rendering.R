# Shared animation rendering helpers.
#
# Notation:
#   X      : reflected state path.
#   U      : cumulative upward reflection at the lower barrier xL.
#   L      : cumulative downward reflection at the upper barrier xU. Animation
#            labels often call this D_t to emphasize downward intervention.
#   J_t    : finite-horizon running-average cost.
#   xL, xU : lower and upper reflecting barriers.
#   xk, xl : drift- and intensity-ambiguity thresholds.
#   gamma  : ergodic value returned by the solver.
#
# Loading this file defines functions only.

.animation_colors <- function() {
  c(
    cost = "#CF2120",
    downward_control = "#009BCE",
    upward_control = "#CE6600",
    drift_threshold = "#A3CF5B",
    intensity_threshold = "#8A2BE2",
    lower_barrier = "#CE6600",
    upper_barrier = "#009BCE",
    state = "black",
    neutral = "gray50"
  )
}

.mask_after <- function(v, i) {
  v <- as.numeric(v)
  if (i < length(v)) v[(i + 1):length(v)] <- NA_real_
  v
}

.frame_name <- function(prefix, digits, i) {
  sprintf("%s_%0*d", prefix, digits, i)
}

.animation_apply_par <- function(margins, axis_mgp = NULL, tcl = NULL,
                                 tick_cex = NULL, axis_title_cex = NULL,
                                 base_cex = NULL, mex = NULL) {
  args <- list(mar = margins)
  if (!is.null(axis_mgp)) args$mgp <- axis_mgp
  if (!is.null(tcl)) args$tcl <- tcl
  if (!is.null(tick_cex)) args$cex.axis <- tick_cex
  if (!is.null(axis_title_cex)) args$cex.lab <- axis_title_cex
  if (!is.null(base_cex)) args$cex <- base_cex
  if (!is.null(mex)) args$mex <- mex
  do.call(par, args)
  invisible()
}

.animation_axis_T <- function(side, t, labels = c("0", "T"), tick_cex = NULL) {
  args <- list(side = side, at = c(0, max(t, finite = TRUE)), labels = labels)
  if (!is.null(tick_cex)) args$cex.axis <- tick_cex
  do.call(axis, args)
}

.animation_control_ylim <- function(y) {
  ylim <- range(y, finite = TRUE)
  if (!all(is.finite(ylim))) ylim <- c(0, 1)
  ylim
}

# PNG-only saver. Keeping animations PNG-only avoids the extra device work done
# by render_and_save(), and keeps frame dimensions stable across a sequence.
render_and_save_png <- function(
    fname_base, plotfun,
    dir = "frames",
    width_in = 6.67, height_in = 4.67, dpi = 200,
    pointsize = NULL, family = NULL,
    match_current = TRUE,
    show = interactive(), save = TRUE
) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  if (match_current && dev.cur() != 1L) {
    sz <- dev.size("in")
    width_in <- sz[1]
    height_in <- sz[2]
    if (is.null(pointsize)) pointsize <- par("ps")
    if (is.null(family)) family <- par("family")
  }
  if (is.null(pointsize)) pointsize <- 12
  fam <- if (is.null(family) || !nzchar(family)) "sans" else family
  
  if (isTRUE(show)) {
    op <- par(family = fam)
    on.exit(par(op), add = TRUE)
    plotfun()
  }
  if (!isTRUE(save)) return(invisible())
  
  png(
    file.path(dir, paste0(fname_base, ".png")),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    pointsize = pointsize,
    type = getOption("bitmapType", "cairo"),
    antialias = "subpixel",
    bg = "transparent"
  )
  op <- par(family = fam)
  on.exit(par(op), add = TRUE)
  tryCatch(plotfun(), finally = dev.off())
  invisible()
}
