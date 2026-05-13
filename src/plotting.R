# Static plotting helpers for solver output and reflected paths.
#
# Code notation:
#   H, Hp  : marginal value derivative and its derivative.
#   X      : reflected state path.
#   U      : cumulative upward reflection at the lower barrier xL.
#   L      : cumulative downward reflection at the upper barrier xU.
#            Some figures label this process D_t to emphasize downward control.
#   xL, xU : lower and upper reflecting barriers.
#   xk, xl : drift- and intensity-ambiguity thresholds.
#   u, l   : upward and downward intervention costs.
#
# Loading this file defines functions only.

# ------------------------------- Shared helpers --------------------------------
.with_plot_margins <- function(expr) {
  op <- par(mar=c(2,2,0.5,0.5)+0.2)
  on.exit(par(op), add=TRUE)
  force(expr)
}

.threshold_style <- function() {
  list(
    cols = c(xL="#B22222", xk="#1E90FF", xl="#8A2BE2", xU="#228B22"),
    ltys = c(xL=1, xk=2, xl=2, xU=1)
  )
}

.sol_thresholds <- function(sol) {
  xs <- unlist(sol$x)[c("xL", "xk", "xl", "xU")]
  if (length(xs) != 4L || any(!is.finite(xs))) {
    stop("sol$x must contain finite xL, xk, xl, and xU values.", call. = FALSE)
  }
  xs
}

.expand_ylim <- function(y, top_blank=0.10, bottom_blank=0.05, fallback=c(-1, 1)) {
  ylim <- range(y, finite=TRUE)
  if (!all(is.finite(ylim))) ylim <- fallback
  if (diff(ylim) == 0) {
    bump <- max(1, abs(ylim[1]))
    ylim <- ylim + c(-1, 1) * 0.05 * bump
  }
  frac_band <- 1 - top_blank - bottom_blank
  ylim + diff(ylim)/frac_band * c(-bottom_blank, top_blank)
}

.legend_widths <- function(widths_in) {
  inch_to_user <- diff(grconvertX(c(0, 1), from="in", to="user"))
  widths_in * inch_to_user
}

.apply_plot_par <- function(margins=NULL, axis_mgp=NULL, tcl=NULL,
                            tick_cex=NULL, axis_title_cex=NULL,
                            base_cex=NULL, mex=NULL, xaxs=NULL,
                            oma=NULL) {
  args <- list()
  if (!is.null(margins)) args$mar <- margins
  if (!is.null(axis_mgp)) args$mgp <- axis_mgp
  if (!is.null(tcl)) args$tcl <- tcl
  if (!is.null(tick_cex)) args$cex.axis <- tick_cex
  if (!is.null(axis_title_cex)) args$cex.lab <- axis_title_cex
  if (!is.null(base_cex)) args$cex <- base_cex
  if (!is.null(mex)) args$mex <- mex
  if (!is.null(xaxs)) args$xaxs <- xaxs
  if (!is.null(oma)) args$oma <- oma
  if (length(args)) do.call(par, args)
  invisible()
}

.as_plot_path <- function(value, t, name) {
  value <- as.numeric(value)
  if (length(value) == 1L) return(rep(value, length(t)))
  if (length(value) != length(t)) {
    stop(sprintf("%s must be scalar or have the same length as the time grid.", name),
         call. = FALSE)
  }
  value
}

.thresholds_for_sim_plot <- function(sol, sim=NULL) {
  xs <- .sol_thresholds(sol)
  out <- list(xL=xs[["xL"]], xU=xs[["xU"]], xk=xs[["xk"]], xl=xs[["xl"]])
  if (!is.null(sim)) {
    for (nm in names(out)) {
      if (!is.null(sim[[nm]])) out[[nm]] <- sim[[nm]]
    }
  }
  out
}

render_and_save <- function(
    fname_base, plotfun, show=interactive(), save=TRUE,
    dir="figures",
    width_in = 6.67, height_in = 4.67, dpi = 300,
    pointsize = NULL, family = NULL,
    match_current = FALSE, use_cairo_pdf = TRUE
) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)

  # Match an existing device only when the caller explicitly asks for it.
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
    op <- par(family=fam)
    on.exit(par(op), add=TRUE)
    .with_plot_margins(plotfun())
  }

  if (!isTRUE(save)) return(invisible())

  png(file.path(dir, paste0(fname_base, ".png")),
      width=width_in, height=height_in, units="in",
      res=dpi, pointsize=pointsize,
      type=getOption("bitmapType", "cairo"), antialias="subpixel")
  tryCatch(.with_plot_margins(plotfun()), finally=dev.off())

  if (isTRUE(use_cairo_pdf) && capabilities("cairo")) {
    cairo_pdf(file.path(dir, paste0(fname_base, ".pdf")),
              width=width_in, height=height_in,
              pointsize=pointsize, family=fam)
  } else {
    pdf(file.path(dir, paste0(fname_base, ".pdf")),
        width=width_in, height=height_in,
        pointsize=pointsize, family=fam, useDingbats=FALSE)
  }
  tryCatch(.with_plot_margins(plotfun()), finally=dev.off())
  invisible()
}

# ----------------------------- H and H' figures --------------------------------
plot_H <- function(sol, params, show=interactive(), save=TRUE, n = 600,
                   pad_x_frac = 0.10, top_blank = 0.10, bottom_blank = 0.05,
                   out_base="H_thresholds", dir="figures",
                   width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                   margins = c(3.2, 3.25, 1.05, 0.9),
                   axis_mgp = c(2.2, 0.7, 0),
                   tcl = -0.25,
                   show_x_axis_title = FALSE,
                   show_y_axis_title = TRUE,
                   x_axis_title = expression(x),
                   y_axis_title = expression(H(x) == V*minute*(x)),
                   show_tick_labels = TRUE,
                   tick_cex = 1.7,
                   axis_title_cex = 1,
                   base_cex = 1, mex = 1) {

  xs <- .sol_thresholds(sol)
  style <- .threshold_style()

  x_max <- max(xs[["xU"]], xs[["xl"]])
  W <- x_max - xs[["xL"]]
  x_from <- xs[["xL"]] - pad_x_frac*W
  x_to <- x_max + pad_x_frac*W
  ylim_y <- .expand_ylim(sol$H(seq(x_from, x_to, length.out=n)),
                         top_blank=top_blank, bottom_blank=bottom_blank)

  plotfun <- function() {
    op <- par(no.readonly=TRUE)
    on.exit(par(op), add=TRUE)

    .apply_plot_par(
      margins=margins, axis_mgp=axis_mgp, tcl=tcl,
      tick_cex=tick_cex, axis_title_cex=axis_title_cex,
      base_cex=base_cex, mex=mex, xaxs="i", oma=c(0,0,0,0)
    )

    curve(sol$H(x), from=x_from, to=x_to, n=n,
          xlab=if (show_x_axis_title) x_axis_title else "",
          ylab=if (show_y_axis_title) y_axis_title else "",
          col="black", lwd=1.4, ylim=ylim_y,
          xaxt="n", yaxt="n", bty="n")

    axis(1, labels=show_tick_labels)
    axis(2, labels=show_tick_labels)

    segments(x0=xs, y0=rep(-params$u, length(xs)), x1=xs, y1=rep(params$l, length(xs)),
             lty=style$ltys, col=style$cols, lwd=2)
    abline(h=-params$u, lty=3, col="#777777")
    abline(h=0, lty=3, col="#BBBBBB")
    abline(h=params$l, lty=4, col="#777777")

    tw <- .legend_widths(c(0.36, rep(0.23, 5), 0.38))
    legend("topleft",
           legend=expression(H(x), underline(x), x^kappa, x^lambda, bar(x), -c[U], c[D]),
           lty=c(1, 1, 2, 2, 1, 3, 4),
           col=c("black", style$cols, "#777777", "#777777"),
           lwd=c(1.4, rep(2, 6)),
           bty="n", horiz=TRUE, seg.len=1.5, x.intersp=0.5,
           text.width=tw)

    box()
  }

  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in=width_in, height_in=height_in, dpi=dpi,
                  match_current=match_current)
}

plot_H_prime <- function(sol, params, show=interactive(), save=TRUE, n = 800,
                         pad_x_frac = 0.10, top_blank = 0.10, bottom_blank = 0.05,
                         out_base="Hprime_thresholds", dir="figures",
                         width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                         margins = NULL,
                         axis_mgp = NULL,
                         tcl = NULL,
                         show_x_axis_title = FALSE,
                         show_y_axis_title = TRUE,
                         x_axis_title = expression(x),
                         y_axis_title = expression(H*minute*(x) == W(x)),
                         show_tick_labels = TRUE,
                         tick_cex = NULL,
                         axis_title_cex = 1,
                         base_cex = NULL, mex = NULL) {
  stopifnot(!is.null(sol$Hp))

  xs <- .sol_thresholds(sol)
  style <- .threshold_style()

  x_max <- max(xs[["xU"]], xs[["xl"]])
  W <- x_max - xs[["xL"]]
  x_from <- xs[["xL"]] - pad_x_frac*W
  x_to <- x_max + pad_x_frac*W

  # Break the derivative line around thresholds to avoid connecting jumps.
  xg <- seq(x_from, x_to, length.out=n)
  yg <- sol$Hp(xg)
  for (b in xs) {
    j <- which.min(abs(xg - b))
    for (k in c(j-1, j, j+1)) {
      if (k >= 1 && k <= length(yg)) yg[k] <- NA_real_
    }
  }
  ylim_y <- .expand_ylim(yg, top_blank=top_blank, bottom_blank=bottom_blank)

  plotfun <- function() {
    op <- par(no.readonly=TRUE)
    on.exit(par(op), add=TRUE)

    .apply_plot_par(
      margins=margins, axis_mgp=axis_mgp, tcl=tcl,
      tick_cex=tick_cex, axis_title_cex=axis_title_cex,
      base_cex=base_cex, mex=mex, xaxs="i"
    )

    plot(xg, yg, type="l",
         xlab=if (show_x_axis_title) x_axis_title else "",
         ylab=if (show_y_axis_title) y_axis_title else "",
         col="black", lwd=1.4, ylim=ylim_y,
         xaxt="n", yaxt="n")
    axis(1, labels=show_tick_labels)
    axis(2, labels=show_tick_labels)

    segments(x0=xs, y0=rep(ylim_y[1], length(xs)),
             x1=xs, y1=rep(ylim_y[2], length(xs)),
             lty=style$ltys, col=style$cols, lwd=2)
    abline(h=0, lty=3, col="#BBBBBB")

    tw <- .legend_widths(c(0.36, rep(0.23, 5), 0.20))
    legend("topleft",
           legend=expression(H*minute*(x), underline(x), x^kappa, x^lambda, bar(x), 0),
           lty=c(1, 1, 2, 2, 1, 3),
           col=c("black", style$cols, "#BBBBBB"),
           lwd=c(1.4, rep(2, 4), 1.2),
           bty="n", horiz=TRUE, seg.len=1.5, x.intersp=0.5,
           text.width=tw)

    box()
  }

  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in=width_in, height_in=height_in, dpi=dpi,
                  match_current=match_current)
}

# ----------------------------- Reflected paths ---------------------------------
.draw_threshold_line <- function(t, y, col, lty, lwd=2) {
  if (length(unique(y)) == 1L) {
    abline(h=y[1], col=col, lty=lty, lwd=lwd)
  } else {
    lines(t, y, col=col, lty=lty, lwd=lwd)
  }
}

.draw_reflected_panel <- function(t, X, xL, xU, xk, xl,
                                  top_blank=0.075, bottom_blank=0.05,
                                  draw_legend=TRUE,
                                  show_x_axis=FALSE,
                                  show_y_axis=TRUE,
                                  show_tick_labels=TRUE) {
  xL <- .as_plot_path(xL, t, "xL")
  xU <- .as_plot_path(xU, t, "xU")
  xk <- .as_plot_path(xk, t, "xk")
  xl <- .as_plot_path(xl, t, "xl")
  if (any(xL >= xU)) stop("xL must be below xU for the whole plotted path.", call. = FALSE)

  band_low <- min(xL, finite=TRUE)
  band_high <- max(xU, finite=TRUE)
  Hcorr <- band_high - band_low
  frac <- 1 - top_blank - bottom_blank
  Htot <- Hcorr/frac
  ylim <- c(band_low - Htot*bottom_blank, band_high + Htot*top_blank)

  style <- .threshold_style()

  plot(t, X, type="n", xlab="", ylab="", ylim=ylim,
       xaxt="n", yaxt="n")
  usr <- par("usr")
  if (length(unique(xL)) == 1L && length(unique(xU)) == 1L) {
    rect(usr[1], xL[1], usr[2], xU[1], col=adjustcolor("gray85", 0.6), border=NA)
  } else {
    polygon(c(t, rev(t)), c(xL, rev(xU)), col=adjustcolor("gray85", 0.6), border=NA)
  }

  lines(t, X, lwd=1.2)
  .draw_threshold_line(t, xL, col=style$cols[["xL"]], lty=style$ltys[["xL"]])
  .draw_threshold_line(t, xU, col=style$cols[["xU"]], lty=style$ltys[["xU"]])
  .draw_threshold_line(t, xk, col=style$cols[["xk"]], lty=style$ltys[["xk"]])
  .draw_threshold_line(t, xl, col=style$cols[["xl"]], lty=style$ltys[["xl"]])

  if (show_x_axis) axis(1, labels=show_tick_labels)
  if (show_y_axis) axis(2, labels=show_tick_labels)

  if (draw_legend) {
    tw <- .legend_widths(c(0.36, rep(0.23, 4)))
    legend("topleft",
           legend=expression(bar(X)[t], underline(x), x^kappa, x^lambda, bar(x)),
           lty=c(1, 1, 2, 2, 1),
           col=c("black", style$cols),
           lwd=c(1.2, rep(2, 4)),
           bty="n", horiz=TRUE, x.intersp=0.5, seg.len=2,
           text.width=tw)
  }
  box()
  invisible(ylim)
}

plot_reflected_jd <- function(sol, params, sim=NULL, seed=123, show=interactive(),
                              save=TRUE, out_base="reflected_jd", dir="figures",
                              width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                              top_blank=0.075, bottom_blank=0.05,
                              margins = NULL,
                              axis_mgp = NULL,
                              tcl = NULL,
                              show_x_axis_title = TRUE,
                              show_y_axis_title = FALSE,
                              x_axis_title = "time",
                              y_axis_title = expression(bar(X)[t]),
                              show_tick_labels = TRUE,
                              tick_cex = NULL,
                              axis_title_cex = 1,
                              base_cex = NULL, mex = NULL) {
  if (is.null(sim)) { set.seed(seed); sim <- simulate_reflected_jd(params=params, thresholds=sol$x) }
  thresholds <- .thresholds_for_sim_plot(sol, sim)
  t <- sim$time
  X <- sim$X

  plotfun <- function() {
    op <- par(no.readonly=TRUE)
    on.exit(par(op), add=TRUE)
    .apply_plot_par(
      margins=margins, axis_mgp=axis_mgp, tcl=tcl,
      tick_cex=tick_cex, axis_title_cex=axis_title_cex,
      base_cex=base_cex, mex=mex
    )
    .draw_reflected_panel(t, X,
                          thresholds$xL, thresholds$xU, thresholds$xk, thresholds$xl,
                          top_blank=top_blank, bottom_blank=bottom_blank,
                          show_x_axis=TRUE, show_y_axis=TRUE,
                          show_tick_labels=show_tick_labels)
    if (show_x_axis_title || show_y_axis_title) {
      title(xlab=if (show_x_axis_title) x_axis_title else "",
            ylab=if (show_y_axis_title) y_axis_title else "",
            cex.lab=axis_title_cex)
    }
  }

  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in=width_in, height_in=height_in, dpi=dpi,
                  match_current=match_current)
}

plot_controls <- function(sim, show=interactive(), save=TRUE,
                          out_base="singular_controls", dir="figures",
                          width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                          margins = NULL,
                          axis_mgp = NULL,
                          tcl = NULL,
                          show_x_axis_title = TRUE,
                          show_y_axis_title = TRUE,
                          x_axis_title = "time",
                          y_axis_title = "cumulative push",
                          show_tick_labels = TRUE,
                          tick_cex = NULL,
                          axis_title_cex = 1,
                          base_cex = NULL, mex = NULL) {
  plotfun <- function() {
    op <- par(no.readonly=TRUE)
    on.exit(par(op), add=TRUE)
    .apply_plot_par(
      margins=margins, axis_mgp=axis_mgp, tcl=tcl,
      tick_cex=tick_cex, axis_title_cex=axis_title_cex,
      base_cex=base_cex, mex=mex
    )

    rng <- range(sim$U, sim$L)
    plot(sim$time, sim$L, type="s",
         xlab=if (show_x_axis_title) x_axis_title else "",
         ylab=if (show_y_axis_title) y_axis_title else "",
         lwd=1.6, ylim=rng, lty=2,
         xaxt=if (show_tick_labels) "s" else "n",
         yaxt=if (show_tick_labels) "s" else "n")
    lines(sim$time, sim$U, type="s", lwd=1.6)
    legend("topleft",
           legend=c(expression(D[t]~"(pushes down)"), expression(U[t]~"(pushes up)")),
           lty=c(2, 1), lwd=2, bty="n")
    box()
  }

  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in=width_in, height_in=height_in, dpi=dpi,
                  match_current=match_current)
}

plot_reflected_with_controls <- function(sol, params, sim=NULL, seed=123,
                                         show=interactive(), save=TRUE,
                                         out_base="reflected_with_controls",
                                         dir="figures", draw_legend = TRUE,
                                         width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                                         heights = c(0.7, 1.7, 0.7),
                                         top_blank=0.075, bottom_blank=0.05,
                                         margins = c(3.2, 3.25, 1.05, 0.9),
                                         axis_mgp = c(2.2, 0.7, 0),
                                         tcl = -0.25,
                                         show_x_axis_title = TRUE,
                                         show_y_axis_title = TRUE,
                                         x_axis_title = "time",
                                         y_axis_title = list(top = expression(L[t]),
                                                             bottom = expression(U[t])),
                                         show_tick_labels = TRUE,
                                         tick_cex = 1.7,
                                         axis_title_cex = 1,
                                         base_cex = 1, mex = 1) {
  if (is.null(sim)) { set.seed(seed); sim <- simulate_reflected_jd(params=params, thresholds=sol$x) }
  thresholds <- .thresholds_for_sim_plot(sol, sim)
  t <- sim$time
  X <- sim$X

  if (is.list(y_axis_title)) {
    y_axis_title_top <- y_axis_title[["top"]]
    y_axis_title_bottom <- y_axis_title[["bottom"]]
    if (is.null(y_axis_title_top) && length(y_axis_title) >= 1) {
      y_axis_title_top <- y_axis_title[[1]]
    }
    if (is.null(y_axis_title_bottom)) {
      y_axis_title_bottom <- y_axis_title_top
    }
  } else {
    y_axis_title_top <- y_axis_title
    y_axis_title_bottom <- y_axis_title
  }

  bottom <- margins[1]
  left <- margins[2]
  top <- margins[3]
  right <- margins[4]

  plotfun <- function() {
    op <- par(no.readonly=TRUE)
    on.exit(par(op), add=TRUE)

    par(oma=c(0,0,0,0), cex=base_cex, mex=mex)
    layout(matrix(1:3, nrow=3), heights=heights)

    # Downward reflection at the upper barrier.
    par(mar=c(0.4, left, top, right),
        mgp=axis_mgp, tcl=tcl,
        cex=base_cex, mex=mex,
        cex.axis=tick_cex, cex.lab=axis_title_cex)
    plot(t, sim$L, type="s",
         xlab="", ylab=if (show_y_axis_title) y_axis_title_top else "",
         xaxt="n", yaxt="n",
         lwd=2, col="#228B22")
    axis(2, labels=show_tick_labels)
    legend("topleft", legend=expression(D[t]), cex=1,
           lty=1, lwd=2, col="#228B22", bty="n")
    box()

    # Reflected state path with barriers and ambiguity thresholds.
    par(mar=c(0.2, left, 0.2, right),
        mgp=axis_mgp, tcl=tcl,
        cex=base_cex, mex=mex,
        cex.axis=tick_cex, cex.lab=axis_title_cex)
    .draw_reflected_panel(t, X,
                          thresholds$xL, thresholds$xU, thresholds$xk, thresholds$xl,
                          top_blank=top_blank, bottom_blank=bottom_blank,
                          draw_legend=draw_legend,
                          show_x_axis=FALSE, show_y_axis=TRUE,
                          show_tick_labels=show_tick_labels)

    # Upward reflection at the lower barrier.
    par(mar=c(bottom, left, 0.4, right),
        mgp=axis_mgp, tcl=tcl,
        cex=base_cex, mex=mex,
        cex.axis=tick_cex, cex.lab=axis_title_cex)
    plot(t, sim$U, type="s",
         xlab=if (show_x_axis_title) x_axis_title else "",
         ylab=if (show_y_axis_title) y_axis_title_bottom else "",
         xaxt="n", yaxt="n",
         lwd=2, col="#B22222")
    axis(1, labels=show_tick_labels)
    axis(2, labels=show_tick_labels)
    legend("topleft", legend=expression(U[t]), cex=1,
           lty=1, lwd=2, col="#B22222", bty="n")
    box()
  }

  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in=width_in, height_in=height_in, dpi=dpi,
                  match_current=match_current)
}
