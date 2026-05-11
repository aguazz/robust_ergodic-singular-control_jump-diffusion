# Static plotting helpers for solver output and reflected paths.
# Source split from legacy/functions_solver_and_images.R.

# ================================ Figures =======================================
.with_plot_margins <- function(expr) { 
  op <- par(mar=c(2,2,0.5,0.5)+0.2)
  on.exit(par(op), add=TRUE); force(expr)
}

render_and_save <- function(
    fname_base, plotfun, show=interactive(), save=TRUE,
    dir="figures",
    width_in = 6.67, height_in = 4.67, dpi = 300,
    pointsize = NULL, family = NULL,
    match_current = FALSE, use_cairo_pdf = TRUE
) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  
  # Derive device settings (size, ps, family) from current device when available
  if (match_current && dev.cur() != 1L) {
    sz <- dev.size("in")
    width_in  <- sz[1]; height_in <- sz[2]
    if (is.null(pointsize)) pointsize <- par("ps")
    if (is.null(family))    family    <- par("family")
  }
  if (is.null(pointsize)) pointsize <- 12
  fam <- if (is.null(family) || !nzchar(family)) "sans" else family
  
  # SHOW: draw on the current device but force the same family for consistency
  if (isTRUE(show)) {
    op <- par(family = fam)               # temporary — reverts on exit
    on.exit(par(op), add = TRUE)
    .with_plot_margins(plotfun())
  }
  
  if (!isTRUE(save)) return(invisible())
  
  # PNG: exact physical size + pointsize
  png(file.path(dir, paste0(fname_base, ".png")),
      width = width_in, height = height_in, units = "in",
      res = dpi, pointsize = pointsize,
      type = getOption("bitmapType", "cairo"), antialias = "subpixel")
  .with_plot_margins(plotfun()); dev.off()
  
  # PDF: same physical size + pointsize + family
  if (isTRUE(use_cairo_pdf) && capabilities("cairo")) {
    cairo_pdf(file.path(dir, paste0(fname_base, ".pdf")),
              width = width_in, height = height_in,
              pointsize = pointsize, family = fam)
  } else {
    pdf(file.path(dir, paste0(fname_base, ".pdf")),
        width = width_in, height = height_in,
        pointsize = pointsize, family = fam, useDingbats = FALSE)
  }
  .with_plot_margins(plotfun()); dev.off()
}

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
  
  xs   <- unlist(sol$x)[c("xL","xk","xl","xU")]
  xL <- xs[1]; xk <- xs[2]; xl <- xs[3]; xU <- xs[4]
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22"); ltys <- c(1,2,2,1)
  
  x_max <- max(xU, xl)
  W <- x_max - xL
  x_from <- xL - pad_x_frac*W
  x_to   <- x_max + pad_x_frac*W
  
  frac_band <- 1 - top_blank - bottom_blank
  ylim_y <- range(sol$H(seq(x_from, x_to, l = n)))
  ylim_y <- ylim_y + diff(ylim_y)/frac_band * c(-bottom_blank, top_blank)
  
  plotfun <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add=TRUE)
    
    par(xaxs="i",
        oma = c(0,0,0,0),
        cex = base_cex,
        mex = mex,
        mar = margins,
        mgp = axis_mgp,
        tcl = tcl,
        cex.axis = tick_cex,
        cex.lab = axis_title_cex)
    
    curve(sol$H(x), from=x_from, to=x_to, n=n,
          xlab=if (show_x_axis_title) x_axis_title else "",
          ylab=if (show_y_axis_title) y_axis_title else "",
          col="black", lwd=1.4, ylim=ylim_y,
          xaxt="n", yaxt="n", bty="n")
    
    axis(1, labels = show_tick_labels)
    axis(2, labels = show_tick_labels)
    
    segments(x0=xs, y0=rep(-params$u, length(xs)), x1=xs, y1=rep(params$l, length(xs)),
             lty=ltys, col=cols, lwd=2)
    abline(h=-params$u, lty=3, col="#777777")
    abline(h=0,       lty=3, col="#BBBBBB")
    abline(h=params$l, lty=4, col="#777777")
    
    inch_to_user <- diff(grconvertX(c(0,1), from="in", to="user"))
    col_w <- 0.23 * inch_to_user
    tw <- c(0.36 * inch_to_user, rep(col_w, 5), 0.38 * inch_to_user)
    
    legend("topleft",
           legend=expression(H(x), underline(x), x^kappa, x^lambda, bar(x), -c[U], c[D]),
           lty=c(1,1,2,2,1,3,4),
           col=c("black", cols, "#777777", "#777777"),
           lwd=c(1.4, rep(2,6)),
           bty="n", horiz=TRUE, seg.len=1.5, x.intersp=0.5,
           text.width=tw)
    
    box()
  }
  
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}


plot_H_prime <- function(sol, params, show=interactive(), save=TRUE, n = 800,
                         pad_x_frac = 0.10, top_blank = 0.10, bottom_blank = 0.05,
                         out_base="Hprime_thresholds", dir="figures",
                         width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                         show_x_axis_title = FALSE,
                         show_y_axis_title = TRUE,
                         x_axis_title = expression(x),
                         y_axis_title = expression(H*minute*(x) == W(x)),
                         axis_title_cex = 1) {
  stopifnot(!is.null(sol$Hp))
  
  xs   <- unlist(sol$x)[c("xL","xk","xl","xU")]
  xL <- xs[1]; xk <- xs[2]; xl <- xs[3]; xU <- xs[4]
  
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22")  # xL,xk,xl,xU
  ltys <- c(1,2,2,1)
  
  x_max <- max(xU, xl)
  W <- x_max - xL
  x_from <- xL - pad_x_frac*W
  x_to   <- x_max + pad_x_frac*W
  
  # sample grid and break the line at thresholds to avoid artificial connections
  xg <- seq(x_from, x_to, length.out = n)
  yg <- sol$Hp(xg)
  
  # Force line breaks around the 4 special points (xL,xk,xl,xU)
  for (b in xs) {
    j <- which.min(abs(xg - b))
    for (k in c(j-1, j, j+1)) {
      if (k >= 1 && k <= length(yg)) yg[k] <- NA_real_
    }
  }
  
  # y-limits with controlled blank margins (like plot_H)
  frac_band <- 1 - top_blank - bottom_blank
  ylim_y <- range(yg, finite = TRUE)
  if (!all(is.finite(ylim_y))) ylim_y <- c(-1, 1)
  if (diff(ylim_y) == 0) {
    bump <- max(1, abs(ylim_y[1]))
    ylim_y <- ylim_y + c(-1, 1) * 0.05 * bump
  }
  ylim_y <- ylim_y + diff(ylim_y)/frac_band * c(-bottom_blank, top_blank)
  
  plotfun <- function() {
    op <- par(xaxs="i", cex.lab = axis_title_cex); on.exit(par(op), add=TRUE)
    
    plot(xg, yg, type="l",
         xlab=if (show_x_axis_title) x_axis_title else "",
         ylab=if (show_y_axis_title) y_axis_title else "",
         col="black", lwd=1.4, ylim=ylim_y)
    
    # thresholds as vertical segments spanning the panel
    segments(x0=xs, y0=rep(ylim_y[1], length(xs)),
             x1=xs, y1=rep(ylim_y[2], length(xs)),
             lty=ltys, col=cols, lwd=2)
    
    # reference line at 0
    abline(h=0, lty=3, col="#BBBBBB")
    
    # legend (same style as plot_H)
    inch_to_user <- diff(grconvertX(c(0,1), from="in", to="user"))
    col_w <- 0.23 * inch_to_user
    tw <- c(0.36 * inch_to_user,
            rep(col_w, 5),
            0.20 * inch_to_user)
    
    legend("topleft",
           legend = expression(H*minute*(x), underline(x), x^kappa, x^lambda, bar(x), 0),
           lty    = c(1, 1, 2, 2, 1, 3),
           col    = c("black", cols, "#BBBBBB"),
           lwd    = c(1.4, rep(2,4), 1.2),
           bty="n", horiz=TRUE, seg.len=1.5, x.intersp=0.5,
           text.width=tw)
    
    box()
  }
  
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}


# --------------------------- Small reusable drawer for the middle panel --------
# Draw the reflected trajectory with thresholds; returns the ylim used (so panels
# can be tuned consistently if needed).
.draw_reflected_panel <- function(t, sim, xL, xU, xk, xl,
                                  top_blank=0.075, bottom_blank=0.05,
                                  draw_legend=TRUE) {
  Hcorr <- xU - xL
  frac  <- 1 - top_blank - bottom_blank
  Htot  <- Hcorr/frac
  ylim  <- c(xL - Htot*bottom_blank, xU + Htot*top_blank)
  
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22") # xL, xk, xl, xU
  ltys <- c(1,2,2,1)
  
  plot(t, sim, type="n", xlab="", ylab="", ylim=ylim, xaxt="n")
  usr <- par("usr")
  rect(usr[1], xL, usr[2], xU, col=adjustcolor("gray85", 0.6), border=NA)
  
  lines(t, sim, lwd=1.2)
  abline(h=c(xL, xU), col=c(cols[1], cols[4]), lty=c(ltys[1], ltys[4]), lwd=2)
  abline(h=c(xk, xl), col=c(cols[2], cols[3]), lty=c(ltys[2], ltys[3]), lwd=2)
  
  if (draw_legend) {
    inch_to_user <- diff(grconvertX(c(0,1), from="in", to="user"))
    col_w <- 0.23 * inch_to_user
    tw <- c(0.36 * inch_to_user, rep(col_w, 4))  # 5 entries total
    
    legend("topleft",
           legend = expression(bar(X)[t], underline(x), x^kappa, x^lambda, bar(x)),
           lty    = c(1, 1, 2, 2, 1),                # 5 entries
           col    = c("black", cols),                # 1 + 4 entries
           lwd    = c(1.2, rep(2,4)),
           bty    = "n", horiz = TRUE, x.intersp = 0.5, seg.len = 2,
           text.width = tw
    )
  }
  box()
  invisible(ylim)
}

# --------------------------- Original single plots (kept, just call helper) ----
plot_reflected_jd <- function(sol, params, sim=NULL, seed=123, show=interactive(), 
                              save=TRUE, out_base="reflected_jd", dir="figures",
                              width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                              top_blank=0.075, bottom_blank=0.05,
                              show_x_axis_title = TRUE,
                              show_y_axis_title = FALSE,
                              x_axis_title = "time",
                              y_axis_title = expression(bar(X)[t]),
                              axis_title_cex = 1) {
  xs <- sol$x
  if (is.null(sim)) { set.seed(seed); sim <- simulate_reflected_jd(params=params, thresholds=xs) }
  t <- sim$time; X <- sim$X; xL <- xs$xL; xU <- xs$xU; xk <- xs$xk; xl <- xs$xl
  plotfun <- function() {
    .draw_reflected_panel(t, X, xL, xU, xk, xl, top_blank=top_blank, bottom_blank=bottom_blank)
    axis(1)
    if (show_x_axis_title || show_y_axis_title) {
      title(xlab = if (show_x_axis_title) x_axis_title else "",
            ylab = if (show_y_axis_title) y_axis_title else "",
            cex.lab = axis_title_cex)
    }
  }
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_controls <- function(sim, show=interactive(), save=TRUE,
                          out_base="singular_controls", dir="figures",
                          width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                          axis_title_cex = 1) {
  plotfun <- function() {
    rng <- range(sim$U, sim$L)
    plot(sim$time, sim$L, type="s", xlab="time", ylab="cumulative push",
         lwd=1.6, ylim=rng, lty=2, cex.lab = axis_title_cex)
    lines(sim$time, sim$U, type="s", lwd=1.6)
    legend("topleft",
           legend=c(expression(D[t]~"(pushes down)"), expression(U[t]~"(pushes up)")),
           lty=c(2,1), lwd=2, bty="n")
    box()
  }
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

# --------------------------- NEW: merged 3-panel figure -------------------------
# Top:   L_t (pushes down)
# Middle: reflected process + thresholds
# Bottom: U_t (pushes up)
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
  xs <- sol$x
  if (is.null(sim)) { set.seed(seed); sim <- simulate_reflected_jd(params=params, thresholds=xs) }
  t <- sim$time; X <- sim$X; xL <- xs$xL; xU <- xs$xU; xk <- xs$xk; xl <- xs$xl
  
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
  
  bottom <- margins[1]; left <- margins[2]; top <- margins[3]; right <- margins[4]
  
  plotfun <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add=TRUE)
    
    par(oma = c(0,0,0,0), cex = base_cex, mex = mex)
    
    layout(matrix(1:3, nrow=3), heights = heights)
    
    # --- TOP: D_t
    par(mar = c(0.4, left, top, right),
        mgp=axis_mgp, tcl=tcl,
        cex = base_cex, mex = mex,
        cex.axis=tick_cex, cex.lab = axis_title_cex)
    plot(t, sim$L, type="s",
         xlab="", ylab=if (show_y_axis_title) y_axis_title_top else "",
         xaxt="n", yaxt="n",
         lwd=2, col = "#228B22")
    axis(2, labels = show_tick_labels)
    legend("topleft", legend = expression(D[t]), cex = 1,
           lty = 1, lwd = 2, col = "#228B22", bty = "n")
    box()
    
    # --- MIDDLE
    par(mar = c(0.2, left, 0.2, right),
        mgp=axis_mgp, tcl=tcl,
        cex = base_cex, mex = mex,
        cex.axis=tick_cex, cex.lab = axis_title_cex)
    .draw_reflected_panel(t, X, xL, xU, xk, xl,
                          top_blank=top_blank, bottom_blank=bottom_blank,
                          draw_legend=draw_legend)
    
    # --- BOTTOM: U_t
    par(mar = c(bottom, left, 0.4, right),
        mgp=axis_mgp, tcl=tcl,
        cex = base_cex, mex = mex,
        cex.axis=tick_cex, cex.lab = axis_title_cex)
    plot(t, sim$U, type="s",
         xlab=if (show_x_axis_title) x_axis_title else "",
         ylab=if (show_y_axis_title) y_axis_title_bottom else "",
         xaxt="n", yaxt="n",
         lwd=2, col = "#B22222")
    axis(1, labels = show_tick_labels)
    axis(2, labels = show_tick_labels)
    legend("topleft", legend = expression(U[t]), cex = 1,
           lty = 1, lwd = 2, col = "#B22222", bty = "n")
    box()
  }
  
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}


