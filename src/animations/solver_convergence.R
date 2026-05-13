# Solver-convergence animation frames.
#
# Loading this file defines functions only.

# Precompute fixed axes so the curve and gamma panels do not jump between
# frames.
solver_anim_limits <- function(history, params,
                               n_grid = 700,
                               pad_x_frac = 0.10,
                               top_blank = 0.10, bottom_blank = 0.05) {
  stopifnot(nrow(history) >= 1)
  
  history$xL <- as.numeric(history$xL)
  history$xU <- as.numeric(history$xU)
  history$xl <- as.numeric(history$xl)
  
  xL_min <- min(history$xL, na.rm=TRUE)
  x_max  <- max(pmax(history$xU, history$xl), na.rm=TRUE)
  
  if (!is.finite(xL_min) || !is.finite(x_max)) {
    stop("Non-finite x-range in iter_history. Check that iter_history has finite xL/xU/xl values.")
  }
  
  W <- x_max - xL_min
  if (!is.finite(W) || W <= 0) W <- 1
  
  x_from <- xL_min - pad_x_frac * W
  x_to   <- x_max  + pad_x_frac * W
  
  xg <- seq(x_from, x_to, length.out = n_grid)
  
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
                                          bty="n", box_lwd=1,
                                          tick_cex = NULL,
                                          h_tick_cex = NULL,
                                          axis_title_cex = NULL,
                                          text_cex = 1.5,
                                          legend_cex = 1.5,
                                          point_cex = 1.5,
                                          legend_point_cex = 1.2,
                                          axis_mgp = NULL,
                                          tcl = NULL,
                                          base_cex = NULL,
                                          mex = NULL) {
  n <- nrow(history)
  stopifnot(i >= 1, i <= n)
  
  fname <- .frame_name(base, digits, i)
  h_tick_cex <- h_tick_cex %||% tick_cex %||% 1.5
  
  plotfun <- function() {
    op <- par(no.readonly=TRUE); on.exit(par(op), add=TRUE)
    
    # Layout:
    #   row 1: gamma (left) + text (right)
    #   row 2: H     (left) + legend (right)
    layout(matrix(c(1,2,
                    3,4), nrow=2, byrow=TRUE),
           widths=c(0.78, 0.22), heights=c(0.42, 0.58))
    
    .animation_apply_par(
      margins = c(3.2, 3.4, 0.6, 0.3),
      axis_mgp = axis_mgp,
      tcl = tcl,
      tick_cex = tick_cex,
      axis_title_cex = axis_title_cex,
      base_cex = base_cex,
      mex = mex
    )
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
    points(i, gy[i], pch=16, cex=point_cex)
    
    if (!identical(bty, "n")) box(bty=bty, lwd=box_lwd)
    
    .animation_apply_par(
      margins = c(3.2, 0.6, 0.4, 0.6),
      axis_mgp = axis_mgp,
      tcl = tcl,
      axis_title_cex = axis_title_cex,
      base_cex = base_cex,
      mex = mex
    )
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    
    x_shift    <- 0.0
    y_shift    <- -0.1
    line_space <- 2.5
    
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
        text(x0, y, labels=e, adj=c(0,1), xpd=NA, cex=text_cex)
        y <- y - dy
      }
    }
    
    sol_i <- build_suboptimal_H(params,
                                history$xL[i], history$xk[i], history$xl[i], history$xU[i])
    
    xs <- sol_i$x
    xg <- lims$xg
    yg <- sol_i$H(xg)
    
    .animation_apply_par(
      margins = c(3.2, 3.4, 0.6, 0.3),
      axis_mgp = axis_mgp,
      tcl = tcl,
      tick_cex = h_tick_cex,
      axis_title_cex = axis_title_cex,
      base_cex = base_cex,
      mex = mex
    )
    plot(xg, yg, type="l", lwd=2,
         xlab=if (isTRUE(show_axes)) "x" else "",
         ylab=if (isTRUE(show_axes)) expression(H(x)) else "",
         xaxt=if (isTRUE(show_axes)) "s" else "n",
         yaxt=if (isTRUE(show_axes)) "s" else "n",
         xlim=c(lims$x_from, lims$x_to), ylim=lims$ylim_H, bty="n")
    
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
    
    .animation_apply_par(
      margins = c(0.8, 0.6, 0.6, 0.6),
      axis_mgp = axis_mgp,
      tcl = tcl,
      base_cex = base_cex,
      mex = mex
    )
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    
    cols_thr <- c("#B22222","#1E90FF","#8A2BE2","#228B22")
    ltys_thr <- c(1,1,1,1)
    
    legend("center", bty="n", xpd=NA, cex=legend_cex, y.intersp=1.05,
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
           pt.cex = c(NA, NA, NA, NA, NA, NA, NA, legend_point_cex))
    
    
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
                                                     bty="n", box_lwd=1,
                                                     tick_cex = NULL,
                                                     h_tick_cex = NULL,
                                                     axis_title_cex = NULL,
                                                     text_cex = 1.5,
                                                     legend_cex = 1.5,
                                                     point_cex = 1.5,
                                                     legend_point_cex = 1.2,
                                                     axis_mgp = NULL,
                                                     tcl = NULL,
                                                     base_cex = NULL,
                                                     mex = NULL,
                                                     frame_pause = 0.5) {
  history <- sol_opt$iter_history
  
  if (is.matrix(history)) history <- as.data.frame(history, stringsAsFactors = FALSE)
  
  if (is.null(history)) stop("sol_opt$iter_history is NULL. Run solve_optimal_barriers(..., record_iterates=TRUE).")
  if (!is.data.frame(history)) stop(sprintf("iter_history must be a data.frame; got class: %s", paste(class(history), collapse=", ")))
  
  req <- c("xL","xk","xl","xU","gamma","fmax","stage","Hp_xL","dHp_xk","dHp_xl","Hp_xU")
  miss <- setdiff(req, names(history))
  if (length(miss)) stop("iter_history is missing columns: ", paste(miss, collapse=", "))
  
  num_cols <- setdiff(req, "stage")
  history[num_cols] <- lapply(history[num_cols], function(v) suppressWarnings(as.numeric(v)))
  
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
                                  show_axes=show_axes, bty=bty, box_lwd=box_lwd,
                                  tick_cex=tick_cex, h_tick_cex=h_tick_cex,
                                  axis_title_cex=axis_title_cex,
                                  text_cex=text_cex, legend_cex=legend_cex,
                                  point_cex=point_cex,
                                  legend_point_cex=legend_point_cex,
                                  axis_mgp=axis_mgp, tcl=tcl,
                                  base_cex=base_cex, mex=mex)
    if (isTRUE(show) && is.finite(frame_pause) && frame_pause > 0) {
      Sys.sleep(frame_pause)
    }
    }
  invisible(lims)
}
