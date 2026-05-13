# Comparative-sweep experiment helpers.
#
# Notation:
#   b, delta, r, eps, sigma, mu, u, l
#     Model parameters passed to make_params().
#   inv_mu
#     Synthetic sweep variable equal to 1 / mu, useful because 1 / mu is the
#     mean jump size under the exponential jump-size convention.
#   xL, xU
#     Lower and upper reflecting barriers.
#   xk, xl
#     Ambiguity thresholds for the drift and jump intensity.
#   gamma
#     Ergodic value returned by solve_optimal_barriers().
#
# Loading this file defines functions only.

.sweep_known_params <- c("b", "delta", "r", "eps", "sigma", "mu", "inv_mu", "u", "l")

.sweep_default_name <- function(sweep_param, sweep_values, prefix = "thresholds_vs") {
  step_str <- if (length(sweep_values) > 1L) {
    diffs <- unique(round(diff(sweep_values), 10))
    if (length(diffs) == 1L) sprintf("by_%g", diffs) else "custom_steps"
  } else {
    "single_value"
  }
  sprintf("%s_%s_%g_to_%g_%s",
          prefix, sweep_param, min(sweep_values), max(sweep_values), step_str)
}

.sweep_fixed_tag <- function(fixed_params, sweep_param) {
  fixed_for_name <- fixed_params
  if (sweep_param %in% names(fixed_for_name)) {
    fixed_for_name[[sweep_param]] <- NULL
  }
  if (identical(sweep_param, "inv_mu")) {
    fixed_for_name$mu <- NULL
  }
  if (!length(fixed_for_name)) {
    return("")
  }
  tag <- paste(
    sprintf(
      "%s_%s",
      names(fixed_for_name),
      vapply(fixed_for_name, function(x) sprintf("%.6g", x), character(1))
    ),
    collapse = "__"
  )
  gsub("\\s+", "", tag)
}

.sweep_name_with_fixed_tag <- function(out_name, fixed_params, sweep_param) {
  fixed_tag <- .sweep_fixed_tag(fixed_params, sweep_param)
  if (nzchar(fixed_tag)) paste0(out_name, "__", fixed_tag) else out_name
}

.sweep_axis_title <- function(sweep_param) {
  switch(
    sweep_param,
    b = expression(b),
    delta = expression(delta),
    r = expression(r),
    eps = expression(epsilon),
    sigma = expression(sigma),
    mu = expression(mu),
    inv_mu = expression(1 / mu),
    u = expression(u),
    l = expression(l),
    sweep_param
  )
}

.sweep_checked_margins <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 4L || any(!is.finite(x))) {
    stop(sprintf("%s must be a numeric vector of length 4.", name))
  }
  x
}

.sweep_apply_plot_par <- function(margins, axis_mgp, tcl,
                                  tick_cex, axis_title_cex,
                                  title_cex, base_cex, mex,
                                  outer_margins = NULL) {
  args <- list(
    xaxs = "i",
    mar = margins,
    mgp = axis_mgp,
    cex.axis = tick_cex,
    cex.lab = axis_title_cex,
    cex.main = title_cex
  )
  if (!is.null(tcl)) args$tcl <- tcl
  if (!is.null(base_cex)) args$cex <- base_cex
  if (!is.null(mex)) args$mex <- mex
  if (!is.null(outer_margins)) args$oma <- outer_margins
  do.call(par, args)
  invisible()
}

.sweep_render_and_save <- function(out_name, plotfun, show, save, out_dir,
                                   width_in, height_in, dpi, match_current) {
  if (exists("render_and_save", mode = "function")) {
    render_and_save(
      fname_base = out_name,
      plotfun = plotfun,
      show = show,
      save = save,
      dir = out_dir,
      width_in = width_in,
      height_in = height_in,
      dpi = dpi,
      match_current = match_current
    )
    return(invisible())
  }
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (isTRUE(show)) plotfun()
  if (isTRUE(save)) {
    grDevices::png(
      file.path(out_dir, paste0(out_name, ".png")),
      width = width_in,
      height = height_in,
      units = "in",
      res = dpi
    )
    tryCatch(plotfun(), finally = grDevices::dev.off())
    
    grDevices::pdf(
      file.path(out_dir, paste0(out_name, ".pdf")),
      width = width_in,
      height = height_in
    )
    tryCatch(plotfun(), finally = grDevices::dev.off())
  }
  invisible()
}

comparative_sweeper <- function(
    sweep_param = "b",
    sweep_values = seq(-10, 10, by = 0.1),
    b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1.0, l = 1.0,
    xL0 = -0.5, xk0 = -0.1, xl0 = 0.3, xU0 = 1.0,
    out_dir  = file.path("figures", "sweeps"),
    out_name = NULL,
    save = FALSE,
    solver_args = list(),
    verbose = TRUE
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  sweep_param <- as.character(sweep_param)
  sweep_values <- as.numeric(sweep_values)
  if (!length(sweep_values) || any(!is.finite(sweep_values))) {
    stop("sweep_values must be a non-empty finite numeric vector.")
  }
  if (!is.list(solver_args)) {
    stop("solver_args must be a list of extra arguments for solve_optimal_barriers().")
  }
  if (!sweep_param %in% .sweep_known_params) {
    stop(sprintf("sweep_param must be one of: %s", paste(.sweep_known_params, collapse = ", ")))
  }
  fixed_params <- list(b = b, delta = delta, r = r, eps = eps, sigma = sigma, mu = mu, u = u, l = l)
  
  res <- data.frame(
    sweep_param = sweep_param,
    sweep_value = sweep_values,
    xL = NA_real_, xk = NA_real_, xl = NA_real_, xU = NA_real_,
    gamma = NA_real_,
    converged = FALSE, error_msg = NA_character_,
    stringsAsFactors = FALSE
  )
  last_sol <- NULL
  
  for (i in seq_along(sweep_values)) {
    val <- sweep_values[i]
    if (isTRUE(verbose)) {
      message(sprintf("[sweep %s] %2d/%d  %s = % .6g ...",
                      sweep_param, i, length(sweep_values), sweep_param, val))
      flush.console()
    }
    
    p <- fixed_params
    if (sweep_param == "inv_mu") {
      if (val <= 0) {
        res$converged[i] <- FALSE
        res$error_msg[i] <- "inv_mu (1/mu) must be > 0"
        next
      }
      p$mu <- 1 / val
    } else {
      p[[sweep_param]] <- val
    }
    
    params <- do.call(make_params, p)
    
    # Adjacent sweep points usually have nearby barriers, so warm-start from the
    # previous successful solve to reduce nonlinear iterations.
    xL_init <- xL0; xk_init <- xk0; xl_init <- xl0; xU_init <- xU0
    if (!is.null(last_sol)) {
      xL_init <- last_sol$x$xL; xk_init <- last_sol$x$xk
      xl_init <- last_sol$x$xl; xU_init <- last_sol$x$xU
    }
    
    sol_opt <- tryCatch(
      do.call(
        solve_optimal_barriers,
        utils::modifyList(
          list(
            p = params,
            xL0 = xL_init,
            xk0 = xk_init,
            xl0 = xl_init,
            xU0 = xU_init,
            verbose = FALSE
          ),
          solver_args
        )
      ),
      error = function(e) e
    )
    
    if (inherits(sol_opt, "error")) {
      res$converged[i] <- FALSE
      res$error_msg[i] <- conditionMessage(sol_opt)
    } else {
      th <- sol_opt$x
      res[i, c("xL", "xk", "xl", "xU")] <- unlist(th[c("xL", "xk", "xl", "xU")])
      res$gamma[i] <- sol_opt$gamma
      res$converged[i] <- TRUE
      res$error_msg[i] <- NA_character_
      last_sol <- sol_opt
    }
  }
  
  if (is.null(out_name)) {
    out_name <- .sweep_default_name(sweep_param, sweep_values)
  }
  
  if (save) {
    csv_base <- .sweep_name_with_fixed_tag(out_name, fixed_params, sweep_param)
    csv_path <- file.path(out_dir, paste0(csv_base, ".csv"))
    write.csv(res, csv_path, row.names = FALSE)
    message(sprintf("Saved sweep metrics to: %s", normalizePath(csv_path)))
  }
  
  list(
    results = res,
    out_dir = out_dir,
    out_name = out_name,
    fixed_params = fixed_params
  )
}

plot_sweep <- function(
    sweep_obj,
    out_dir  = file.path("figures", "sweeps"),
    out_name = NULL,
    show = interactive(),
    save = FALSE,
    width_in = 6.67,
    height_in = 4.67,
    dpi = 300,
    match_current = FALSE,
    cols = c("#B22222", "#1E90FF", "#8A2BE2", "#228B22"),
    ltys = c(1, 1, 1, 1),
    lwds = 2,
    title = TRUE,
    show_x_axis_title = TRUE,
    show_y_axis_title = TRUE,
    x_axis_title = NULL,
    y_axis_title = NULL,
    show_tick_labels = TRUE,
    tick_cex = 1.4,
    axis_title_cex = 1,
    title_cex = 1,
    base_cex = NULL,
    mex = NULL,
    margins = c(2.2, 2.2, 1.5, 1.2),
    outer_margins = NULL,
    axis_mgp = c(3, 1, 0),
    tcl = NULL,
    legend_position = "topleft",
    legend_cex = NULL,
    legend_bty = "n",
    legend_horiz = TRUE,
    legend_seg_len = 2,
    legend_x_intersp = 0.5,
    plot_gamma = FALSE,
    gamma_layout = c("stacked", "separate")
) {
  gamma_layout <- match.arg(gamma_layout)
  margins <- .sweep_checked_margins(margins, "margins")
  if (!is.null(outer_margins)) {
    outer_margins <- .sweep_checked_margins(outer_margins, "outer_margins")
  }
  bottom <- margins[1]; left <- margins[2]; top <- margins[3]; right <- margins[4]
  
  res <- if (is.data.frame(sweep_obj)) {
    sweep_obj
  } else {
    sweep_obj$results %||% sweep_obj
  }
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  needed <- c("sweep_param", "sweep_value", "xL", "xk", "xl", "xU")
  if (!all(needed %in% names(res))) {
    stop("sweep_obj must contain sweep_param, sweep_value, xL, xk, xl, and xU.")
  }
  ok <- if ("converged" %in% names(res)) {
    as.logical(res$converged)
  } else {
    stats::complete.cases(res[, c("xL", "xk", "xl", "xU")])
  }
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) {
    warning("No successful solutions; nothing to plot.")
    return(invisible(NULL))
  }
  
  sweep_param <- as.character(unique(res$sweep_param))
  if (length(sweep_param) != 1L) {
    warning("Multiple sweep_param values detected; using the first.")
    sweep_param <- sweep_param[1L]
  }
  
  sweep_values <- as.numeric(res$sweep_value)
  if (is.null(out_name)) {
    out_name <- .sweep_default_name(sweep_param, sweep_values)
  }
  
  fixed_for_name <- NULL
  if (!is.data.frame(sweep_obj) && !is.null(sweep_obj$fixed_params)) {
    fixed_for_name <- sweep_obj$fixed_params
  }
  if (!is.null(fixed_for_name)) {
    out_name <- .sweep_name_with_fixed_tag(out_name, fixed_for_name, sweep_param)
  }
  
  param_label_str <- if (sweep_param == "inv_mu") "1/mu" else sweep_param
  if (is.null(x_axis_title)) {
    x_axis_title <- .sweep_axis_title(sweep_param)
  }
  
  default_y_axis_title_thresholds <- "threshold value"
  default_y_axis_title_gamma <- expression(gamma)
  if (is.null(y_axis_title)) {
    y_axis_title_thresholds <- default_y_axis_title_thresholds
    y_axis_title_gamma <- default_y_axis_title_gamma
  } else if (is.list(y_axis_title)) {
    if (is.null(names(y_axis_title))) {
      y_axis_title_thresholds <- y_axis_title[[1]]
      y_axis_title_gamma <- if (length(y_axis_title) >= 2L) y_axis_title[[2]] else y_axis_title[[1]]
    } else {
      y_axis_title_thresholds <- y_axis_title[["thresholds"]] %||% default_y_axis_title_thresholds
      y_axis_title_gamma <- y_axis_title[["gamma"]] %||% default_y_axis_title_gamma
    }
  } else {
    y_axis_title_thresholds <- y_axis_title
    y_axis_title_gamma <- y_axis_title
  }
  
  ylm_th <- range(as.numeric(unlist(res[ok, c("xL", "xk", "xl", "xU")])), na.rm = TRUE)
  has_gamma <- "gamma" %in% names(res) && any(is.finite(res$gamma[ok]))
  if (!has_gamma) plot_gamma <- FALSE
  if (plot_gamma) {
    ylm_g <- range(res$gamma[ok], na.rm = TRUE)
  }
  
  cols <- rep_len(cols, 4L)
  ltys <- rep_len(ltys, 4L)
  lwds <- rep_len(lwds, 4L)
  threshold_names <- c("xL", "xk", "xl", "xU")
  threshold_labels <- expression(underline(x), x^kappa, x^lambda, bar(x))
  
  title_str <- if (isTRUE(title)) {
    sprintf("Ambiguity thresholds & barriers vs %s", param_label_str)
  } else {
    ""
  }
  title_gamma <- if (isTRUE(title)) {
    sprintf("Ergodic value %s vs %s", "\u03b3", param_label_str)
  } else {
    ""
  }
  
  thresholds_panel <- function(show_x_title = TRUE, mar_override = NULL) {
    op <- par(c("mar", "mgp", "tcl", "cex.axis", "cex.lab", "cex.main",
                "cex", "mex", "xaxs"))
    on.exit(par(op), add = TRUE)
    mar <- if (is.null(mar_override)) margins else mar_override
    .sweep_apply_plot_par(
      margins = mar,
      axis_mgp = axis_mgp,
      tcl = tcl,
      tick_cex = tick_cex,
      axis_title_cex = axis_title_cex,
      title_cex = title_cex,
      base_cex = base_cex,
      mex = mex,
      outer_margins = outer_margins
    )
    
    plot(
      sweep_values, res$xL,
      type = "n",
      xlab = if (show_x_title && show_x_axis_title) x_axis_title else "",
      ylab = if (show_y_axis_title) y_axis_title_thresholds else "",
      ylim = ylm_th,
      xaxt = if (isTRUE(show_tick_labels)) "s" else "n",
      yaxt = if (isTRUE(show_tick_labels)) "s" else "n"
    )
    grid()
    for (i in seq_along(threshold_names)) {
      lines(sweep_values, res[[threshold_names[i]]],
            lwd = lwds[i], lty = ltys[i], col = cols[i])
    }
    
    legend_args <- list(
      legend = threshold_labels,
      col = cols,
      lty = ltys,
      lwd = lwds,
      bty = legend_bty,
      horiz = legend_horiz,
      x.intersp = legend_x_intersp,
      seg.len = legend_seg_len
    )
    if (is.numeric(legend_position) && length(legend_position) == 2L &&
        all(is.finite(legend_position))) {
      legend_args$x <- legend_position[1]
      legend_args$y <- legend_position[2]
    } else {
      legend_args$x <- legend_position
    }
    if (!is.null(legend_cex)) legend_args$cex <- legend_cex
    do.call(legend, legend_args)
    title(title_str)
  }
  
  gamma_panel <- function(show_x_title = FALSE, show_x_axis_ticks = TRUE,
                          mar_override = NULL) {
    op <- par(c("mar", "mgp", "tcl", "cex.axis", "cex.lab", "cex.main",
                "cex", "mex", "xaxs"))
    on.exit(par(op), add = TRUE)
    mar <- if (is.null(mar_override)) margins else mar_override
    .sweep_apply_plot_par(
      margins = mar,
      axis_mgp = axis_mgp,
      tcl = tcl,
      tick_cex = tick_cex,
      axis_title_cex = axis_title_cex,
      title_cex = title_cex,
      base_cex = base_cex,
      mex = mex,
      outer_margins = outer_margins
    )
    
    plot(
      sweep_values, res$gamma,
      type = "l",
      lwd = 2,
      xlab = if (show_x_title && show_x_axis_title) x_axis_title else "",
      ylab = if (show_y_axis_title) y_axis_title_gamma else "",
      ylim = ylm_g,
      xaxt = if (show_x_axis_ticks && show_tick_labels) "s" else "n",
      yaxt = if (show_tick_labels) "s" else "n"
    )
    grid()
    legend("topleft", legend = expression(gamma), lwd = 2, bty = "n")
    title(title_gamma)
  }
  
  if (!plot_gamma) {
    plotfun <- function() thresholds_panel(show_x_title = TRUE)
    .sweep_render_and_save(out_name, plotfun, show, save, out_dir,
                           width_in, height_in, dpi, match_current)
  } else if (gamma_layout == "separate") {
    .sweep_render_and_save(
      out_name,
      function() thresholds_panel(show_x_title = TRUE),
      show, save, out_dir, width_in, height_in, dpi, match_current
    )
    .sweep_render_and_save(
      paste0(out_name, "_gamma"),
      function() gamma_panel(show_x_title = FALSE, show_x_axis_ticks = TRUE),
      show, save, out_dir, width_in, height_in, dpi, match_current
    )
  } else {
    plotfun <- function() {
      layout(matrix(1:2, nrow = 2), heights = c(0.35, 0.65))
      on.exit(layout(1), add = TRUE)
      
      gamma_panel(
        show_x_title = FALSE,
        show_x_axis_ticks = FALSE,
        mar_override = c(0.4, left, top, right)
      )
      thresholds_panel(
        show_x_title = TRUE,
        mar_override = c(bottom, left, 0.4, right)
      )
    }
    .sweep_render_and_save(out_name, plotfun, show, save, out_dir,
                           width_in, height_in, dpi, match_current)
  }
  
  invisible(NULL)
}
