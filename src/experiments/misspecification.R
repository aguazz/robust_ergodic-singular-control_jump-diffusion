# Misspecification-grid experiment helpers.
#
# Notation:
#   delta, eps
#     Ambiguity radii used to define the model solved at each grid point.
#   gamma_with_ambiguity
#     Ergodic value of the robust policy optimized for the ambiguous model.
#   gamma_zero_ambiguity_policy
#     Ergodic value obtained by fixing the certainty-model reflecting barriers,
#     re-solving only the ambiguity thresholds, and evaluating that policy under
#     the ambiguous model.
#   relative_cost_pct
#     Public reporting surface for the robust misspecification cost (RMC),
#     computed as 100 * (gamma_zero_ambiguity_policy - gamma_with_ambiguity) /
#     gamma_with_ambiguity.
#
# Loading this file defines functions only.
.safe_relative_change <- function(num, den, tol = 1e-12) {
  out <- rep(NA_real_, length(num))
  ok <- is.finite(num) & is.finite(den) & (abs(den) > tol)
  out[ok] <- (num[ok] - den[ok]) / den[ok]
  out
}

# Public-facing reporting now defines RMC in percentage terms. Keep the
# unscaled relative-cost ratio internally, but silently map public reporting
# requests onto the percentage-scaled surface.
.misspecification_reporting_surface <- function(surface) {
  if (identical(surface, "relative_cost")) "relative_cost_pct" else surface
}

.ensure_misspecification_surface_has_finite_values <- function(grid_obj, surface) {
  zmat <- grid_obj[[surface]]
  if (is.null(zmat) || !is.matrix(zmat)) {
    stop(sprintf("Surface '%s' is not available in grid_obj.", surface))
  }
  if (any(is.finite(zmat))) {
    return(invisible(zmat))
  }
  
  msg <- sprintf("Surface '%s' does not contain any finite values.", surface)
  certainty_error_msg <- grid_obj$certainty_solution_error
  certainty_failed <- is.null(grid_obj$certainty_solution) &&
    is.character(certainty_error_msg) &&
    length(certainty_error_msg) >= 1L &&
    !is.na(certainty_error_msg[1]) &&
    nzchar(certainty_error_msg[1])
  
  if (certainty_failed && surface %in% c(
    "relative_cost_pct",
    "relative_cost",
    "gamma_difference",
    "gamma_zero_ambiguity_policy"
  )) {
    msg <- paste0(
      msg,
      " The certainty baseline could not be solved, so the RMC comparison surface is undefined. ",
      "certainty_solution_error: ", certainty_error_msg[1]
    )
    
    cp <- grid_obj$certainty_params
    if (!is.null(cp$b) && isTRUE(all.equal(as.numeric(cp$b), 0))) {
      msg <- paste0(
        msg,
        " At certainty (delta = 0, epsilon = 0), this means a_{i,1} = mu * b = 0, ",
        "which is exactly the unimplemented S1 cubic branch."
      )
    }
  }
  
  stop(msg)
}

.misspecification_rdata_path <- function(out_dir, out_name) {
  file.path(out_dir, paste0(out_name, ".RData"))
}

save_misspecification_cost_rdata <- function(
    grid_obj,
    out_dir = grid_obj$out_dir %||% file.path("figures", "misspecification"),
    out_name = NULL,
    object_name = "grid_obj"
) {
  if (is.null(out_name)) {
    out_name <- grid_obj$out_name %||% "misspecification"
  }
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_path <- .misspecification_rdata_path(out_dir, out_name)
  grid_obj$out_dir <- out_dir
  grid_obj$out_name <- out_name
  grid_obj$rdata_path <- out_path
  
  save_env <- new.env(parent = emptyenv())
  assign(object_name, grid_obj, envir = save_env)
  save(list = object_name, file = out_path, envir = save_env)
  
  message(sprintf("Saved misspecification grid to: %s", normalizePath(out_path)))
  invisible(out_path)
}

load_misspecification_cost_rdata <- function(
    path = NULL,
    out_dir = file.path("figures", "misspecification"),
    out_name = NULL,
    object_name = NULL
) {
  if (is.null(path)) {
    if (is.null(out_name)) {
      stop("Provide either 'path' or 'out_name' to load a misspecification grid.")
    }
    path <- .misspecification_rdata_path(out_dir, out_name)
  }
  
  if (!file.exists(path)) {
    stop(sprintf("Misspecification grid file not found: %s", path))
  }
  
  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = load_env)
  
  if (is.null(object_name)) {
    if (length(loaded_names) != 1L) {
      stop("The .RData file contains multiple objects; please specify 'object_name'.")
    }
    object_name <- loaded_names[[1]]
  }
  
  if (!exists(object_name, envir = load_env, inherits = FALSE)) {
    stop(sprintf("Object '%s' was not found in %s", object_name, path))
  }
  
  grid_obj <- get(object_name, envir = load_env, inherits = FALSE)
  if (!is.list(grid_obj)) {
    stop(sprintf("Object '%s' in %s is not a misspecification grid object.", object_name, path))
  }
  
  if (is.null(grid_obj$out_dir) || !nzchar(grid_obj$out_dir)) {
    grid_obj$out_dir <- dirname(path)
  }
  if (is.null(grid_obj$out_name) || !nzchar(grid_obj$out_name)) {
    grid_obj$out_name <- tools::file_path_sans_ext(basename(path))
  }
  grid_obj$rdata_path <- path
  class(grid_obj) <- unique(c("misspecification_cost_grid", class(grid_obj)))
  
  grid_obj
}

misspecification_cost_grid <- function(
    delta_values,
    eps_values,
    b = 0, r = 1, sigma = 1, mu = 1, u = 1.0, l = 1.0,
    certainty_delta = 0, certainty_eps = 0,
    xL0 = -0.5, xk0 = -0.1, xl0 = 0.3, xU0 = 1.0,
    method = "Broyden",
    control = list(ftol=1e-11, xtol=1e-11, maxit=1e6, stepmax=1, trace=0),
    tol_build = 1e-10,
    solver_verbose = FALSE,
    tol_regime_switch = 1e-3,
    switch_to_newton = FALSE,
    switch_ftol = 5e-2,
    switch_iter_min = 2L,
    newton_control = list(ftol=1e-12, xtol=1e-12, maxit=1e5, stepmax=1, trace=0),
    fd_step = 1e-6,
    relative_tol = 1e-12,
    verbose = TRUE,
    out_dir = file.path("figures", "misspecification"),
    out_name = NULL,
    save = FALSE
) {
  delta_values <- as.numeric(delta_values)
  eps_values <- as.numeric(eps_values)
  
  if (!length(delta_values) || !length(eps_values)) {
    stop("delta_values and eps_values must both be non-empty.")
  }
  if (any(!is.finite(delta_values)) || any(delta_values < 0)) {
    stop("delta_values must be finite and non-negative.")
  }
  if (any(!is.finite(eps_values)) || any(eps_values < 0 | eps_values > 1)) {
    stop("eps_values must be finite and lie in [0, 1].")
  }
  
  certainty_params <- make_params(
    b = b, delta = certainty_delta, r = r, eps = certainty_eps,
    sigma = sigma, mu = mu, u = u, l = l
  )
  
  certainty_sol_try <- tryCatch(
    solve_optimal_barriers(
      certainty_params, xL0, xk0, xl0, xU0,
      method = method,
      control = control,
      tol_build = tol_build,
      verbose = solver_verbose,
      tol_regime_switch = tol_regime_switch,
      switch_to_newton = switch_to_newton,
      switch_ftol = switch_ftol,
      switch_iter_min = switch_iter_min,
      newton_control = newton_control,
      fd_step = fd_step
    ),
    error = function(e) e
  )
  certainty_sol <- if (inherits(certainty_sol_try, "error")) NULL else certainty_sol_try
  certainty_error_msg <- if (inherits(certainty_sol_try, "error")) {
    conditionMessage(certainty_sol_try)
  } else {
    NA_character_
  }
  
  dim_labs <- function(x) format(signif(x, 8), trim = TRUE, scientific = FALSE)
  dims <- list(delta = dim_labs(delta_values), eps = dim_labs(eps_values))
  
  gamma_with_ambiguity <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                                 dimnames = dims)
  gamma_zero_ambiguity_policy <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                                        dimnames = dims)
  gamma_difference <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                             dimnames = dims)
  relative_cost <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                          dimnames = dims)
  relative_cost_pct <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                              dimnames = dims)
  robust_converged <- matrix(FALSE, nrow = length(delta_values), ncol = length(eps_values),
                             dimnames = dims)
  zero_policy_converged <- matrix(FALSE, nrow = length(delta_values), ncol = length(eps_values),
                                  dimnames = dims)
  error_msg <- matrix("", nrow = length(delta_values), ncol = length(eps_values),
                      dimnames = dims)
  
  append_error <- function(existing, msg) {
    if (!nzchar(existing)) msg else paste(existing, msg, sep = " | ")
  }
  start_vals <- function(sol) {
    if (!is.null(sol) && !is.null(sol$x)) {
      return(sol$x[c("xL", "xk", "xl", "xU")])
    }
    c(xL = xL0, xk = xk0, xl = xl0, xU = xU0)
  }
  
  total <- length(delta_values) * length(eps_values)
  counter <- 0L
  
  for (i in seq_along(delta_values)) {
    last_robust_sol <- certainty_sol
    last_zero_sol <- certainty_sol
    
    for (j in seq_along(eps_values)) {
      counter <- counter + 1L
      delta_i <- delta_values[i]
      eps_j <- eps_values[j]
      
      if (isTRUE(verbose)) {
        message(sprintf("[misspecification] %3d/%d  delta = %.6g, eps = %.6g",
                        counter, total, delta_i, eps_j))
        flush.console()
      }
      
      amb_params <- make_params(
        b = b, delta = delta_i, r = r, eps = eps_j,
        sigma = sigma, mu = mu, u = u, l = l
      )
      
      same_as_certainty <-
        isTRUE(all.equal(delta_i, certainty_delta)) &&
        isTRUE(all.equal(eps_j, certainty_eps))
      robust_start <- start_vals(last_robust_sol)
      
      robust_try <- tryCatch(
        if (same_as_certainty && !is.null(certainty_sol)) {
          certainty_sol
        } else {
          solve_optimal_barriers(
            amb_params,
            robust_start[["xL"]],
            robust_start[["xk"]],
            robust_start[["xl"]],
            robust_start[["xU"]],
            method = method,
            control = control,
            tol_build = tol_build,
            verbose = solver_verbose,
            tol_regime_switch = tol_regime_switch,
            switch_to_newton = switch_to_newton,
            switch_ftol = switch_ftol,
            switch_iter_min = switch_iter_min,
            newton_control = newton_control,
            fd_step = fd_step
          )
        },
        error = function(e) e
      )
      
      if (inherits(robust_try, "error")) {
        error_msg[i, j] <- append_error(
          error_msg[i, j],
          paste0("robust policy: ", conditionMessage(robust_try))
        )
      } else {
        gamma_with_ambiguity[i, j] <- robust_try$gamma
        robust_converged[i, j] <- TRUE
        last_robust_sol <- robust_try
      }
      
      if (is.null(certainty_sol)) {
        zero_try <- structure(
          list(message = paste0(
            "certainty baseline unavailable: ", certainty_error_msg
          )),
          class = c("simpleError", "error", "condition")
        )
      } else {
        zero_start <- start_vals(last_zero_sol)
        zero_try <- tryCatch(
          if (same_as_certainty) {
            certainty_sol
          } else {
            solve_ambiguity_thresholds_fixed_barriers(
              amb_params,
              certainty_sol$x$xL,
              certainty_sol$x$xU,
              zero_start[["xk"]],
              zero_start[["xl"]],
              method = method,
              control = control,
              tol_build = tol_build,
              verbose = solver_verbose,
              tol_regime_switch = tol_regime_switch,
              switch_to_newton = switch_to_newton,
              switch_ftol = switch_ftol,
              switch_iter_min = switch_iter_min,
              newton_control = newton_control,
              fd_step = fd_step
            )
          },
          error = function(e) e
        )
      }
      
      if (inherits(zero_try, "error")) {
        error_msg[i, j] <- append_error(
          error_msg[i, j],
          paste0("certainty-barriers policy: ", conditionMessage(zero_try))
        )
      } else {
        gamma_zero_ambiguity_policy[i, j] <- zero_try$gamma
        zero_policy_converged[i, j] <- TRUE
        last_zero_sol <- zero_try
      }
      
      if (robust_converged[i, j] && zero_policy_converged[i, j]) {
        gamma_difference[i, j] <-
          gamma_zero_ambiguity_policy[i, j] - gamma_with_ambiguity[i, j]
        relative_cost[i, j] <-
          .safe_relative_change(gamma_zero_ambiguity_policy[i, j],
                                gamma_with_ambiguity[i, j],
                                tol = relative_tol)
        relative_cost_pct[i, j] <- 100 * relative_cost[i, j]
      }
    }
  }
  
  long_results <- data.frame(
    delta = rep(delta_values, times = length(eps_values)),
    eps = rep(eps_values, each = length(delta_values)),
    gamma_robust = as.vector(gamma_with_ambiguity),
    gamma_nonrobust_worst_case = as.vector(gamma_zero_ambiguity_policy),
    gamma_with_ambiguity = as.vector(gamma_with_ambiguity),
    gamma_zero_ambiguity_policy = as.vector(gamma_zero_ambiguity_policy),
    gamma_difference = as.vector(gamma_difference),
    relative_cost = as.vector(relative_cost),
    relative_cost_pct = as.vector(relative_cost_pct),
    robust_converged = as.vector(robust_converged),
    zero_policy_converged = as.vector(zero_policy_converged),
    converged = as.vector(robust_converged & zero_policy_converged),
    error_msg = as.vector(error_msg),
    stringsAsFactors = FALSE
  )
  long_results$error_msg[!nzchar(long_results$error_msg)] <- NA_character_
  
  if (is.null(out_name)) {
    out_name <- sprintf("misspecification_delta_%g_to_%g__eps_%g_to_%g",
                        min(delta_values), max(delta_values),
                        min(eps_values), max(eps_values))
  }
  
  out <- list(
    delta_values = delta_values,
    eps_values = eps_values,
    certainty_params = certainty_params,
    certainty_solution = certainty_sol,
    certainty_solution_error = certainty_error_msg,
    gamma_robust = gamma_with_ambiguity,
    gamma_nonrobust_worst_case = gamma_zero_ambiguity_policy,
    gamma_with_ambiguity = gamma_with_ambiguity,
    gamma_zero_ambiguity_policy = gamma_zero_ambiguity_policy,
    gamma_difference = gamma_difference,
    relative_cost = relative_cost,
    relative_cost_pct = relative_cost_pct,
    robust_converged = robust_converged,
    zero_policy_converged = zero_policy_converged,
    error_msg = error_msg,
    long_results = long_results,
    out_dir = out_dir,
    out_name = out_name,
    rdata_path = NULL
  )
  class(out) <- c("misspecification_cost_grid", class(out))
  
  if (save) {
    out$rdata_path <- save_misspecification_cost_rdata(
      out,
      out_dir = out_dir,
      out_name = out_name
    )
  }
  
  out
}

export_misspecification_cost_latex <- function(
    grid_obj,
    surface = c("relative_cost_pct", "relative_cost"),
    out_dir = file.path("figures", "misspecification"),
    out_name = NULL,
    save = TRUE,
    digits = NULL,
    format = "f",
    na_string = "--",
    include_table_env = TRUE,
    placement = "htbp",
    centering = TRUE,
    booktabs = TRUE,
    caption = NULL,
    label = NULL,
    row_header = "\\diagbox[width=1.6cm,height=0.7cm]{$\\delta$}{$\\epsilon$}",
    column_align = NULL
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  
  if (is.null(grid_obj[[surface]]) || !is.matrix(grid_obj[[surface]])) {
    stop(sprintf("Surface '%s' is not available in grid_obj.", surface))
  }
  
  zmat <- grid_obj[[surface]]
  if (is.null(digits)) {
    digits <- 2L
  }
  
  if (is.null(caption)) {
    caption <- "RMC matrix"
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_table")
  }
  
  if (is.null(label) && isTRUE(include_table_env)) {
    label <- paste0("tab:", gsub("[^A-Za-z0-9]+", "_", out_name))
  }
  
  row_labels <- rownames(zmat)
  if (is.null(row_labels)) {
    row_labels <- format(signif(seq_len(nrow(zmat)), 8), trim = TRUE, scientific = FALSE)
  }
  col_labels <- colnames(zmat)
  if (is.null(col_labels)) {
    col_labels <- format(signif(seq_len(ncol(zmat)), 8), trim = TRUE, scientific = FALSE)
  }
  
  if (is.null(column_align)) {
    column_align <- paste0("r|", paste(rep("r", ncol(zmat)), collapse = ""))
  }
  
  format_cell <- function(x) {
    out <- rep(na_string, length(x))
    ok <- is.finite(x)
    if (any(ok)) {
      xx <- x[ok]
      xx[abs(xx) < 0.5 * 10^(-digits)] <- 0
      out[ok] <- formatC(xx, format = format, digits = digits)
    }
    out
  }
  
  header_row <- paste0(
    paste(c(row_header, col_labels), collapse = " & "),
    " \\\\"
  )
  body_rows <- vapply(
    seq_len(nrow(zmat)),
    function(i) {
      paste0(
        paste(c(row_labels[i], format_cell(zmat[i, ])), collapse = " & "),
        " \\\\"
      )
    },
    character(1)
  )
  
  lines <- character(0)
  if (isTRUE(include_table_env)) {
    lines <- c(lines, paste0("\\begin{table}[", placement, "]"))
    if (isTRUE(centering)) lines <- c(lines, "\\centering")
    if (!is.null(caption)) lines <- c(lines, paste0("\\caption{", caption, "}"))
    if (!is.null(label))   lines <- c(lines, paste0("\\label{", label, "}"))
  }
  
  lines <- c(lines, paste0("\\begin{tabular}{", column_align, "}"))
  lines <- c(lines, if (isTRUE(booktabs)) "\\toprule" else "\\hline")
  lines <- c(lines, header_row)
  lines <- c(lines, if (isTRUE(booktabs)) "\\midrule" else "\\hline")
  lines <- c(lines, body_rows)
  lines <- c(lines, if (isTRUE(booktabs)) "\\bottomrule" else "\\hline")
  lines <- c(lines, "\\end{tabular}")
  
  if (isTRUE(include_table_env)) {
    lines <- c(lines, "\\end{table}")
  }
  
  latex <- paste(lines, collapse = "\n")
  
  out_path <- NULL
  if (isTRUE(save)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_path <- file.path(out_dir, paste0(out_name, ".tex"))
    writeLines(lines, out_path)
    message(sprintf("Saved LaTeX table to: %s", normalizePath(out_path)))
  }
  
  invisible(list(
    latex = latex,
    path = out_path,
    matrix = zmat,
    surface = surface
  ))
}

.format_misspecification_param_value <- function(x) {
  gsub("\\s+", "", format(signif(as.numeric(x), 8), trim = TRUE, scientific = FALSE))
}

.format_misspecification_tag_value <- function(x) {
  x <- .format_misspecification_param_value(x)
  x <- gsub("-", "minus_", x, fixed = TRUE)
  x <- gsub("\\.", "p", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

.misspecification_example_tag <- function(b_value, r) {
  paste0(
    "b_", .format_misspecification_tag_value(b_value),
    "__r_", .format_misspecification_tag_value(r)
  )
}

.resolve_misspecification_example_tag <- function(tag = NULL, b_value = NULL, r = NULL) {
  if (!is.null(tag)) {
    tag <- as.character(tag)[1]
    if (nzchar(tag)) {
      return(tag)
    }
  }
  if (is.null(b_value) || is.null(r)) {
    stop("Provide either 'tag' or both 'b_value' and 'r' to identify the misspecification example.")
  }
  .misspecification_example_tag(b_value = b_value, r = r)
}

.misspecification_param_suffix <- function(sigma, mu, u, l) {
  paste0(
    "__sigma_", .format_misspecification_param_value(sigma),
    "__mu_", .format_misspecification_param_value(mu),
    "__u_", .format_misspecification_param_value(u),
    "__l_", .format_misspecification_param_value(l)
  )
}

.misspecification_cost_example_spec <- function(
    tag = NULL, b_value = NULL,
    r = 1, sigma = 1, mu = 1, u = 1, l = 1
) {
  tag <- .resolve_misspecification_example_tag(tag = tag, b_value = b_value, r = r)
  base_name <- paste0("misspec_", tag, .misspecification_param_suffix(sigma, mu, u, l))
  base_out_dir <- file.path("figures", "misspecification", tag)
  
  list(
    out_dir = base_out_dir,
    base_name = base_name,
    main_grid = list(
      delta_values = seq(0, 1, by = 0.01),
      eps_values = seq(0, 1, by = 0.01),
      out_name = paste0(base_name, "_surface")
    ),
    grid_vs_delta = list(
      delta_values = seq(0, 1, by = 0.025),
      eps_values = seq(0, 1, by = 0.2),
      out_name = paste0(base_name, "_vs_delta_data")
    ),
    grid_vs_epsilon = list(
      delta_values = seq(0, 1, by = 0.2),
      eps_values = seq(0, 1, by = 0.025),
      out_name = paste0(base_name, "_vs_epsilon_data")
    ),
    text_grid = list(
      delta_values = seq(0, 1, by = 0.2),
      eps_values = seq(0, 1, by = 0.2),
      out_name = paste0(base_name, "_table_data")
    )
  )
}

generate_misspecification_cost_example_data <- function(
    b_value, tag = NULL,
    r = 1, sigma = 1, mu = 1, u = 1, l = 1,
    force = FALSE
) {
  tag <- .resolve_misspecification_example_tag(tag = tag, b_value = b_value, r = r)
  spec <- .misspecification_cost_example_spec(
    tag = tag,
    b_value = b_value,
    r = r, sigma = sigma, mu = mu, u = u, l = l
  )
  grid_names <- c("main_grid", "grid_vs_delta", "grid_vs_epsilon", "text_grid")
  
  grid_common_args <- list(
    b = b_value, r = r, sigma = sigma, mu = mu, u = u, l = l,
    xL0 = -0.5, xk0 = -0.1, xl0 = 0.5, xU0 = 1,
    verbose = FALSE,
    out_dir = spec$out_dir,
    save = TRUE
  )
  
  out_paths <- setNames(vector("list", length(grid_names)), grid_names)
  
  for (grid_name in grid_names) {
    grid_args <- spec[[grid_name]]
    out_path <- .misspecification_rdata_path(spec$out_dir, grid_args$out_name)
    
    if (file.exists(out_path) && !isTRUE(force)) {
      message(sprintf("Using cached misspecification grid: %s", normalizePath(out_path)))
      out_paths[[grid_name]] <- out_path
      next
    }
    
    grid_obj <- do.call(
      misspecification_cost_grid,
      c(grid_args, grid_common_args)
    )
    out_paths[[grid_name]] <- grid_obj$rdata_path
  }
  
  invisible(out_paths)
}

load_misspecification_cost_example_data <- function(
    tag = NULL,
    b_value = NULL,
    r = 1, sigma = 1, mu = 1, u = 1, l = 1
) {
  tag <- .resolve_misspecification_example_tag(tag = tag, b_value = b_value, r = r)
  spec <- .misspecification_cost_example_spec(
    tag = tag,
    b_value = b_value,
    r = r, sigma = sigma, mu = mu, u = u, l = l
  )
  load_one <- function(grid_args) {
    load_misspecification_cost_rdata(
      out_dir = spec$out_dir,
      out_name = grid_args$out_name
    )
  }
  
  list(
    main_grid = load_one(spec$main_grid),
    grid_vs_delta = load_one(spec$grid_vs_delta),
    grid_vs_epsilon = load_one(spec$grid_vs_epsilon),
    text_grid = load_one(spec$text_grid)
  )
}

.misspecification_grid_model_params <- function(grid_obj) {
  cp <- grid_obj$certainty_params
  needed <- c("b", "r", "sigma", "mu", "u", "l")
  if (is.null(cp) || !all(needed %in% names(cp))) {
    return(NULL)
  }
  as.list(cp[needed])
}

# Build the expensive grids once; rerun with force = TRUE when either the grid
# setup or the model parameters change.
misspec_example_common_params <- list(
  sigma = 1,
  mu = 1,
  u = 1,
  l = 1
)

misspec_example_specs <- list(
  list(
    b_value = 2,
    r = 1,
    vs_delta_label_x = 0.9,
    canvas_panel_layout = "3x1",
    canvas_row_heights = c(1, 1, 1)
  ),
  list(
    b_value = -2,
    r = 1,
    vs_delta_label_x = 0.2,
    canvas_panel_layout = "3x1",
    canvas_row_heights = c(1, 1, 1)
  )
)

run_misspecification_cost_example <- function(
    example_spec,
    common_params = misspec_example_common_params,
    force = FALSE
) {
  stopifnot(is.list(example_spec), !is.null(example_spec$b_value))
  
  model_param_names <- c("r", "sigma", "mu", "u", "l")
  plot_param_names <- c(
    "surface", "table_surface", "vs_delta_label_x",
    "canvas_panel_layout", "canvas_column_widths", "canvas_row_heights",
    "canvas_width_in", "canvas_height_in"
  )
  
  example_model_params <- utils::modifyList(
    common_params,
    example_spec[intersect(names(example_spec), model_param_names)]
  )
  example_model_params$r <- example_model_params$r %||% 1
  example_plot_params <- example_spec[intersect(names(example_spec), plot_param_names)]
  tag <- .resolve_misspecification_example_tag(
    tag = example_spec$tag %||% NULL,
    b_value = example_spec$b_value,
    r = example_model_params$r
  )
  
  paths <- do.call(
    generate_misspecification_cost_example_data,
    c(
      list(
        b_value = example_spec$b_value,
        tag = tag,
        force = force
      ),
      example_model_params
    )
  )
  
  # Reload cached .RData grids when only plot styling changes.
  result <- do.call(
    plot_misspecification_cost_example_set,
    c(
      list(tag = tag, b_value = example_spec$b_value),
      example_model_params,
      example_plot_params
    )
  )
  
  list(
    spec = utils::modifyList(
      example_model_params,
      utils::modifyList(example_spec, list(tag = tag))
    ),
    paths = paths,
    result = result
  )
}
