# Table-formatting helpers for the companion Shiny app.

app_fmt_num <- function(x, digits = 7) {
  ifelse(
    is.finite(as.numeric(x)),
    formatC(as.numeric(x), format = "fg", digits = digits),
    ifelse(is.na(x), NA_character_, as.character(x))
  )
}

app_fmt_sci <- function(x, digits = 3) {
  ifelse(
    is.finite(as.numeric(x)),
    formatC(as.numeric(x), format = "e", digits = digits),
    ifelse(is.na(x), "n/a", as.character(x))
  )
}

app_fmt_key_num <- function(x, digits = 6, sci_digits = 3) {
  numeric_x <- as.numeric(x)
  out <- rep(NA_character_, length(numeric_x))
  finite <- is.finite(numeric_x)
  tiny <- finite & numeric_x != 0 & abs(numeric_x) < 10^(-digits)
  ordinary <- finite & !tiny

  out[tiny] <- formatC(numeric_x[tiny], format = "e", digits = sci_digits)
  out[ordinary] <- sub(
    "\\.?0+$",
    "",
    formatC(numeric_x[ordinary], format = "f", digits = digits)
  )
  out[!finite] <- ifelse(is.na(x[!finite]), NA_character_, as.character(x[!finite]))
  out
}

app_fmt_sci_math <- function(x, digits = 3) {
  numeric_x <- as.numeric(x)
  out <- rep("n/a", length(numeric_x))
  finite <- is.finite(numeric_x)
  parts <- strsplit(formatC(numeric_x[finite], format = "e", digits = digits), "e", fixed = TRUE)
  out[finite] <- vapply(parts, function(part) {
    sprintf("\\(%s\\times 10^{%d}\\)", part[[1]], as.integer(part[[2]]))
  }, character(1))
  out[!finite & !is.na(x)] <- as.character(x[!finite & !is.na(x)])
  out
}

app_params_table <- function(params) {
  data.frame(
    parameter = c("b", "delta", "r", "epsilon", "sigma", "mu", "c_U", "c_D"),
    value = app_fmt_num(c(
      params$b, params$delta, params$r, params$eps,
      params$sigma, params$mu, params$u, params$l
    )),
    stringsAsFactors = FALSE
  )
}

app_solution_table <- function(sol) {
  xs <- unlist(sol$x)[c("xL", "xk", "xl", "xU")]
  fmax <- app_max_residual(sol)

  data.frame(
    quantity = c(
      "lower barrier underline{x}",
      "drift threshold xk",
      "intensity threshold xl",
      "upper barrier overline{x}",
      "ergodic value gamma",
      "regime",
      "chosen solver",
      "max residual"
    ),
    value = c(
      app_fmt_key_num(xs[["xL"]]),
      app_fmt_key_num(xs[["xk"]]),
      app_fmt_key_num(xs[["xl"]]),
      app_fmt_key_num(xs[["xU"]]),
      app_fmt_key_num(sol$gamma),
      if (isTRUE(sol$regime2)) "Regime 2: x^lambda > overline{x}" else "Regime 1: x^lambda <= overline{x}",
      sol$chosen_solver %||% "unknown",
      app_fmt_key_num(fmax)
    ),
    stringsAsFactors = FALSE
  )
}

app_html_table <- function(headers, rows) {
  shiny::tags$table(
    class = "table table-sm app-data-table",
    shiny::tags$thead(
      shiny::tags$tr(lapply(headers, shiny::tags$th))
    ),
    shiny::tags$tbody(
      lapply(seq_len(nrow(rows)), function(i) {
        shiny::tags$tr(
          lapply(rows[i, , drop = TRUE], function(cell) {
            shiny::tags$td(shiny::HTML(as.character(cell)))
          })
        )
      })
    )
  )
}

app_params_table_ui <- function(params) {
  rows <- data.frame(
    parameter = c("\\(b\\)", "\\(\\delta\\)", "\\(r\\)", "\\(\\epsilon\\)",
                  "\\(\\sigma\\)", "\\(\\mu\\)", "\\(c_U\\)", "\\(c_D\\)"),
    value = app_fmt_num(c(
      params$b, params$delta, params$r, params$eps,
      params$sigma, params$mu, params$u, params$l
    )),
    stringsAsFactors = FALSE
  )
  app_html_table(c("Parameter", "Value"), rows)
}

app_solution_table_ui <- function(sol) {
  xs <- unlist(sol$x)[c("xL", "xk", "xl", "xU")]
  fmax <- app_max_residual(sol)
  regime <- if (isTRUE(sol$regime2)) {
    "Regime 2: \\(x^\\lambda > \\overline{x}\\)"
  } else {
    "Regime 1: \\(x^\\lambda \\le \\overline{x}\\)"
  }

  rows <- data.frame(
    quantity = c("\\(\\underline{x}\\)", "\\(x^\\kappa\\)", "\\(x^\\lambda\\)",
                 "\\(\\overline{x}\\)", "\\(\\gamma\\)", "Regime", "Chosen solver",
                 "Max residual"),
    value = c(
      app_fmt_key_num(xs[["xL"]]),
      app_fmt_key_num(xs[["xk"]]),
      app_fmt_key_num(xs[["xl"]]),
      app_fmt_key_num(xs[["xU"]]),
      app_fmt_key_num(sol$gamma),
      regime,
      sol$chosen_solver %||% "unknown",
      app_fmt_key_num(fmax)
    ),
    stringsAsFactors = FALSE
  )
  app_html_table(c("Quantity", "Value"), rows)
}

app_diagnostic_check_labels <- function(checks, math = FALSE) {
  if (isTRUE(math)) {
    labels <- c(
      "continuity@xL" = "\\(H(\\underline{x})=-c_U\\)",
      "continuity@xk" = "\\(H_1(x^\\kappa)=H_2(x^\\kappa)\\)",
      "continuity@xl" = "\\(H_2(x^\\lambda)=H_3(x^\\lambda)\\)",
      "continuity@xU" = "\\(H(\\overline{x})=c_D\\)",
      "H(xk)=0" = "\\(H(x^\\kappa)=0\\)",
      "IH(xl)=0" = "\\(\\mathcal{I}H(x^\\lambda)=0\\)",
      "H'(jump)@xL" = "\\(H'(\\underline{x})=0\\)",
      "H'(jump)@xk" = "\\(H_2'(x^\\kappa)=H_1'(x^\\kappa)\\)",
      "H'(jump)@xl" = "\\(H_3'(x^\\lambda)=H_2'(x^\\lambda)\\)",
      "H'(jump)@xU" = "\\(H'(\\overline{x})=0\\)"
    )
  } else {
    labels <- c(
      "continuity@xL" = "H(underline{x})=-c_U",
      "continuity@xk" = "H1(x^kappa)=H2(x^kappa)",
      "continuity@xl" = "H2(x^lambda)=H3(x^lambda)",
      "continuity@xU" = "H(overline{x})=c_D",
      "H(xk)=0" = "H(x^kappa)=0",
      "IH(xl)=0" = "IH(x^lambda)=0",
      "H'(jump)@xL" = "H'(underline{x})=0",
      "H'(jump)@xk" = "H2'(x^kappa)=H1'(x^kappa)",
      "H'(jump)@xl" = "H3'(x^lambda)=H2'(x^lambda)",
      "H'(jump)@xU" = "H'(overline{x})=0"
    )
  }
  out <- unname(labels[checks])
  out[is.na(out)] <- checks[is.na(out)]
  out
}

app_diagnostic_condition_labels <- function(types, checks = NULL) {
  labels <- c(
    continuity = "Continuity",
    ambiguity_condition = "Switch",
    derivative_jump = "Smooth fit"
  )
  out <- unname(labels[types])
  if (!is.null(checks)) {
    out[types == "ambiguity_condition" & checks == "H(xk)=0"] <- "Drift switch"
    out[types == "ambiguity_condition" & checks == "IH(xl)=0"] <- "Intensity switch"
  }
  out[is.na(out)] <- types[is.na(out)]
  out
}

app_diagnostic_pass_labels <- function(pass) {
  ifelse(is.na(pass), "n/a", ifelse(pass, "pass", "fail"))
}

app_diagnostic_pass_html <- function(pass) {
  ifelse(
    is.na(pass),
    "<span class='diagnostic-pass diagnostic-pass-na'>n/a</span>",
    ifelse(
      pass,
      "<span class='diagnostic-pass diagnostic-pass-ok' role='img' aria-label='pass'>&#10003;</span>",
      "<span class='diagnostic-pass diagnostic-pass-fail' role='img' aria-label='fail'>&#10007;</span>"
    )
  )
}

app_diagnostic_checks_table <- function(sol) {
  diag <- diagnose(sol)
  checks <- diag$checks
  data.frame(
    Check = app_diagnostic_check_labels(checks$check, math = FALSE),
    Value = app_fmt_sci(checks$value, digits = 3),
    Condition = app_diagnostic_condition_labels(checks$type, checks$check),
    Pass = app_diagnostic_pass_labels(checks$pass),
    stringsAsFactors = FALSE
  )
}

app_diagnostic_checks_table_ui <- function(sol) {
  diag <- diagnose(sol)
  checks <- diag$checks
  rows <- data.frame(
    Check = app_diagnostic_check_labels(checks$check, math = TRUE),
    Value = app_fmt_sci_math(checks$value, digits = 3),
    Condition = app_diagnostic_condition_labels(checks$type, checks$check),
    Pass = app_diagnostic_pass_html(checks$pass),
    stringsAsFactors = FALSE
  )
  app_html_table(c("Check", "Value", "Condition", "Pass"), rows)
}

app_sweep_summary_table <- function(sweep_obj) {
  res <- sweep_obj$results
  n_total <- nrow(res)
  n_ok <- sum(isTRUE(res$converged) | (is.logical(res$converged) & res$converged), na.rm = TRUE)
  data.frame(
    metric = c("sweep parameter", "grid points", "successful solves", "failed solves"),
    value = c(
      unique(res$sweep_param)[1],
      n_total,
      n_ok,
      n_total - n_ok
    ),
    stringsAsFactors = FALSE
  )
}

app_sweep_results_table <- function(sweep_obj) {
  res <- sweep_obj$results
  numeric_cols <- intersect(c("sweep_value", "xL", "xk", "xl", "xU", "gamma"), names(res))
  res[numeric_cols] <- lapply(res[numeric_cols], app_fmt_num, digits = 7)
  names(res)[names(res) == "xL"] <- "underline{x}"
  names(res)[names(res) == "xU"] <- "overline{x}"
  res
}

app_sweep_results_table_ui <- function(sweep_obj) {
  res <- app_sweep_results_table(sweep_obj)
  headers <- names(res)
  headers[headers == "underline{x}"] <- "\\(\\underline{x}\\)"
  headers[headers == "xk"] <- "\\(x^\\kappa\\)"
  headers[headers == "xl"] <- "\\(x^\\lambda\\)"
  headers[headers == "overline{x}"] <- "\\(\\overline{x}\\)"
  headers[headers == "gamma"] <- "\\(\\gamma\\)"
  headers[headers == "sweep_value"] <- "Sweep value"
  headers[headers == "sweep_param"] <- "Sweep parameter"
  app_html_table(headers, res)
}
