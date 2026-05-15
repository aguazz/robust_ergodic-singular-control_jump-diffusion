required_packages <- c("shiny", "bslib", "nleqslv")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    sprintf("Install missing package(s): %s", paste(missing_packages, collapse = ", ")),
    call. = FALSE
  )
}

.app_file_path <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(normalizePath(ofile, winslash = "/", mustWork = TRUE))
    }
  }
  if (file.exists("app.R")) {
    return(normalizePath("app.R", winslash = "/", mustWork = TRUE))
  }
  normalizePath(file.path("app", "app.R"), winslash = "/", mustWork = TRUE)
}

app_dir <- dirname(.app_file_path())
app_repo_root <- normalizePath(file.path(app_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(app_repo_root)
shiny::addResourcePath("app-www", file.path(app_dir, "www"))

source(file.path(app_repo_root, "src", "load.R"), chdir = TRUE)
source(file.path(app_repo_root, "src", "app", "validation.R"), chdir = TRUE)
source(file.path(app_repo_root, "src", "app", "solve_wrappers.R"), chdir = TRUE)
source(file.path(app_repo_root, "src", "app", "plot_wrappers.R"), chdir = TRUE)
source(file.path(app_repo_root, "src", "app", "tables.R"), chdir = TRUE)

model_defaults <- app_default_model_params()
guess_defaults <- app_default_initial_guess()
solver_defaults <- app_default_solver_options()
sim_defaults <- app_default_simulation_options()

sweep_default_ranges <- list(
  b = c(-5, 5),
  delta = c(0, 5),
  r = c(0.05, 10),
  eps = c(0, 1),
  sigma = c(0.1, 5),
  mu = c(0.2, 2),
  inv_mu = c(0.05, 3),
  u = c(0.1, 5),
  l = c(0.1, 5)
)

theme <- bslib::bs_theme(
  version = 5,
  primary = "#2F6F73",
  secondary = "#6B7280",
  success = "#228B22",
  danger = "#B22222",
  base_font = bslib::font_google("Source Sans 3"),
  code_font = bslib::font_google("Source Code Pro")
)

app_math <- function(tex) {
  shiny::span(shiny::HTML(sprintf("\\(%s\\)", tex)))
}

app_section_help <- function(label, help) {
  shiny::span(
    class = "section-title-with-help",
    shiny::span(label),
    shiny::span(
      class = "section-info-icon",
      role = "button",
      tabindex = "0",
      `data-bs-toggle` = "popover",
      `data-bs-trigger` = "hover focus",
      `data-bs-placement` = "right",
      `data-bs-html` = "true",
      `data-bs-title` = paste(label, "parameters"),
      `data-bs-content` = help,
      `aria-label` = paste(label, "information"),
      "i"
    )
  )
}

model_help <- paste0(
  "<ul class='section-help-list'>",
  "<li>\\(b\\): baseline drift.</li>",
  "<li>\\(\\delta\\): drift-ambiguity radius.</li>",
  "<li>\\(r\\): reference jump intensity.</li>",
  "<li>\\(\\epsilon\\): intensity-ambiguity radius in \\([0,1]\\).</li>",
  "<li>\\(\\sigma\\): diffusion volatility.</li>",
  "<li>\\(\\mu\\): negative exponential jump-size rate, with \\(E[Y]=-1/\\mu\\).</li>",
  "<li>\\(c_U,c_D\\): upward and downward intervention costs.</li>",
  "</ul>"
)

solver_help <- paste0(
  "<ul class='section-help-list'>",
  "<li>\\(\\underline{x}_0,x^\\kappa_0,x^\\lambda_0,\\overline{x}_0\\): initial guesses for \\(\\underline{x},x^\\kappa,x^\\lambda,\\overline{x}\\).</li>",
  "<li>\\(f_{\\mathrm{tol}},x_{\\mathrm{tol}}\\): nonlinear-solver tolerances.</li>",
  "<li>\\(\\mathrm{maxit},\\mathrm{stepmax}\\): limits for the Broyden/Newton root search.</li>",
  "<li>\\(\\mathrm{tol}_{\\mathrm{regime}}\\): near-boundary tolerance for classifying Regime 2.</li>",
  "</ul>"
)

simulation_help <- paste0(
  "<ul class='section-help-list'>",
  "<li>\\(T\\): final simulated time.</li>",
  "<li>\\(\\Delta t\\): Euler time step.</li>",
  "<li>Seed: reproducibility control.</li>",
  "<li>Max steps: safety cap for long simulations.</li>",
  "<li>\\(x_0\\): initial state, used when midpoint start is disabled.</li>",
  "</ul>"
)

app_number_text <- function(inputId, label, value, help = NULL) {
  shiny::textInput(
    inputId,
    label = label,
    value = format(value, scientific = FALSE, trim = TRUE)
  )
}

app_plot_download <- function(pngOutputId, pdfOutputId = paste0(pngOutputId, "_pdf")) {
  shiny::div(
    class = "plot-actions",
    shiny::downloadButton(
      pngOutputId,
      "Download PNG",
      class = "btn-sm"
    ),
    shiny::downloadButton(
      pdfOutputId,
      "Download PDF",
      class = "btn-sm"
    )
  )
}

parameter_grid <- shiny::div(
  class = "parameter-grid",
  app_number_text("b", app_math("b"), model_defaults$b),
  app_number_text("delta", app_math("\\delta"), model_defaults$delta),
  app_number_text("r", app_math("r"), model_defaults$r),
  app_number_text("eps", app_math("\\epsilon"), model_defaults$eps),
  app_number_text("sigma", app_math("\\sigma"), model_defaults$sigma),
  app_number_text("mu", app_math("\\mu"), model_defaults$mu),
  app_number_text("u", app_math("c_U"), model_defaults$u),
  app_number_text("l", app_math("c_D"), model_defaults$l)
)

solver_controls <- bslib::accordion(
  open = FALSE,
  bslib::accordion_panel(
    title = app_section_help("Solver", solver_help),
    value = "solver",
    shiny::div(
      class = "parameter-grid",
      app_number_text("xL0", app_math("\\underline{x}_0"), guess_defaults$xL0),
      app_number_text("xk0", app_math("x_{\\kappa,0}"), guess_defaults$xk0),
      app_number_text("xl0", app_math("x_{\\lambda,0}"), guess_defaults$xl0),
      app_number_text("xU0", app_math("\\overline{x}_0"), guess_defaults$xU0),
      app_number_text("ftol", app_math("f_{\\mathrm{tol}}"), solver_defaults$ftol),
      app_number_text("xtol", app_math("x_{\\mathrm{tol}}"), solver_defaults$xtol),
      app_number_text("maxit", app_math("\\mathrm{maxit}"), solver_defaults$maxit),
      app_number_text("stepmax", app_math("\\mathrm{stepmax}"), solver_defaults$stepmax),
      app_number_text(
        "tol_regime_switch",
        app_math("\\mathrm{tol}_{\\mathrm{regime}}"),
        solver_defaults$tol_regime_switch
      )
    ),
    shiny::div(
      class = "section-checkboxes",
      shiny::checkboxInput("switch_to_newton", "Newton refinement", solver_defaults$switch_to_newton),
      shiny::checkboxInput("record_iterates", "Record iterates", solver_defaults$record_iterates)
    )
  ),
  bslib::accordion_panel(
    title = app_section_help("Simulation", simulation_help),
    value = "simulation",
    shiny::div(
      class = "parameter-grid",
      app_number_text("sim_T", app_math("T"), sim_defaults$T),
      app_number_text("sim_dt", app_math("\\Delta t"), sim_defaults$dt),
      app_number_text("sim_seed", "seed", sim_defaults$seed),
      app_number_text("sim_max_steps", "max steps", sim_defaults$max_steps)
    ),
    shiny::div(
      class = "section-checkboxes",
      shiny::checkboxInput("sim_midpoint", "Start at midpoint", TRUE),
      shiny::conditionalPanel(
        condition = "!input.sim_midpoint",
        app_number_text("sim_x0", app_math("x_0"), 0)
      )
    ),
    shiny::actionButton(
      "simulate",
      shiny::tagList(shiny::icon("wave-square"), "Simulate"),
      class = "btn-outline-primary full-width"
    )
  )
)

ui <- bslib::page_sidebar(
  title = shiny::div(
    class = "title-wrap",
    shiny::span(class = "title-main", "Robust Ergodic Control of Jump-Diffusion Systems"),
    shiny::span(class = "title-sub", "under Drift and Intensity Uncertainty - Companion App")
  ),
  theme = theme,
  fillable = FALSE,
  shiny::withMathJax(),
  shiny::tags$head(
    shiny::tags$link(rel = "stylesheet", type = "text/css", href = "app-www/styles.css"),
    shiny::tags$script(shiny::HTML(
      "
      document.addEventListener('DOMContentLoaded', function() {
        if (window.bootstrap && bootstrap.Popover) {
          document.querySelectorAll('[data-bs-toggle=\"popover\"]').forEach(function(el) {
            if (!bootstrap.Popover.getInstance(el)) {
              new bootstrap.Popover(el);
              el.addEventListener('shown.bs.popover', function() {
                if (window.MathJax && window.MathJax.typesetPromise) {
                  window.MathJax.typesetPromise();
                } else if (window.MathJax && window.MathJax.Hub && window.MathJax.Hub.Queue) {
                  window.MathJax.Hub.Queue(['Typeset', window.MathJax.Hub]);
                }
              });
            }
          });
        }

        var solveOnEnterIds = new Set([
          'b', 'delta', 'r', 'eps', 'sigma', 'mu', 'u', 'l',
          'xL0', 'xk0', 'xl0', 'xU0',
          'ftol', 'xtol', 'maxit', 'stepmax', 'tol_regime_switch',
          'sim_T', 'sim_dt', 'sim_seed', 'sim_max_steps', 'sim_x0'
        ]);

        document.addEventListener('keydown', function(event) {
          if (event.key !== 'Enter' || event.defaultPrevented || event.isComposing) return;
          var target = event.target;
          if (!target || !solveOnEnterIds.has(target.id)) return;
          var solveButton = document.getElementById('solve');
          if (!solveButton || solveButton.disabled) return;
          event.preventDefault();
          solveButton.click();
        });
      });
      "
    ))
  ),
  sidebar = bslib::sidebar(
    width = 350,
    shiny::div(class = "sidebar-heading", app_section_help("Model", model_help)),
    parameter_grid,
    shiny::actionButton(
      "solve",
      shiny::tagList(shiny::icon("calculator"), "Solve"),
      class = "btn-primary full-width solve-button"
    ),
    solver_controls
  ),
  shiny::div(
    class = "app-content",
  bslib::navset_tab(
    id = "main_tabs",
    bslib::nav_panel(
      "Solution",
      shiny::uiOutput("solution_status"),
      bslib::layout_columns(
        col_widths = c(4, 8),
        bslib::card(
          class = "stacked-table-card",
          bslib::card_header("Key quantities"),
          shiny::uiOutput("solution_table"),
          shiny::div(
            class = "card-bottom-actions",
            shiny::downloadButton("download_solution", "Results CSV", class = "btn-sm")
          )
        ),
        bslib::card(
          class = "stacked-table-card",
          bslib::card_header("Residual checks"),
          shiny::uiOutput("diagnostic_table"),
          shiny::div(
            class = "card-bottom-actions",
            shiny::downloadButton("download_diagnostics", "Diagnostics CSV", class = "btn-sm")
          )
        )
      ),
      bslib::card(
        bslib::card_header("Value diagnostics"),
        bslib::navset_tab(
          id = "solution_value_tabs",
          bslib::nav_panel(
            title = shiny::tagList("Marginal value ", app_math("H")),
            value = "marginal_value",
            shiny::plotOutput("plot_H", height = 500),
            app_plot_download("download_plot_H")
          ),
          bslib::nav_panel(
            title = shiny::tagList("Derivative diagnostic ", app_math("H'(x)")),
            value = "derivative_diagnostic",
            shiny::plotOutput("plot_H_prime", height = 500),
            app_plot_download("download_plot_H_prime")
          ),
          bslib::nav_panel(
            title = "Recorded solver path",
            value = "solver_path",
            shiny::plotOutput("plot_solver_history", height = 500),
            app_plot_download("download_plot_solver_history")
          )
        )
      )
    ),
    bslib::nav_panel(
      "Simulation",
      shiny::uiOutput("simulation_status"),
      bslib::card(
        bslib::card_header("Simulation diagnostics"),
        bslib::navset_tab(
          id = "simulation_plot_tabs",
          bslib::nav_panel(
            title = shiny::tagList("Reflected state and singular controls ", app_math("X_t,U_t,D_t")),
            value = "reflected_controls",
            shiny::plotOutput("plot_simulation", height = 720),
            app_plot_download("download_plot_simulation")
          ),
          bslib::nav_panel(
            title = shiny::tagList("Running-average cost ", app_math("J_t")),
            value = "running_cost",
            shiny::plotOutput("plot_cost", height = 460),
            app_plot_download("download_plot_cost")
          )
        )
      )
    ),
    bslib::nav_panel(
      "Comparative Statics",
      bslib::card(
        bslib::card_header("Sweep setup"),
        bslib::layout_columns(
          col_widths = c(3, 2, 2, 2, 3),
          shiny::selectInput(
            "sweep_param",
            "parameter",
            choices = c(
              "b" = "b",
              "delta" = "delta",
              "r" = "r",
              "epsilon" = "eps",
              "sigma" = "sigma",
              "mu" = "mu",
              "1 / mu" = "inv_mu",
              "c_U" = "u",
              "c_D" = "l"
            ),
            selected = "b"
          ),
          app_number_text("sweep_from", "from", sweep_default_ranges$b[1], "Lower endpoint of the one-dimensional parameter sweep."),
          app_number_text("sweep_to", "to", sweep_default_ranges$b[2], "Upper endpoint of the one-dimensional parameter sweep."),
          app_number_text("sweep_n", "points", 25, "Number of grid points in the sweep."),
          shiny::div(
            class = "sweep-actions",
            shiny::checkboxInput("sweep_gamma", "Include gamma", TRUE),
            shiny::actionButton(
              "run_sweep",
              shiny::tagList(shiny::icon("chart-line"), "Run sweep"),
              class = "btn-primary"
            )
          )
        )
      ),
      shiny::uiOutput("sweep_status"),
      bslib::card(
        bslib::card_header(shiny::tagList("Thresholds and value ", app_math("\\underline{x},x^\\kappa,x^\\lambda,\\overline{x},\\gamma"))),
        shiny::plotOutput("plot_sweep", height = 700),
        shiny::div(
          class = "plot-actions",
          shiny::downloadButton("download_plot_sweep", "Download PNG", class = "btn-sm"),
          shiny::downloadButton("download_plot_sweep_pdf", "Download PDF", class = "btn-sm"),
          shiny::downloadButton("download_sweep", "Download data", class = "btn-sm")
        )
      )
    )
  )
  )
)

server <- function(input, output, session) {
  write_plot_png <- function(file, width = 1400, height = 900, res = 144, expr) {
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
    force(expr)
  }

  write_plot_pdf <- function(file, width = 9, height = 6, expr) {
    grDevices::pdf(file, width = width, height = height, onefile = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    force(expr)
  }

  current_params <- function() {
    app_make_params(
      b = input$b,
      delta = input$delta,
      r = input$r,
      eps = input$eps,
      sigma = input$sigma,
      mu = input$mu,
      u = input$u,
      l = input$l
    )
  }

  current_guesses <- function() {
    app_validate_initial_guess(input$xL0, input$xk0, input$xl0, input$xU0)
  }

  current_options <- function(record_iterates = input$record_iterates) {
    app_validate_solver_options(
      ftol = input$ftol,
      xtol = input$xtol,
      maxit = input$maxit,
      stepmax = input$stepmax,
      tol_regime_switch = input$tol_regime_switch,
      switch_to_newton = input$switch_to_newton,
      record_iterates = record_iterates
    )
  }

  solution_result <- shiny::eventReactive(
    input$solve,
    {
      tryCatch(
        shiny::withProgress(message = "Solving free-boundary system", value = 0, {
          params <- current_params()
          guesses <- current_guesses()
          options <- current_options()
          shiny::incProgress(0.25, detail = "Building model")
          sol <- app_solve_solution(params, guesses, options)
          shiny::incProgress(0.65, detail = "Checking diagnostics")
          diag <- diagnose(sol)
          shiny::incProgress(0.10)
          list(
            ok = TRUE,
            params = params,
            guesses = guesses,
            options = options,
            sol = sol,
            diagnostics = diag,
            message = app_solver_message(sol)
          )
        }),
        error = function(e) {
          list(ok = FALSE, message = conditionMessage(e))
        }
      )
    },
    ignoreNULL = FALSE
  )

  solution_ok <- shiny::reactive({
    res <- solution_result()
    shiny::validate(shiny::need(isTRUE(res$ok), res$message))
    res
  })

  simulation_result <- shiny::eventReactive(
    list(solution_result(), input$simulate),
    {
      res <- solution_result()
      if (!isTRUE(res$ok)) return(list(ok = FALSE, message = res$message))

      tryCatch(
        shiny::withProgress(message = "Simulating reflected path", value = 0, {
          sim_options <- app_validate_simulation_options(
            T = input$sim_T,
            dt = input$sim_dt,
            seed = input$sim_seed,
            use_midpoint = input$sim_midpoint,
            x0 = input$sim_x0,
            max_steps = input$sim_max_steps
          )
          shiny::incProgress(0.35, detail = "Drawing shocks")
          sim <- app_simulate_solution(res$params, res$sol, sim_options)
          shiny::incProgress(0.45, detail = "Computing pathwise cost")
          cost <- finite_horizon_ergodic_cost(sim, res$params)
          shiny::incProgress(0.20)
          list(ok = TRUE, sim = sim, cost = cost, options = sim_options)
        }),
        error = function(e) {
          list(ok = FALSE, message = conditionMessage(e))
        }
      )
    },
    ignoreNULL = FALSE
  )

  simulation_ok <- shiny::reactive({
    res <- simulation_result()
    shiny::validate(shiny::need(isTRUE(res$ok), res$message))
    res
  })

  output$solution_status <- shiny::renderUI({
    res <- solution_result()
    if (is.null(res)) return(NULL)
    if (!isTRUE(res$ok)) {
      return(shiny::div(class = "app-alert danger", res$message))
    }
    shiny::div(
      class = "app-alert success",
      sprintf(
        "Solved with %s; max residual %s.",
        res$sol$chosen_solver %||% "solver",
        app_fmt_sci(app_max_residual(res$sol), digits = 3)
      )
    )
  })

  output$solution_table <- shiny::renderUI({
    shiny::withMathJax(app_solution_table_ui(solution_ok()$sol))
  })

  output$plot_H <- shiny::renderPlot({
    res <- solution_ok()
    app_plot_H(res$sol, res$params)
  }, res = 96)

  output$plot_H_prime <- shiny::renderPlot({
    res <- solution_ok()
    app_plot_H_prime(res$sol, res$params)
  }, res = 96)

  output$simulation_status <- shiny::renderUI({
    sim <- simulation_result()
    if (is.null(sim)) return(NULL)
    if (!isTRUE(sim$ok)) {
      return(shiny::div(class = "app-alert danger", sim$message))
    }
    shiny::div(
      class = "app-alert success",
      sprintf(
        "Simulated %s steps with seed %s.",
        format(sim$options$n_steps, big.mark = ","),
        sim$options$seed
      )
    )
  })

  output$plot_simulation <- shiny::renderPlot({
    sol <- solution_ok()
    sim <- simulation_ok()
    app_plot_reflected_with_controls(sol$sol, sol$params, sim$sim)
  }, res = 96)

  output$plot_cost <- shiny::renderPlot({
    sol <- solution_ok()
    sim <- simulation_ok()
    app_plot_cost(sim$cost, sol$sol$gamma)
  }, res = 96)

  observeEvent(input$sweep_param, {
    rng <- sweep_default_ranges[[input$sweep_param]]
    shiny::updateTextInput(session, "sweep_from", value = rng[1])
    shiny::updateTextInput(session, "sweep_to", value = rng[2])
  }, ignoreInit = TRUE)

  sweep_result <- shiny::eventReactive(input$run_sweep, {
    tryCatch(
      shiny::withProgress(message = "Running comparative statics", value = 0, {
        sol_res <- solution_ok()
        values <- app_validate_sweep_values(input$sweep_from, input$sweep_to, input$sweep_n)
        guesses <- app_solution_as_initial_guess(sol_res$sol, sol_res$guesses)
        options <- current_options(record_iterates = FALSE)
        shiny::incProgress(0.15, detail = "Preparing grid")
        sweep <- app_run_sweep(sol_res$params, guesses, options, input$sweep_param, values)
        shiny::incProgress(0.85)
        list(ok = TRUE, sweep = sweep)
      }),
      error = function(e) {
        list(ok = FALSE, message = conditionMessage(e))
      }
    )
  }, ignoreNULL = TRUE)

  sweep_ok <- shiny::reactive({
    res <- sweep_result()
    shiny::validate(shiny::need(!is.null(res), "Run a sweep to populate this panel."))
    shiny::validate(shiny::need(isTRUE(res$ok), res$message))
    res$sweep
  })

  output$sweep_status <- shiny::renderUI({
    res <- sweep_result()
    if (is.null(res)) {
      return(shiny::div(class = "app-alert neutral", "No sweep has been run in this session."))
    }
    if (!isTRUE(res$ok)) {
      return(shiny::div(class = "app-alert danger", res$message))
    }
    ss <- app_sweep_summary_table(res$sweep)
    shiny::div(
      class = "app-alert success",
      sprintf(
        "%s successful solves over %s grid points.",
        ss$value[ss$metric == "successful solves"],
        ss$value[ss$metric == "grid points"]
      )
    )
  })

  output$plot_sweep <- shiny::renderPlot({
    app_plot_sweep_thresholds(sweep_ok(), include_gamma = input$sweep_gamma)
  }, res = 96)

  output$diagnostic_table <- shiny::renderUI({
    shiny::withMathJax(app_diagnostic_checks_table_ui(solution_ok()$sol))
  })

  output$plot_solver_history <- shiny::renderPlot({
    app_plot_solver_history(solution_ok()$sol)
  }, res = 96)

  output$download_solution <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_solution_", app_filename_tag(Sys.Date()), ".csv")
    },
    content = function(file) {
      utils::write.csv(app_solution_table(solution_ok()$sol), file, row.names = FALSE)
    }
  )

  output$download_plot_H <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_H_", app_filename_tag(Sys.Date()), ".png")
    },
    content = function(file) {
      res <- solution_ok()
      write_plot_png(file, width = 1300, height = 780, expr = {
        app_plot_H(res$sol, res$params)
      })
    }
  )

  output$download_plot_H_pdf <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_H_", app_filename_tag(Sys.Date()), ".pdf")
    },
    content = function(file) {
      res <- solution_ok()
      write_plot_pdf(file, width = 9.0, height = 5.4, expr = {
        app_plot_H(res$sol, res$params)
      })
    }
  )

  output$download_plot_H_prime <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_H_prime_", app_filename_tag(Sys.Date()), ".png")
    },
    content = function(file) {
      res <- solution_ok()
      write_plot_png(file, width = 1300, height = 780, expr = {
        app_plot_H_prime(res$sol, res$params)
      })
    }
  )

  output$download_plot_H_prime_pdf <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_H_prime_", app_filename_tag(Sys.Date()), ".pdf")
    },
    content = function(file) {
      res <- solution_ok()
      write_plot_pdf(file, width = 9.0, height = 5.4, expr = {
        app_plot_H_prime(res$sol, res$params)
      })
    }
  )

  output$download_plot_simulation <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_simulation_", app_filename_tag(Sys.Date()), ".png")
    },
    content = function(file) {
      sol <- solution_ok()
      sim <- simulation_ok()
      write_plot_png(file, width = 1400, height = 1050, expr = {
        app_plot_reflected_with_controls(sol$sol, sol$params, sim$sim)
      })
    }
  )

  output$download_plot_simulation_pdf <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_simulation_", app_filename_tag(Sys.Date()), ".pdf")
    },
    content = function(file) {
      sol <- solution_ok()
      sim <- simulation_ok()
      write_plot_pdf(file, width = 9.7, height = 7.3, expr = {
        app_plot_reflected_with_controls(sol$sol, sol$params, sim$sim)
      })
    }
  )

  output$download_plot_cost <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_running_cost_", app_filename_tag(Sys.Date()), ".png")
    },
    content = function(file) {
      sol <- solution_ok()
      sim <- simulation_ok()
      write_plot_png(file, width = 1300, height = 700, expr = {
        app_plot_cost(sim$cost, sol$sol$gamma)
      })
    }
  )

  output$download_plot_cost_pdf <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_running_cost_", app_filename_tag(Sys.Date()), ".pdf")
    },
    content = function(file) {
      sol <- solution_ok()
      sim <- simulation_ok()
      write_plot_pdf(file, width = 9.0, height = 4.9, expr = {
        app_plot_cost(sim$cost, sol$sol$gamma)
      })
    }
  )

  output$download_plot_sweep <- shiny::downloadHandler(
    filename = function() {
      sweep <- sweep_ok()
      paste0("robust_sweep_plot_", app_filename_tag(unique(sweep$results$sweep_param)[1]), ".png")
    },
    content = function(file) {
      write_plot_png(file, width = 1400, height = 1050, expr = {
        app_plot_sweep_thresholds(sweep_ok(), include_gamma = input$sweep_gamma)
      })
    }
  )

  output$download_plot_sweep_pdf <- shiny::downloadHandler(
    filename = function() {
      sweep <- sweep_ok()
      paste0("robust_sweep_plot_", app_filename_tag(unique(sweep$results$sweep_param)[1]), ".pdf")
    },
    content = function(file) {
      write_plot_pdf(file, width = 9.7, height = 7.3, expr = {
        app_plot_sweep_thresholds(sweep_ok(), include_gamma = input$sweep_gamma)
      })
    }
  )

  output$download_plot_solver_history <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_solver_history_", app_filename_tag(Sys.Date()), ".png")
    },
    content = function(file) {
      write_plot_png(file, width = 1300, height = 700, expr = {
        app_plot_solver_history(solution_ok()$sol)
      })
    }
  )

  output$download_plot_solver_history_pdf <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_solver_history_", app_filename_tag(Sys.Date()), ".pdf")
    },
    content = function(file) {
      write_plot_pdf(file, width = 9.0, height = 4.9, expr = {
        app_plot_solver_history(solution_ok()$sol)
      })
    }
  )

  output$download_diagnostics <- shiny::downloadHandler(
    filename = function() {
      paste0("robust_diagnostics_", app_filename_tag(Sys.Date()), ".csv")
    },
    content = function(file) {
      checks <- app_diagnostic_checks_table(solution_ok()$sol)
      utils::write.csv(checks, file, row.names = FALSE)
    }
  )

  output$download_sweep <- shiny::downloadHandler(
    filename = function() {
      sweep <- sweep_ok()
      paste0("robust_sweep_", app_filename_tag(unique(sweep$results$sweep_param)[1]), ".csv")
    },
    content = function(file) {
      utils::write.csv(sweep_ok()$results, file, row.names = FALSE)
    }
  )
}

shiny::shinyApp(ui, server)
