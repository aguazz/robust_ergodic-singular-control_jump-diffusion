.script_path <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(normalizePath(ofile, winslash = "/", mustWork = TRUE))
    }
  }
  normalizePath(file.path("run", "launch_app.R"), winslash = "/", mustWork = TRUE)
}

.parse_port <- function(value, source) {
  port <- suppressWarnings(as.integer(value))
  if (length(port) != 1L || is.na(port) || port < 0L || port > 65535L) {
    stop(sprintf("%s must be an integer port between 0 and 65535.", source), call. = FALSE)
  }
  port
}

repo_root <- normalizePath(file.path(dirname(.script_path()), ".."), winslash = "/", mustWork = TRUE)
app_dir <- normalizePath(file.path(repo_root, "app"), winslash = "/", mustWork = TRUE)
setwd(repo_root)

host <- Sys.getenv("SHINY_HOST", "127.0.0.1")
port_env <- Sys.getenv("SHINY_PORT", "")
option_port <- getOption("shiny.port", NULL)
port_is_forced <- nzchar(port_env) || !is.null(option_port)
port <- if (nzchar(port_env)) {
  .parse_port(port_env, "SHINY_PORT")
} else if (!is.null(option_port)) {
  .parse_port(option_port, "options('shiny.port')")
} else {
  3838L
}

launch_env <- tolower(Sys.getenv("SHINY_LAUNCH_BROWSER", "true"))
launch_browser <- launch_env %in% c("1", "true", "yes", "y")

.run_app <- function(port) {
  shiny::runApp(
    app_dir,
    host = host,
    port = port,
    launch.browser = launch_browser
  )
}

tryCatch(
  .run_app(port),
  error = function(e) {
    server_error <- grepl("Failed to create server|address already in use|socket", conditionMessage(e), ignore.case = TRUE)
    if (port_is_forced || !server_error) {
      stop(e)
    }

    fallback_port <- httpuv::randomPort(host = host)
    message(sprintf(
      "Port %s is unavailable on %s; retrying on free port %s. Set SHINY_PORT to force a specific port.",
      port,
      host,
      fallback_port
    ))
    .run_app(fallback_port)
  }
)
