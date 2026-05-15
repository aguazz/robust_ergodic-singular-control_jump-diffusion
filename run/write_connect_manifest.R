.script_path <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(normalizePath(ofile, winslash = "/", mustWork = TRUE))
    }
  }
  normalizePath(file.path("run", "write_connect_manifest.R"), winslash = "/", mustWork = TRUE)
}

required_packages <- c("shiny", "bslib", "nleqslv", "rsconnect")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    sprintf("Install missing package(s): %s", paste(missing_packages, collapse = ", ")),
    call. = FALSE
  )
}

repo_root <- normalizePath(file.path(dirname(.script_path()), ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

app_files <- c(
  "app.R",
  "app/app.R",
  "app/www/styles.css",
  list.files("src", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
)
app_files <- sort(unique(gsub("\\\\", "/", app_files)))

rsconnect::writeManifest(
  appDir = ".",
  appFiles = app_files,
  appPrimaryDoc = "app.R",
  appMode = "shiny"
)

message(sprintf("Wrote manifest.json with %d deployment files.", length(app_files)))
