# Load the project modules in dependency order.

.loader_path <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(normalizePath(ofile, winslash = "/", mustWork = TRUE))
    }
  }
  normalizePath(file.path("src", "load.R"), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(
  file.path(dirname(.loader_path()), ".."),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(repo_root, "src", "solver.R"), chdir = TRUE)
source(file.path(repo_root, "src", "simulation.R"), chdir = TRUE)
source(file.path(repo_root, "src", "plotting.R"), chdir = TRUE)
source(file.path(repo_root, "src", "experiments.R"), chdir = TRUE)
source(file.path(repo_root, "src", "animations.R"), chdir = TRUE)

invisible(repo_root)
