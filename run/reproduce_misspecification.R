source(file.path("src", "load.R"))

runs <- lapply(
  misspec_example_specs,
  run_misspecification_cost_example
)

names(runs) <- vapply(
  runs,
  function(run) run$spec$tag,
  character(1)
)

cat(paste(
  vapply(runs, function(run) run$result$tex$latex, character(1)),
  collapse = "\n\n"
))
