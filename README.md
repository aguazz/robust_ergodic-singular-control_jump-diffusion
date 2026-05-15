# Robust ergodic singular control for jump diffusions

This repository contains the R code used to solve and reproduce numerical
experiments for the paper `paper.pdf`, which studies a robust ergodic singular
control problem for a jump-diffusion model with drift and jump-intensity
ambiguity.

The code is intentionally kept as a lightweight research repository rather than
an R package. The main entry point is `src/load.R`.

## Repository layout

```text
src/
  load.R          Load all project functions in dependency order.
  solver.R        Numerical solver, diagnostics, and barrier equations.
  simulation.R    Reflected jump-diffusion simulation and pathwise costs.
  plotting.R      Static figures for the value derivative, paths, and controls.
  app/            Helpers used by the Shiny companion app.
  experiments/
    sweeps.R                    Parameter sweeps and sweep plots.
    misspecification.R          Misspecification grids, cache, tables, runners.
    misspecification_plotting.R Misspecification surface and slice plots.
  animations/
    rendering.R                 Shared PNG-frame rendering helpers.
    solver_convergence.R        Solver-convergence animation frames.
    control_canvas.R            Control/state/cost animation frames.

run/
  launch_app.R                 Launch the Shiny companion app.
  smoke.R                      Fast numerical check.
  example_basic.R              Minimal solve-and-print example.
  reproduce_figures.R          Rebuild the main static figures.
  reproduce_sweeps.R           Rebuild the parameter-sweep outputs.
  reproduce_misspecification.R Rebuild misspecification surfaces and tables.
  render_animations.R          Render animation frames under frames/.

figures/       Existing generated figures and cached grids.
frames/        Generated animation frames.
docs/          Lightweight function references and notes.
app/           Shiny companion app.
app.R          Top-level Shiny app entry point.
legacy/        Historical split source files kept for reference.
```

## Requirements

The core scripts use base R plus `nleqslv` for nonlinear systems. The optional
companion app also uses `shiny` and `bslib`.

Install missing packages with:

```r
install.packages(c("nleqslv", "shiny", "bslib"))
```

## Quick start

From the repository root:

```sh
Rscript run/smoke.R
```

or work directly inside R:

```r
source(file.path("src", "load.R"))

p <- make_params(
  b = 0, delta = 1, r = 1, eps = 0.5,
  sigma = 1, mu = 1, u = 1, l = 1
)

sol <- solve_optimal_barriers(
  p, xL0 = -0.5, xk0 = -0.1, xl0 = 0.5, xU0 = 1,
  switch_to_newton = FALSE,
  tol_regime_switch = 1e-3
)

unlist(sol$x)
sol$gamma
```

The smoke check currently verifies the baseline solution
`xL = -0.9900086`, `xk = -0.2214323`, `xl = 0.5642733`,
`xU = 0.7511408`, and `gamma = 2.48011711465`.

## Shiny companion app

Launch the interactive app from the repository root with:

```sh
Rscript run/launch_app.R
```

The app lets users choose model parameters, solve the free-boundary problem,
inspect the optimal barriers and ambiguity thresholds, simulate reflected
paths, and run modest comparative statics. See
[`docs/shiny-app.md`](docs/shiny-app.md) for details.

For compact guides to the main APIs, see
[`docs/solver-reference.md`](docs/solver-reference.md),
[`docs/simulation-reference.md`](docs/simulation-reference.md), and
[`docs/plotting-reference.md`](docs/plotting-reference.md), plus
[`docs/experiments-reference.md`](docs/experiments-reference.md) for sweeps and
misspecification grids and [`docs/animations-reference.md`](docs/animations-reference.md)
for animation-frame helpers.

## Reproducing outputs

Static figures:

```sh
Rscript run/reproduce_figures.R
```

Misspecification figures and tables:

```sh
Rscript run/reproduce_misspecification.R
```

Full parameter sweeps:

```sh
Rscript run/reproduce_sweeps.R
```

The sweep and misspecification scripts can be computationally expensive. They
are now explicit scripts rather than top-level code that runs when the solver is
loaded.

Animation frames:

```sh
Rscript run/render_animations.R
```

Frames are written under `frames/`.
