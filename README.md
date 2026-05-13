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
  experiments.R   Parameter sweeps, misspecification grids, and table export.
  animations.R    Frame-by-frame animation helpers.

run/
  smoke.R                      Fast numerical check.
  example_basic.R              Minimal solve-and-print example.
  reproduce_figures.R          Rebuild the main static figures.
  reproduce_sweeps.R           Rebuild the parameter-sweep outputs.
  reproduce_misspecification.R Rebuild misspecification surfaces and tables.
  render_animations.R          Render animation frames under frames/.

figures/       Existing generated figures and cached grids.
frames/        Generated animation frames.
legacy/        Backup copy of the original monolithic scripts.
```

## Requirements

The code uses base R plus:

- `nleqslv` for nonlinear systems.
- `magick` only for optional PNG frame optimization.

Install missing packages with:

```r
install.packages(c("nleqslv", "magick"))
```

## Quick start

From the repository root:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" run\smoke.R
```

or, inside R:

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

## Reproducing outputs

Static figures:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" run\reproduce_figures.R
```

Misspecification figures and tables:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" run\reproduce_misspecification.R
```

Full parameter sweeps:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" run\reproduce_sweeps.R
```

The sweep and misspecification scripts can be computationally expensive. They
are now explicit scripts rather than top-level code that runs when the solver is
loaded.

Animation frames:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" run\render_animations.R
```

Frames are written under `frames/`.

## Notes for future edits

- Add numerical changes in `src/solver.R`.
- Add plotting changes in `src/plotting.R` or `src/animations.R`.
- Add new experiments as opt-in scripts under `run/`.
- Keep generated files out of source modules. Loading `src/load.R` should
  define functions only.
- The original files are preserved in `legacy/` for comparison and rollback.
