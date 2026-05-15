# Shiny App Reference

The companion Shiny app lives under `app/` and is launched through
`run/launch_app.R`. It is a thin interactive layer over the existing research
code; the numerical solver, simulation, plotting, and sweep helpers remain in
`src/`.

## Launch

From the repository root:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" run\launch_app.R
```

or inside R:

```r
shiny::runApp("app")
```

Required packages are `nleqslv`, `shiny`, and `bslib`.

## Structure

```text
app/
  app.R             Shiny UI and server.
  www/styles.css    App-specific visual styling.

src/app/
  validation.R      Input defaults and validation helpers.
  solve_wrappers.R  Solver, simulation, sweep, and cache wrappers.
  plot_wrappers.R   App-specific plot wrappers.
  tables.R          Display/download table formatting.
```

The app sources `src/load.R` and then the app helper files. `src/load.R` is not
modified to load Shiny-specific code, so script workflows stay lightweight.

## Main Panels

- **Solution** solves the free-boundary problem and reports barriers,
  ambiguity thresholds, regime, ergodic value, residuals, `H`, `H'`, and the
  worst-case drift/intensity distortions.
- **Simulation** simulates the reflected jump-diffusion path under the solved
  policy and plots the state, singular controls, and running-average cost.
- **Comparative Statics** runs small explicit one-dimensional sweeps using
  `comparative_sweeper()`.
- **Diagnostics** exposes `diagnose(sol)` and the optional recorded solver path.

## Notes

The app intentionally uses a button-triggered solve. Changing a parameter does
not recompute the free-boundary system until the user clicks **Solve**. Sweeps
are also explicit and limited to modest grid sizes because full paper-scale
experiments are computationally expensive.
