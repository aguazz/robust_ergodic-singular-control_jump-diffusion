# Experiments Reference

This note explains the experiment helpers under `src/experiments/`. These
functions are intentionally lighter than R package documentation, but they give
the main arguments, outputs, and intended workflow.

For normal use, read the files in this order:

1. `src/experiments/sweeps.R`
2. `src/experiments/misspecification.R`
3. `src/experiments/misspecification_plotting.R`

The scripts in `run/` call these helpers through `src/load.R`.

## Notation

- `delta`: drift-ambiguity radius.
- `eps`: jump-intensity ambiguity radius.
- `xL`, `xU`: lower and upper reflecting barriers.
- `xk`, `xl`: drift- and intensity-ambiguity thresholds.
- `gamma`: ergodic value returned by the solver.
- `RMC`: robust misspecification cost, reported as `relative_cost_pct`.

## Sweep Helpers

### `comparative_sweeper(...)`

Solves the optimal-barrier problem over a one-dimensional parameter grid.

Important arguments:

- `sweep_param`: one of `b`, `delta`, `r`, `eps`, `sigma`, `mu`, `inv_mu`,
  `u`, or `l`. `inv_mu` sweeps `1 / mu` but converts back to `mu` before
  calling `make_params()`.
- `sweep_values`: numeric vector of parameter values.
- `b`, `delta`, `r`, `eps`, `sigma`, `mu`, `u`, `l`: baseline parameters;
  all except `sweep_param` are held fixed.
- `xL0`, `xk0`, `xl0`, `xU0`: initial guesses for the first solve.
- `solver_args`: optional list of extra arguments passed to
  `solve_optimal_barriers()`.
- `save`: if `TRUE`, write the sweep table as CSV under `out_dir`.

Returns:

A list with `results`, `out_dir`, `out_name`, and `fixed_params`. `results` is a
data frame with one row per sweep value and columns for `xL`, `xk`, `xl`, `xU`,
`gamma`, convergence status, and any error message.

### `plot_sweep(sweep_obj, ...)`

Plots the barriers and ambiguity thresholds from `comparative_sweeper()`.

Useful plot controls:

- `show`, `save`, `out_dir`, `out_name`, `width_in`, `height_in`, `dpi`.
- `show_x_axis_title`, `show_y_axis_title`, `x_axis_title`, `y_axis_title`.
- `show_tick_labels`, `tick_cex`, `axis_title_cex`, `title_cex`.
- `margins`, `outer_margins`, `axis_mgp`, `tcl`, `base_cex`, `mex`.
- `plot_gamma`: if `TRUE`, include the ergodic value `gamma`.
- `gamma_layout`: `"stacked"` or `"separate"`.

Defaults preserve the visual style used before the split.

## Misspecification Grids

### `misspecification_cost_grid(delta_values, eps_values, ...)`

Computes the robust misspecification-cost grid. At each `(delta, eps)` point it
solves:

- the robust policy optimized for that ambiguous model;
- the certainty-barrier policy evaluated under the same ambiguous model.

The function warm-starts along each row of the grid, which is important for
long runs.

Returns:

An object of class `misspecification_cost_grid`, including:

- `gamma_with_ambiguity`
- `gamma_zero_ambiguity_policy`
- `gamma_difference`
- `relative_cost`
- `relative_cost_pct`
- convergence matrices and a long-form `long_results` data frame

### Cache and Table Helpers

- `save_misspecification_cost_rdata(grid_obj, ...)`: saves a grid as `.RData`.
- `load_misspecification_cost_rdata(...)`: reloads a cached grid.
- `export_misspecification_cost_latex(grid_obj, ...)`: exports a LaTeX matrix
  for a selected surface.

## Misspecification Plots

The public plotting helpers in `misspecification_plotting.R` are:

- `plot_misspecification_cost_3d()`
- `plot_misspecification_cost_contour()`
- `plot_misspecification_cost_vs_delta()`
- `plot_misspecification_cost_vs_epsilon()`
- `plot_misspecification_cost_canvas()`

These functions share the same plot-control style as `src/plotting.R`:
`show`, `save`, device size, axis-title visibility, tick sizes, margins, grid
controls, and optional legends. The two slice plots are wrappers for reading the
same grid in different directions: one fixes epsilon values and varies `delta`;
the other fixes delta values and varies `eps`.

## Reproduction Wrappers

`misspec_example_specs` stores the two paper-example configurations currently
used by `run/reproduce_misspecification.R`.

`run_misspecification_cost_example(example_spec, force = FALSE)` generates or
loads the cached grids, draws all misspecification figures, exports the LaTeX
table, and returns the paths/results invisibly.

Use `force = TRUE` only when the model parameters or grid definitions have
changed and cached `.RData` files should be rebuilt.
