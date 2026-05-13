# Plotting Reference

This note explains the static plotting helpers in `src/plotting.R`: what each
plot function draws, which objects it expects, and how to control the visual
settings without changing the defaults.

For normal use, read the functions in this order:

1. `plot_H()`
2. `plot_H_prime()`
3. `plot_reflected_jd()`
4. `plot_controls()`
5. `plot_reflected_with_controls()`

`render_and_save()` is the shared device wrapper used by these plotting
functions and by several experiment helpers.

## Notation

The plotting code uses the same notation as the solver and simulation modules:

- `H`: marginal value derivative returned by the solver.
- `Hp`: derivative of `H`.
- `X`: reflected state path.
- `U`: cumulative upward reflection at the lower barrier `xL`.
- `L`: cumulative downward reflection at the upper barrier `xU`.
- `xL`, `xU`: lower and upper reflecting barriers.
- `xk`: drift-ambiguity threshold.
- `xl`: intensity-ambiguity threshold.
- `u`, `l`: upward and downward intervention costs.

Some plot labels use `D_t` for the downward reflection process. The simulation
object itself stores this process as `L`.

## Common Plot Arguments

Most public plotting functions share these arguments:

- `show`: draw on the active graphics device. Defaults to `interactive()`.
- `save`: save PNG and PDF files. Defaults to `TRUE`.
- `out_base`: output filename without extension.
- `dir`: output directory.
- `width_in`, `height_in`, `dpi`: PNG/PDF size settings.
- `match_current`: if `TRUE`, reuse size and text settings from the active
  device.
- `margins`: passed to `par(mar = ...)`.
- `axis_mgp`: passed to `par(mgp = ...)`.
- `tcl`: tick length.
- `show_x_axis_title`, `show_y_axis_title`: control axis-title visibility.
- `x_axis_title`, `y_axis_title`: axis labels.
- `show_tick_labels`: hide or show tick labels.
- `tick_cex`: axis tick-label size.
- `axis_title_cex`: axis-title size.
- `base_cex`, `mex`: base graphics text scaling parameters.

Defaults were chosen to preserve the pre-refactor visual appearance. Some newer
controls default to `NULL` in functions that previously relied on base-R
defaults; passing a value opts into explicit styling.

## Public Functions

### `plot_H(sol, params, ...)`

Plots the solved marginal value derivative `H(x)` together with the four
thresholds/barriers and the intervention-cost levels.

Core arguments:

- `sol`: solution object from `solve_optimal_barriers()` or compatible object
  from `build_suboptimal_H()`.
- `params`: parameter list with costs `u` and `l`.
- `n`: number of grid points used to draw the curve.
- `pad_x_frac`: horizontal padding around the plotted threshold range.
- `top_blank`, `bottom_blank`: extra vertical whitespace inside the plotting
  range.

Returns:

The function is used for its plotting side effect and returns invisibly.

Example:

```r
plot_H(
  sol, p,
  show = interactive(),
  save = TRUE,
  show_x_axis_title = TRUE,
  x_axis_title = expression(x)
)
```

### `plot_H_prime(sol, params, ...)`

Plots `Hp(x)`, the derivative of `H`, with line breaks around threshold points
so derivative jumps are not connected visually.

Core arguments:

- `sol`: solution object with an `Hp` function.
- `params`: accepted for API symmetry with `plot_H()`.
- `n`, `pad_x_frac`, `top_blank`, `bottom_blank`: grid and layout controls.

Notes:

The visual style now accepts the same common controls as `plot_H()`. By default,
the function still follows the earlier base-R look unless styling arguments are
provided.

### `plot_reflected_jd(sol, params, sim = NULL, seed = 123, ...)`

Plots one reflected state path with barriers and ambiguity thresholds.

Arguments:

- `sol`: solution object with `xL`, `xk`, `xl`, and `xU` in `sol$x`.
- `params`: model parameters, used only when `sim` is `NULL`.
- `sim`: optional simulation object from `simulate_reflected_jd()`.
- `seed`: seed used if the function needs to simulate a path.
- `top_blank`, `bottom_blank`: vertical whitespace around the reflecting band.

Threshold paths:

If `sim` contains `xL`, `xU`, `xk`, or `xl` paths, those paths are drawn. If not,
the scalar values from `sol$x` are used. This keeps old behavior for static
thresholds and supports time-varying simulation output.

Example:

```r
sim <- simulate_reflected_jd(params = p, thresholds = sol$x, seed = 123)

plot_reflected_jd(
  sol, p, sim = sim,
  show = interactive(),
  save = TRUE
)
```

### `plot_controls(sim, ...)`

Plots the cumulative singular controls from a simulation object.

Arguments:

- `sim`: simulation object with `time`, `U`, and `L`.
- Common plotting arguments listed above.

Output:

- `sim$L` is drawn as the downward push and labelled `D_t` in the legend.
- `sim$U` is drawn as the upward push.

Example:

```r
plot_controls(
  sim,
  show = interactive(),
  save = TRUE,
  show_x_axis_title = TRUE,
  show_y_axis_title = TRUE
)
```

### `plot_reflected_with_controls(sol, params, sim = NULL, seed = 123, ...)`

Builds a three-panel figure:

1. cumulative downward reflection;
2. reflected state path with barriers and ambiguity thresholds;
3. cumulative upward reflection.

Arguments:

- `sol`: solution object.
- `params`: model parameters, used only when `sim` is `NULL`.
- `sim`: optional simulation object.
- `heights`: relative panel heights.
- `draw_legend`: whether to draw the threshold legend in the middle panel.
- `y_axis_title`: either a single label used for both control panels, or a list
  with `top` and `bottom` labels.
- Common plotting arguments listed above.

Example:

```r
plot_reflected_with_controls(
  sol, p, sim = sim,
  heights = c(0.7, 1.7, 0.85),
  draw_legend = FALSE,
  show_x_axis_title = TRUE,
  x_axis_title = "t",
  show_y_axis_title = FALSE
)
```

## Device Helper

### `render_and_save(fname_base, plotfun, ...)`

Draws a plot function on screen, saves it, or both.

Arguments:

- `fname_base`: output filename without extension.
- `plotfun`: zero-argument function that creates the plot.
- `show`: if `TRUE`, draw on the active device.
- `save`: if `TRUE`, save both PNG and PDF files.
- `dir`: output directory.
- `width_in`, `height_in`, `dpi`: output size for saved files.
- `pointsize`, `family`: graphics-device text settings.
- `match_current`: if `TRUE`, infer size and text settings from the active
  device.
- `use_cairo_pdf`: prefer `cairo_pdf()` when available.

This helper is intentionally public because `experiments.R` also uses it for
custom plots.

## Internal Helpers

The remaining functions are implementation details:

- `.threshold_style()`: shared colors and line types for `xL`, `xk`, `xl`,
  and `xU`.
- `.sol_thresholds()`: extracts and checks `sol$x`.
- `.expand_ylim()`: adds controlled blank space to y-limits.
- `.legend_widths()`: converts legend widths from inches to user coordinates.
- `.apply_plot_par()`: applies optional base-graphics styling arguments.
- `.as_plot_path()`: accepts scalar or full-length threshold paths.
- `.thresholds_for_sim_plot()`: chooses simulation threshold paths when
  available, falling back to `sol$x`.
- `.draw_reflected_panel()`: shared state-path panel used by both reflected-path
  figures.

User scripts should normally call the public plotting functions instead of
these helpers.

## Minimal Workflow

```r
source(file.path("src", "load.R"))

p <- make_params(
  b = 0, delta = 1, r = 1, eps = 0.5,
  sigma = 1, mu = 1, u = 1, l = 1
)

sol <- solve_optimal_barriers(
  p, xL0 = -0.5, xk0 = -0.1, xl0 = 0.5, xU0 = 1,
  switch_to_newton = FALSE
)

sim <- simulate_reflected_jd(params = p, thresholds = sol$x, seed = 123)

plot_H(sol, p)
plot_H_prime(sol, p)
plot_reflected_jd(sol, p, sim = sim)
plot_controls(sim)
plot_reflected_with_controls(sol, p, sim = sim)
```
