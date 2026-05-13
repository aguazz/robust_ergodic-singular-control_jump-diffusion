# Animations Reference

This note explains the frame-rendering helpers under `src/animations/`. They are
for producing PNG sequences used in slides, not for building videos directly.

The modules are loaded through `src/load.R` in this order:

1. `src/animations/rendering.R`
2. `src/animations/solver_convergence.R`
3. `src/animations/control_canvas.R`

## Notation

- `X`: reflected state path.
- `U`: cumulative upward reflection at the lower barrier `xL`.
- `L`: cumulative downward reflection at the upper barrier `xU`. Animation
  labels often write this as `D_t`.
- `J_t`: finite-horizon running-average cost.
- `xL`, `xU`: lower and upper reflecting barriers.
- `xk`, `xl`: drift- and intensity-ambiguity thresholds.
- `gamma`: ergodic value returned by the solver.

## Shared Helpers

### `render_and_save_png(fname_base, plotfun, ...)`

Saves one PNG frame. It is intentionally PNG-only so animation frames keep a
stable size and do not incur the PDF-writing overhead used by static figures.

Common arguments:

- `fname_base`: filename without extension.
- `plotfun`: zero-argument plotting function.
- `dir`: output folder.
- `width_in`, `height_in`, `dpi`: frame size and resolution.
- `show`, `save`: draw interactively and/or save to disk.

## Solver-Convergence Frames

### `save_animation_frames_solver_convergence(sol_opt, params, ...)`

Exports one frame per recorded nonlinear-solver iterate. The input solution must
have been created with `record_iterates = TRUE`.

Important arguments:

- `every`: keep every `every`-th iterate.
- `show_axes`, `bty`, `box_lwd`: panel visibility and box styling.
- `tick_cex`, `h_tick_cex`, `axis_title_cex`: axis text controls.
- `text_cex`, `legend_cex`, `point_cex`: controls for the diagnostics column,
  legend, and current-iterate marker.
- `frame_pause`: pause between frames when `show = TRUE`.

The default text sizes reproduce the old animation styling.

### `plot_frame_solver_convergence(history, params, i, lims, ...)`

Draws a single solver-convergence frame. This is mainly useful for testing one
frame before exporting a full sequence.

`lims` should usually come from `solver_anim_limits(history, params)`.

## Control Canvas Frames

### `save_animation_frames_canvas(sol, params, sim, ...)`

Exports the three-left plus one-right Beamer canvas:

- top-left: downward intervention process `D_t` (`sim$L`);
- middle-left: reflected state path `X_t`;
- bottom-left: upward intervention process `U_t`;
- right: running-average cost `J_t`.

Important arguments:

- `show_ambiguity`, `show_barriers`, `show_gamma`: toggle threshold/barrier and
  ergodic-value overlays.
- `show_axis_labels`, `show_x_axes`, `show_y_axes`: global axis controls.
- `bottom_x_axis`, `cost_x_axis`: preserve the old behavior where the bottom
  control panel and cost panel can show their x-axes even when upper panels do
  not.
- `widths`, `heights`, `mar_D`, `mar_X`, `mar_U`, `mar_J`: layout controls.
- `tick_cex`, `axis_title_cex`, `cost_scale_cex`, `legend_cex`,
  `legend_min_cex`, `legend_ncol`, `marker_cex`: text, marker-size, and legend
  controls. The canvas legend keeps the legacy one-column labels and line
  swatches, with a single default scale for both robust and non-robust frames.

Defaults preserve the visual appearance used by `run/render_animations.R`.

### `plot_frame_canvas_3left_1right(sol, params, sim, i, ...)`

Draws one canvas frame. Use this for quick layout checks before producing a full
sequence.

## Reproduction Script

`run/render_animations.R` currently builds:

- `frames/erg_sing_control/`
- `frames/robust_erg_sing_control/`
- `frames/solver_convergence/`

The script does not require `magick`; PNG optimization support was removed
because it was not part of the active rendering workflow.
