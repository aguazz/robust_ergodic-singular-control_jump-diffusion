# Solver Reference

This note explains the main functions in `src/solver.R`: what they are for,
which arguments matter, and what they return. It is intentionally lighter than
R package documentation, but follows the same spirit.

For a first pass through the code, read the functions in this order:

1. `make_params()`
2. `solve_optimal_barriers()`
3. `diagnose()`
4. `build_suboptimal_H()`
5. `IH_at_xlambda_stable()`

The rest of the functions are either fixed-barrier variants or internal
numerical helpers.

## Notation

The code uses short names for the four free-boundary points:

- `xL`: lower reflecting barrier.
- `xk`: drift-ambiguity threshold, written as `x^kappa` in the paper.
- `xl`: intensity-ambiguity threshold, written as `x^lambda` in the paper.
- `xU`: upper reflecting barrier.

The model has two regimes:

- **Regime 1**: `xl <= xU`. All three ambiguity regions appear inside the
  continuation band.
- **Regime 2**: `xl > xU`. The intensity threshold lies above the upper
  reflecting barrier, so Region 3 is outside the band.

## Public Functions

### `make_params(b, delta, r, eps, sigma, mu, u, l)`

Creates the parameter list expected by all solver functions.

Arguments:

- `b`: baseline drift.
- `delta`: radius of drift ambiguity.
- `r`: baseline jump intensity.
- `eps`: relative radius of jump-intensity ambiguity.
- `sigma`: diffusion volatility. Must be positive.
- `mu`: exponential jump-size parameter. Must be positive.
- `u`: upward intervention cost.
- `l`: downward intervention cost.

Returns:

A list containing the inputs plus `EY = -1 / mu`, the mean jump size under the
exponential downward-jump convention used by the model.

Example:

```r
p <- make_params(
  b = 1, delta = 0.1, r = 1, eps = 0.2,
  sigma = 1, mu = 1, u = 1, l = 1
)
```

### `solve_optimal_barriers(p, xL0, xk0, xl0, xU0, ...)`

Solves the full free-boundary problem: lower and upper reflecting barriers plus
the two ambiguity thresholds.

Core arguments:

- `p`: parameter list from `make_params()`.
- `xL0`, `xk0`, `xl0`, `xU0`: initial guesses. They must satisfy
  `xL0 < xk0`, `xk0 < xl0`, and `xk0 < xU0`. `xl0` may be above `xU0`.
- `method`: nonlinear-system method passed to `nleqslv`; defaults to
  `"Broyden"`.
- `control`: control list for the first `nleqslv` stage.
- `tol_build`: tolerance passed to the inner candidate builder for small-gap
  warnings.
- `tol_regime_switch`: tolerance used to treat `xl` very close to `xU` as
  Regime 2, avoiding a nearly degenerate Region 3.
- `switch_to_newton`: if `TRUE`, optionally refines a good Broyden solution
  with a finite-difference Newton step.
- `record_iterates`: if `TRUE`, stores an iteration history useful for
  convergence plots and animations.

Returns:

A list with:

- `H`: function evaluating the solved marginal value derivative.
- `Hp`: function evaluating the derivative of `H`.
- `gamma`: ergodic value.
- `regime2`: logical flag for the solved regime.
- `pieces`: region-by-region coefficients, roots, anchors, and matching
  matrices.
- `params`: the parameter list used in the solve.
- `x`: list with `xL`, `xk`, `xl`, and `xU`.
- `nleqslv_primary`: raw result from the first nonlinear solve.
- `nleqslv_refined`: raw result from the optional Newton refinement, or `NULL`.
- `chosen_solver`: `"broyden"` or `"newton"`.
- `iter_history`: data frame of recorded iterates, or `NULL`.

Example:

```r
sol <- solve_optimal_barriers(
  p, xL0 = -0.5, xk0 = -0.1, xl0 = 0.5, xU0 = 1,
  switch_to_newton = FALSE,
  tol_regime_switch = 1e-3
)

unlist(sol$x)
sol$gamma
```

### `diagnose(sol, tol = 1e-8)`

Checks whether a candidate or solved object satisfies the main matching,
ambiguity, and smooth-fit conditions.

Arguments:

- `sol`: object returned by `build_suboptimal_H()`,
  `solve_optimal_barriers()`, or
  `solve_ambiguity_thresholds_fixed_barriers()`.
- `tol`: pass/fail tolerance.

Returns:

A list with:

- `checks`: data frame with condition names, residual values, type labels, and
  pass flags.
- `condition_numbers`: condition-number estimates for the region matching
  matrices.

Notes:

In Regime 2, continuity and derivative checks at `xl` are marked `NA` because
`xl` is outside the continuation band. The relevant intensity condition is still
reported as `IH(xl)=0`.

### `build_suboptimal_H(p, xL, xk, xl, xU, ...)`

Builds the candidate function `H` for fixed barriers and ambiguity thresholds.
This is the inner solver used repeatedly by the outer nonlinear solve.

Arguments:

- `p`: parameter list from `make_params()`.
- `xL`, `xk`, `xl`, `xU`: fixed barriers and ambiguity thresholds.
- `tol_gap`: warns when a region is very narrow.
- `tol_a1`: tolerance for detecting the degenerate polynomial branch that is
  not implemented.
- `tol_lambda`: tolerance for detecting nearly repeated characteristic roots.
- `kappa3_scale_thresh`: condition-number threshold above which row scaling is
  applied to the small matching systems.
- `tol_regime_switch`: near-boundary tolerance for Regime 2 classification.

Returns:

A list with the same core fields as the solved object:

- `H`, `Hp`, `gamma`, `regime2`, `pieces`, `params`, and `x`.

Use this function when you want to inspect the candidate associated with a
specific set of boundaries, not when you want to solve for the optimal ones.

### `IH_at_xlambda_stable(sol, xla = NULL)`

Evaluates the intensity-ambiguity residual `IH(xla)`. If `xla` is omitted, it
uses `sol$x$xl`.

Arguments:

- `sol`: object returned by `build_suboptimal_H()` or a solver.
- `xla`: optional point at which to evaluate the residual.

Returns:

A single numeric value. At a correctly solved intensity threshold, this should
be close to zero.

Notes:

The implementation uses the anchored exponential representation to avoid
unstable exponentials. In Regime 2, it also includes the constant tail where
`H = l` between `xU` and `xl`.

### `solve_ambiguity_thresholds_fixed_barriers(p, xL, xU, xk0, xl0, ...)`

Solves only for the ambiguity thresholds `xk` and `xl`, holding the reflecting
barriers `xL` and `xU` fixed.

Core arguments:

- `p`: parameter list from `make_params()`.
- `xL`, `xU`: fixed reflecting barriers.
- `xk0`, `xl0`: initial threshold guesses. They must satisfy
  `xL < xk0 < xU` and `xk0 < xl0`.
- Other numerical arguments mirror `solve_optimal_barriers()`.

Returns:

A list with `H`, `Hp`, `gamma`, `regime2`, `pieces`, `params`, `x`,
`nleqslv_primary`, `nleqslv_refined`, and `chosen_solver`.

Use this function for misspecification or comparative-statics exercises where
the controller's reflecting band is held fixed.

## Return Object Anatomy

Most solver objects have the same shape:

```r
list(
  H = function(x) ...,
  Hp = function(x) ...,
  gamma = ...,
  regime2 = TRUE/FALSE,
  pieces = list(region1 = ..., region2 = ..., region3 = ...),
  params = p,
  x = list(xL = ..., xk = ..., xl = ..., xU = ...)
)
```

The `pieces` entries expose the region-level coefficients:

- `u_minus`, `u_plus`: homogeneous coefficients.
- `a`: ODE coefficients and worst-case ambiguity parameters.
- `lam`: characteristic roots.
- `q`: polynomial coefficients.
- `anchor_minus`, `anchor_plus`: per-root anchors used for stable exponentials.
- `M`: matching matrix used to solve the region coefficients.

These fields are mainly for diagnostics and plotting. User scripts usually need
only `sol$x`, `sol$gamma`, `sol$H`, and `sol$Hp`.

## Parameter Transforms

The nonlinear solvers work in unconstrained variables to preserve ordering.
These helpers are useful when debugging residuals:

- `z_from_x(xL, xk, xl, xU)` and `x_from_z(z)` transform the full four-point
  free-boundary problem.
- `z_from_x_fixed_barriers(xL, xk, xl, xU)` and
  `x_from_z_fixed_barriers(z, xL, xU)` transform the fixed-barrier problem.

## Residual Functions

These functions are mostly internal, but useful for debugging:

- `opt_conditions_residuals(z, p, ...)`: four residuals for the full
  free-boundary problem.
- `ambiguity_threshold_residuals_fixed_barriers(z, p, xL, xU, ...)`: two
  residuals for fixed barriers.

In Regime 1, the last threshold residual is a smooth-fit condition at `xl`.
In Regime 2, it becomes `IH(xl)=0`.

## Internal Helpers

The remaining functions are implementation details. The most important ones are:

- `region_ambiguity()`: worst-case drift and intensity distortions by region.
- `region_a_coeffs()`, `lambda_pm()`, and `q_poly_coeffs()`: ODE ingredients.
- `anchors_for_region()`, `g_eval()`, and `gp_eval()`: anchored exponential
  basis functions.
- `solve_scaled_2x2()`: row-scaled 2x2 solves used in the matching systems.
- `anchored_exp_integral()` and `Q_poly()`: homogeneous and polynomial pieces
  of the jump-convolution residual.
- `region_value()` and `region_slope()`: evaluate a solved region and its
  derivative.
- `fd_jacobian()`: central finite-difference Jacobian for optional Newton
  refinement.

These helpers are kept in `solver.R` rather than separate files because they
are small and tightly tied to the solver formulas.
