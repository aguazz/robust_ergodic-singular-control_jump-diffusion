# Simulation Reference

This note explains the simulation helpers in `src/simulation.R`: how to simulate
reflected jump-diffusion paths, how threshold inputs are interpreted, and how to
compute finite-horizon running-average costs.

For normal use, read the functions in this order:

1. `simulate_reflected_jd()`
2. `finite_horizon_ergodic_cost()`

The only internal helper is `.threshold_path()`, which standardizes constant,
vector-valued, and time-dependent threshold inputs.

## Notation

The simulation code uses the same boundary notation as the solver:

- `X`: reflected state path.
- `U`: cumulative upward reflection at the lower barrier `xL`.
- `L`: cumulative downward reflection at the upper barrier `xU`.
- `xL`, `xU`: lower and upper reflecting barriers.
- `xk`: drift-ambiguity threshold.
- `xl`: intensity-ambiguity threshold.
- `dJ`: uncompensated compound-Poisson jump increment.

Some figures label the downward reflection process as `D_t`. The simulation
object keeps the field name `L` for backward compatibility with the rest of the
repo.

## Public Functions

### `simulate_reflected_jd(T = 8, dt = 0.001, x0 = NULL, params, thresholds, seed = NULL)`

Simulates one reflected jump-diffusion path under bang-bang worst-case ambiguity.

Arguments:

- `T`: final simulation time. Must be positive.
- `dt`: time step. Must be positive. The number of steps is `round(T / dt)`.
- `x0`: initial state. If `NULL`, the simulation starts at the midpoint between
  `xL` and `xU` at time zero.
- `params`: parameter list, usually created with `make_params()`.
- `thresholds`: list with `xL`, `xU`, `xk`, and `xl`.
- `seed`: optional random seed for reproducible paths.

Threshold inputs:

Each entry of `thresholds` may be one of:

- a scalar, used as a constant path;
- a vector of length `round(T / dt) + 1`;
- a function of time.

Time-dependent thresholds may be vectorized, for example:

```r
thresholds <- list(
  xL = function(t) sol$x$xL + 0.05 * sin(t / 5),
  xU = function(t) sol$x$xU + 0.10 * cos(t / 5),
  xk = sol$x$xk,
  xl = sol$x$xl
)
```

If a threshold function is not vectorized, the helper falls back to evaluating
it one time point at a time.

Returns:

A list with:

- `time`: simulation grid.
- `X`: reflected state path.
- `U`: cumulative upward reflection at `xL`.
- `L`: cumulative downward reflection at `xU`.
- `jumped`: logical vector indicating whether at least one jump occurred during
  each step.
- `dJ`: compound-Poisson jump increments.
- `xL`, `xU`: lower and upper barrier paths.
- `xk`, `xl`: ambiguity-threshold paths.

Example:

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

sim <- simulate_reflected_jd(
  T = 8, dt = 0.001,
  params = p,
  thresholds = sol$x,
  seed = 123
)
```

Implementation notes:

- The drift ambiguity is selected from the current state relative to `xk`.
- The jump-intensity ambiguity is selected from the current state relative to
  `xl`.
- Jumps are downward exponential jumps.
- Reflection is applied after the Euler/jump update, using the next-time-step
  barrier values.
- The jump increment `dJ` is uncompensated; the drift includes the compensation
  term for the reference intensity `params$r`.

### `finite_horizon_ergodic_cost(sim, params, running = function(x) x^2)`

Computes the pathwise running-average cost

```text
(running-cost integral + u * U_t + l * L_t) / t
```

for one simulated path.

Arguments:

- `sim`: object returned by `simulate_reflected_jd()`. It must contain `time`,
  `X`, `U`, and `L`.
- `params`: parameter list with costs `u` and `l`.
- `running`: running-cost function applied to the simulated state path. The
  default is quadratic cost, `x^2`.

Returns:

A data frame with:

- `time`: simulation grid.
- `cum`: cumulative cost up to each time.
- `avg`: running-average cost. The first value is `NA` to avoid division by
  zero at `t = 0`.

Example:

```r
cost_df <- finite_horizon_ergodic_cost(sim, p)

tail(cost_df)
```

## Internal Helper

### `.threshold_path(value, tt, name)`

Standardizes one threshold input into a numeric vector on the simulation grid
`tt`.

Accepted inputs are scalar constants, full-length vectors, or functions of
time. This helper is internal because user scripts should pass thresholds
through `simulate_reflected_jd()` rather than call it directly.

## Common Usage Patterns

Use solver output directly:

```r
sim <- simulate_reflected_jd(params = p, thresholds = sol$x, seed = 123)
```

Use a fixed reflecting band with moving ambiguity thresholds:

```r
sim <- simulate_reflected_jd(
  params = p,
  thresholds = list(
    xL = sol$x$xL,
    xU = sol$x$xU,
    xk = function(t) sol$x$xk + 0.05 * sin(t),
    xl = function(t) sol$x$xl + 0.05 * cos(t)
  )
)
```

Compute the cost along an already simulated path:

```r
cost_df <- finite_horizon_ergodic_cost(sim, p)
```
