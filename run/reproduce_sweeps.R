source(file.path("src", "load.R"))

# ---------------------- EXAMPLE SWEEPER ---------------------------------------
cost_matrix <- matrix(c(1, 1,
                        2, 1,
                        1, 2), nrow = 3, byrow = T)

## 1) Sweep b
sweep_b <- comparative_sweeper(
  sweep_param  = "b",
  sweep_values = seq(-10, 10, by = 0.01),
  delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1, l = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_b,
  title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
  margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
  axis_title_cex = 1.4,
  plot_gamma = TRUE, gamma_layout = "stacked",
  save = TRUE 
)

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 1) Sweep b
  sweep_b <- comparative_sweeper(
    sweep_param  = "b",
    sweep_values = seq(-20, 20, by = 0.01),
    delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_b,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE  
  )
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 2) Sweep delta
  sweep_delta <- comparative_sweeper(
    sweep_param  = "delta",
    sweep_values = seq(0, 100, by = 0.01),
    b = 0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_delta,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE 
  )
  # Fitting from family of functions
  res <- sweep_delta$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 3) Sweep r
  sweep_r <- comparative_sweeper(
    sweep_param  = "r",
    sweep_values = seq(0.05, 100, by = 0.05),
    b = 0, delta = 1.0, eps = 0.5, sigma = 1, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_r,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE 
  )
  # Fitting from family of functions
  res <- sweep_r$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 4) Sweep eps (epsilon)
  sweep_eps <- comparative_sweeper(
    sweep_param  = "eps",
    sweep_values = seq(0, 1, by = 0.001),
    b = 0, delta = 1.0, r = 1, sigma = 1, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_eps,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE 
  )
  # Fitting from family of functions
  res <- sweep_eps$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 5) Sweep sigma
  sweep_sigma <- comparative_sweeper(
    sweep_param  = "sigma",
    sweep_values = seq(0.1, 50, by = 0.01),
    b = 0, delta = 1.0, r = 1, eps = 0.5, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_sigma,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE 
  )
  # Fitting from family of functions
  res <- sweep_sigma$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 6) Sweep mu
  sweep_mu <- comparative_sweeper(
    sweep_param  = "mu",
    sweep_values = 1/seq(0.51, 5, by = 0.01),
    b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_mu,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE 
  )
  # Fitting from family of functions
  res <- sweep_mu$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 6.1) Sweep 1/mu (E[Y])
  sweep_inv_mu <- comparative_sweeper(
    sweep_param  = "inv_mu",
    sweep_values = seq(0, 3, by = 0.01),
    b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_inv_mu,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE 
  )
  # Fitting from family of functions
  res <- sweep_inv_mu$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

## 7) Sweep u
# b = 2
sweep_u <- comparative_sweeper(
  sweep_param  = "u",
  sweep_values = seq(0.01, 10, by = 0.01),
  b = 2, delta = 1, r = 1, eps = 0.5, sigma = 1, mu = 1, l = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_u,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE 
)
# b = -2
sweep_u <- comparative_sweeper(
  sweep_param  = "u",
  sweep_values = seq(0.01, 10, by = 0.01),
  b = -2, delta = 1, r = 1, eps = 0.5, sigma = 1, mu = 1, l = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_u,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE  
)
# Fitting from family of functions
res <- sweep_u$results
ok  <- res$converged
x <- res$sweep_value[ok]
y <- res$gamma[ok]
fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
fit_obj$best_family
coef(fit_obj$best_fit)
sum(residuals(fit_obj$best_fit)^2)
# Make a smooth curve over the same range
yy <- fit_obj$best_fit$m$fitted()
plot(x, y, pch = 16)
lines(x, yy, lwd = 2, col = "red")

## 8) Sweep l
# b = 2
sweep_l <- comparative_sweeper(
  sweep_param  = "l",
  sweep_values = seq(1, 10, by = 0.01),
  b = 2, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_l,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE 
)
# b = -2
sweep_l <- comparative_sweeper(
  sweep_param  = "l",
  sweep_values = seq(1, 10, by = 0.01),
  b = -2, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_l,
    title = FALSE, show_x_axis_title = TRUE, show_y_axis_title = FALSE,
    margins = c(3.7, 2.2, 1.5, 1.2), axis_mgp = c(2.8, 1, 0), 
    axis_title_cex = 1.4,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE
)
# Fitting from family of functions
res <- sweep_l$results
ok  <- res$converged
x <- res$sweep_value[ok]
y <- res$gamma[ok]
fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
fit_obj$best_family
coef(fit_obj$best_fit)
sum(residuals(fit_obj$best_fit)^2)
# Make a smooth curve over the same range
yy <- fit_obj$best_fit$m$fitted()
plot(x, y, pch = 16)
lines(x, yy, lwd = 2, col = "red")

