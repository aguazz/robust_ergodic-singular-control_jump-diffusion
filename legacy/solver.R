# Core numerical solver for the robust ergodic singular-control problem.
# Source split from legacy/functions_solver_and_images.R; no examples are executed here.

# ================================== STABLE CORE ==================================
# ============================== Sections 5.1 -- 5.3 ==============================
# Implements the anchored inner solver (u_i^±) and optimal solver with robust numerics.
# Differences vs. previous template:
#   • PER-ROOT anchoring (x_{a,i,−}, x_{a,i,+}) from §5.1: anchor left if λ≥0, right if λ<0.
#   • Region 2: IH(x_λ)=0 uses the anchored + stable form (Eq. (49)) via φ_1mexp.
#   • All H, H' constructions and optimality residuals updated to use per-root anchors.
# References: eqs. (25)–(31), (32), (35), (37)–(41) and their STABLE forms (46)–(51),
#             plus optimality conditions (42)–(45) and stable analogues (52)–(55).
# Paper: "Ergodic singular control for ambiguous compound-Poisson jump diffusion processes"
#        (Sections 4–5; esp. 5.1–5.3).  [Anchoring, stability, and γ updates]
# ---------------------------------------------------------------------------------

# --------------------------------- Utilities -------------------------------------
nz <- function(x, tol=1e-16) abs(x) > tol

# Stable "relative exponential": exprel(z) = (exp(z) - 1)/z with series near z=0
exprel <- function(z) {
  out <- z
  small <- abs(z) < 1e-6
  out[!small] <- expm1(z[!small]) / z[!small]
  # 1 + z/2 + z^2/6  when |z| is tiny (guarding cancellation)
  out[ small] <- 1 + 0.5*z[ small] + (z[ small]^2)/6
  out
}

# Stable φ(a, d) = (1 - exp(-a d)) / a  with series near a=0  (used in Eq. (49))
phi_1mexp <- function(a, d) {
  ad <- a * d
  d * exprel(-ad)
}

# Stable safe division
safe_div <- function(num, den, tol=1e-12, msg="Division by ~0") {
  if (any(abs(den) < tol)) stop(msg)
  num / den
}

# Stable Horner polynomial evaluation for quadratic: c2 x^2 + c1 x + d0
poly2_eval  <- function(c2, c1, d0, x) (c2 * x + c1) * x + d0
poly2_prime <- function(c2, c1, x) 2*c2 * x + c1

# ---------------------------- Parameters & regions -------------------------------
make_params <- function(b, delta, r, eps, sigma, mu, u, l) {
  stopifnot(mu > 0, sigma > 0, r > 0, eps >= 0, eps <= 1, delta >= 0, u >= 0, l >= 0)
  EY <- -1/mu
  list(b=b, delta=delta, r=r, eps=eps, sigma=sigma, mu=mu, EY=EY, u=u, l=l)
}

# Worst-case ambiguity per region (Sec. 4.2, eq. (21))
region_ambiguity <- function(i, p) {
  if (i == 1)      list(kappa=-p$delta, lambda=p$r*(1+p$eps))
  else if (i == 2) list(kappa=+p$delta, lambda=p$r*(1+p$eps))
  else if (i == 3) list(kappa=+p$delta, lambda=p$r*(1-p$eps))
  else stop("Region must be 1,2,3.")
}

# a_{i,·} (Eq. (26))
region_a_coeffs <- function(i, p) {
  rc <- region_ambiguity(i, p)
  a_st <- p$b + p$sigma*rc$kappa + p$r/p$mu      # a* = b + σκ + r/μ
  list(
    a1 = p$mu*a_st - rc$lambda,                   # a_{i,1}
    a2 = a_st + 0.5*p$mu*p$sigma^2,               # a_{i,2}
    a3 = 0.5*p$sigma^2,                           # a_{i,3} = σ^2/2
    a_star = a_st,
    lambda_wc = rc$lambda, kappa_wc = rc$kappa
  )
}

# λ_i^± (28) with robust checks (S2)
lambda_pm <- function(a_coeffs, tol_lambda=1e-12) {
  a1 <- a_coeffs$a1; a2 <- a_coeffs$a2; a3 <- a_coeffs$a3
  disc <- a2^2 - 4*a3*a1
  if (disc <= 0) stop("Discriminant ≤ 0: parameters yield no real roots (check σ>0, λ*>0).")
  root <- sqrt(disc)
  lam_plus  <- (a2 + root)/(2*a3)
  lam_minus <- (a2 - root)/(2*a3)
  if (abs(lam_plus - lam_minus) <= tol_lambda) {
    # (S2) repeated/near-coincident roots: not implemented here
    stop(paste0("S2: near-coincident roots in λ± (|Δ| small). ",
                "Use repeated-root branch or relax tol_lambda; ",
                "current |λ+−λ−|=", signif(abs(lam_plus-lam_minus), 6)))
  }
  list(lam_plus=lam_plus, lam_minus=lam_minus, disc=disc)
}

# ------------------------ q_i(x) coefficients (Sec. 4.1.1) ----------------------
# Assumption 2: a_{i,1} != 0 -> quadratic polynomial; otherwise (S1) abort with message
q_poly_coeffs <- function(a_coeffs, p, tol_a1=1e-12) {
  a1 <- a_coeffs$a1; a2 <- a_coeffs$a2; a3 <- a_coeffs$a3; mu <- p$mu
  if (abs(a1) <= tol_a1) {
    # (S1) cubic fallback (Remark 1) not implemented in anchored γ-free decomposition
    stop("S1: |a_{i,1}|≈0 detected. Use the cubic (Remark 1) branch for p_i(x).")
  }
  # (29) with γ-free shift qi = pi - (μ/a1)γ   (qi independent of γ)
  c2 <- - mu / a1
  c1 <-  2*(mu*a2 - a1)/a1^2
  d0 <- ( 2*a1*a2 + 2*mu*a1*a3 - 2*mu*a2^2 )/a1^3
  list(c2=c2, c1=c1, d0=d0, gamma_fac = mu/a1)
}
q_eval      <- function(qc, x) poly2_eval(qc$c2, qc$c1, qc$d0, x)
qprime_eval <- function(qc, x) poly2_prime(qc$c2, qc$c1, x)

# ----------------------------- Anchoring (Sec. 5.1) -----------------------------
# PER-ROOT anchors: for each root λ_i^±, choose left endpoint if λ≥0, right endpoint if λ<0.
anchors_for_region <- function(lam_minus, lam_plus, x_left, x_right) {
  list(
    minus = if (lam_minus >= 0) x_left else x_right,
    plus  = if (lam_plus  >= 0) x_left else x_right
  )
}

# Helpers: anchored exponentials g(x) = exp(-λ (x - x_a)), g'(x) = -λ * g(x)
g_eval  <- function(lambda, x, xa)  exp(-lambda * (x - xa))
gp_eval <- function(lambda, x, xa) -lambda * exp(-lambda * (x - xa))

# ----------------------------- γ via (46) (stable) ------------------------------
# γ = 0.5 σ^2 H1'(x) + A   with A from §4.2.1; derivative uses PER-ROOT anchors (Eq. (46))
gamma_from_xL_stable <- function(xL, xa1m, xa1p, u1m, u1p, lam1, q1c, p) {
  sigma <- p$sigma; mu <- p$mu; r <- p$r; u <- p$u; delta <- p$delta; eps <- p$eps; b <- p$b
  A <- -u*(b - delta*p$sigma + r/mu) + u*(1+eps)*r/mu + xL^2
  H1p <- gp_eval(lam1$lam_minus, xL, xa1m) * u1m +
    gp_eval(lam1$lam_plus , xL, xa1p) * u1p +
    qprime_eval(q1c, xL)
  0.5*sigma^2 * H1p + A
}

# ---------------------------- Matrix solve (2x2, LU) ----------------------------
solve_2x2_LU <- function(M, b, tag="") {
  # Condition estimates (W2)
  k2 <- tryCatch(kappa(M), error=function(e) Inf)
  if (is.finite(k2) && k2 >= 1e8 && k2 <= 1e12)
    warning(sprintf("W2: cond2(M)≈%.3e for %s (ill-conditioned)", k2, tag))
  if (is.infinite(k2) || k2 > 1e12)
    warning(sprintf("W2 (strong): cond2(M)≈%s for %s", format(k2), tag))
  as.numeric(solve(M, b))
}

# -------------- Region 2 integrals for IH(xλ) (stable, eq. (49)) ---------------
# Qi(y) helper: polynomial contribution used inside (49); depends on q_i and γ via (μ/a_{i,1})γ term
Q_poly <- function(y, qi, ai1, gamma, mu) {
  c2 <- qi$c2; c1 <- qi$c1; d0 <- qi$d0
  term2 <- c2*( y^2/mu - 2*y/(mu^2) + 2/(mu^3) )
  term1 <- c1*( y/mu  - 1/(mu^2) )
  term0 <- ( d0 + (mu/ai1)*gamma )*(1/mu)
  term2 + term1 + term0
}

# ============================== Inner solver (stable) ===========================
# Builds suboptimal H using PER-ROOT anchored u_i^± per Sec. 5.2 (46)–(51)
build_suboptimal_H <- function(p, xL, xk, xl, xU,
                               tol_gap=1e-10, tol_a1=1e-12, tol_lambda=1e-12,
                               kappa3_scale_thresh = 50,
                               tol_regime_switch = 1e-3) {
  # NEW: allow crossing xl ~ xU safely (avoid degenerate Region 3)
  tol_sw <- tol_regime_switch * (1 + abs(xU) + abs(xl))
  regime2 <- (xl > xU)
  if (!regime2 && (xU - xl) <= tol_sw) regime2 <- TRUE
  
  # effective split point for the *piecewise H construction*
  xlam_eff <- if (regime2) xU else xl
  
  
  # Gaps and warnings (S3/W1)
  d1 <- xk - xL
  d2 <- xl - xk
  d3eff <- xU - xlam_eff
  
  if (d1 <= tol_gap) warning("S3/W1: small gap on I1 (xk - xL) → near-collinearity.")
  if (d2 <= tol_gap) warning("S3/W1: small gap on I2 (xl - xk) → near-collinearity.")
  if (!regime2 && d3eff <= tol_gap) warning("S3/W1: small gap on I3 (xU - xl) → near-collinearity.")
  
  
  # Per-region coefficients (a_{i,·}), λ±, and q_i
  a1c <- region_a_coeffs(1, p)
  a2c <- region_a_coeffs(2, p)
  a3c <- region_a_coeffs(3, p)
  
  lam1 <- lambda_pm(a1c, tol_lambda=tol_lambda)
  lam2 <- lambda_pm(a2c, tol_lambda=tol_lambda)
  lam3 <- lambda_pm(a3c, tol_lambda=tol_lambda)
  
  q1c <- q_poly_coeffs(a1c, p, tol_a1)
  q2c <- q_poly_coeffs(a2c, p, tol_a1)
  q3c <- q_poly_coeffs(a3c, p, tol_a1)
  
  # PER-ROOT anchors
  anc1 <- anchors_for_region(lam1$lam_minus, lam1$lam_plus, xL, xk)
  anc2 <- anchors_for_region(lam2$lam_minus, lam2$lam_plus, xk, xlam_eff)
  
  # Region 3 only exists in Regime 1 (xl <= xU)
  if (!regime2) {
    anc3 <- anchors_for_region(lam3$lam_minus, lam3$lam_plus, xl, xU)
  } else {
    anc3 <- list(minus = NA_real_, plus = NA_real_)
  }
  
  # Convenience closures for anchored exponentials
  g1m  <- function(x) g_eval (lam1$lam_minus, x, anc1$minus);  gp1m <- function(x) gp_eval(lam1$lam_minus, x, anc1$minus)
  g1p  <- function(x) g_eval (lam1$lam_plus , x, anc1$plus );  gp1p <- function(x) gp_eval(lam1$lam_plus , x, anc1$plus )
  g2m  <- function(x) g_eval (lam2$lam_minus, x, anc2$minus);  gp2m <- function(x) gp_eval(lam2$lam_minus, x, anc2$minus)
  g2p  <- function(x) g_eval (lam2$lam_plus , x, anc2$plus );  gp2p <- function(x) gp_eval(lam2$lam_plus , x, anc2$plus )
  g3m  <- function(x) g_eval (lam3$lam_minus, x, anc3$minus);  gp3m <- function(x) gp_eval(lam3$lam_minus, x, anc3$minus)
  g3p  <- function(x) g_eval (lam3$lam_plus , x, anc3$plus );  gp3p <- function(x) gp_eval(lam3$lam_plus , x, anc3$plus )
  
  # ------------------------ Region 1: solve u1^± (stable Eq. (47)) ----------------
  u <- p$u; mu <- p$mu; sigma2 <- p$sigma^2; a11 <- a1c$a1
  
  # m rows follow §5.2.2 with derivative term evaluated at xL in both equations
  m11m <- g1m(xL) + (mu/a11)*(sigma2/2)*gp1m(xL)
  m11p <- g1p(xL) + (mu/a11)*(sigma2/2)*gp1p(xL)
  m12m <- g1m(xk) + (mu/a11)*(sigma2/2)*gp1m(xL)
  m12p <- g1p(xk) + (mu/a11)*(sigma2/2)*gp1p(xL)
  
  A   <- -u*(p$b - p$delta*p$sigma + p$r/mu) + u*(1+p$eps)*p$r/mu + xL^2
  b11 <- -u - q_eval(q1c, xL) - (mu/a11)*((sigma2/2)*qprime_eval(q1c, xL) + A)
  b12 <-    0 - q_eval(q1c, xk) - (mu/a11)*((sigma2/2)*qprime_eval(q1c, xL) + A)
  
  M1_raw <- rbind(c(m11m, m11p), c(m12m, m12p))
  b1_raw <- c(b11, b12)
  
  # Condition number before scaling (to decide whether to scale)
  k1_raw <- tryCatch(kappa(M1_raw), error = function(e) Inf)
  
  if (is.finite(k1_raw) && k1_raw > kappa3_scale_thresh) {
    # Simple and effective: row scaling by row norms (or max-abs)
    s1r <- max(1, sqrt(sum(M1_raw[1,]^2)))  # guard with max(1, ·) to avoid over-scaling tiny rows
    s2r <- max(1, sqrt(sum(M1_raw[2,]^2)))
    M1  <- rbind(M1_raw[1,] / s1r, M1_raw[2,] / s2r)
    b1  <- c(b1_raw[1] / s1r, b1_raw[2] / s2r)
    s1  <- solve_2x2_LU(M1, b1, tag = "Region 1 (u1±) [row-scaled]")
  } else {
    M1  <- M1_raw
    b1  <- b1_raw
    s1  <- solve_2x2_LU(M1, b1, tag = "Region 1 (u1±)")
  }
  
  u1m <- s1[1]; u1p <- s1[2]
  
  # γ from (46) using anchored H1′ at xL
  gamma <- gamma_from_xL_stable(xL, anc1$minus, anc1$plus, u1m, u1p, lam1, q1c, p)
  
  # ------------------------ Region 2: solve u2^± ----------------
  if (!regime2) {
    # ===== Regime 1 (xl <= xU): your original IH(xλ)=0 system =====
    
    # Row 1: H-continuity at xκ
    m21m <- g2m(xk)
    m21p <- g2p(xk)
    b21  <- -(mu/a2c$a1)*gamma - q_eval(q2c, xk)
    
    # Row 2: IH(xλ)=0 (anchored, stable form, Eq. (49))
    m22m <- mu * exp(-lam2$lam_minus * (xl - anc2$minus)) * phi_1mexp(mu - lam2$lam_minus, (xl - xk))
    m22p <- mu * exp(-lam2$lam_plus  * (xl - anc2$plus )) * phi_1mexp(mu - lam2$lam_plus , (xl - xk))
    
    # RHS β(2)_2: region 1 homogeneous + polynomials (anchored per root)
    if (abs(mu - lam1$lam_minus) < 10^-6) {
      term_c1_minus <-
        mu * u1m * exp(-mu*(xl - xk)) * exp(-lam1$lam_minus*(xk - anc1$minus)) *
        phi_1mexp(mu - lam1$lam_minus, (xk - xL))
    } else {
      term_c1_minus <-
        mu * u1m * (
          exp(-mu*(xl - xk)) * exp(-lam1$lam_minus*(xk - anc1$minus)) -
            exp(-mu*(xl - xL)) * exp(-lam1$lam_minus*(xL - anc1$minus))
        ) / (mu - lam1$lam_minus)
    }
    
    if (abs(mu - lam1$lam_plus) < 10^-6) {
      term_c1_plus <-
        mu * u1p * exp(-mu*(xl - xk)) * exp(-lam1$lam_plus*(xk - anc1$plus)) *
        phi_1mexp(mu - lam1$lam_plus, (xk - xL))
    } else {
      term_c1_plus <-
        mu * u1p * (
          exp(-mu*(xl - xk)) * exp(-lam1$lam_plus*(xk - anc1$plus)) -
            exp(-mu*(xl - xL)) * exp(-lam1$lam_plus*(xL - anc1$plus))
        ) / (mu - lam1$lam_plus)
    }
    
    term_c1 <- term_c1_minus + term_c1_plus
    
    Q1_xk <- Q_poly(xk, q1c, a1c$a1, gamma, mu)
    Q1_xL <- Q_poly(xL, q1c, a1c$a1, gamma, mu)
    Q2_xl <- Q_poly(xl, q2c, a2c$a1, gamma, mu)
    Q2_xk <- Q_poly(xk, q2c, a2c$a1, gamma, mu)
    
    term_poly <- mu*exp(-mu*xl) * (exp(mu*xk)*Q1_xk - exp(mu*xL)*Q1_xL + exp(mu*xl)*Q2_xl - exp(mu*xk)*Q2_xk)
    b22 <- + u*exp(mu*(xL - xl)) - term_c1 - term_poly
    
    M2_raw <- rbind(c(m21m, m21p), c(m22m, m22p))
    b2_raw <- c(b21, b22)
    
    k2_raw <- tryCatch(kappa(M2_raw), error = function(e) Inf)
    if (is.finite(k2_raw) && k2_raw > kappa3_scale_thresh) {
      s1r <- max(1, sqrt(sum(M2_raw[1,]^2)))
      s2r <- max(1, sqrt(sum(M2_raw[2,]^2)))
      M2  <- rbind(M2_raw[1,] / s1r, M2_raw[2,] / s2r)
      b2  <- c(b2_raw[1] / s1r, b2_raw[2] / s2r)
      s2  <- solve_2x2_LU(M2, b2, tag = "Region 2 (u2±) [row-scaled]")
    } else {
      M2 <- M2_raw
      b2 <- b2_raw
      s2 <- solve_2x2_LU(M2, b2, tag = "Region 2 (u2±)")
    }
    
    u2m <- s2[1]; u2p <- s2[2]
    
  } else {
    # ===== Regime 2 (xl > xU): NO region 3 inside the band. Solve u2± by H2(xk)=0 and H2(xU)=l =====
    
    # Row 1: H2(xk)=0
    m21m <- g2m(xk)
    m21p <- g2p(xk)
    b21  <- -(mu/a2c$a1)*gamma - q_eval(q2c, xk)
    
    # Row 2: H2(xU)=l  (since H is constant = l for x>=xU)
    m22m <- g2m(xU)
    m22p <- g2p(xU)
    b22  <- p$l - (mu/a2c$a1)*gamma - q_eval(q2c, xU)
    
    M2_raw <- rbind(c(m21m, m21p), c(m22m, m22p))
    b2_raw <- c(b21, b22)
    
    k2_raw <- tryCatch(kappa(M2_raw), error = function(e) Inf)
    if (is.finite(k2_raw) && k2_raw > kappa3_scale_thresh) {
      s1r <- max(1, sqrt(sum(M2_raw[1,]^2)))
      s2r <- max(1, sqrt(sum(M2_raw[2,]^2)))
      M2  <- rbind(M2_raw[1,] / s1r, M2_raw[2,] / s2r)
      b2  <- c(b2_raw[1] / s1r, b2_raw[2] / s2r)
      s2  <- solve_2x2_LU(M2, b2, tag = "Region 2 (u2±) [Regime 2, row-scaled]")
    } else {
      M2 <- M2_raw
      b2 <- b2_raw
      s2 <- solve_2x2_LU(M2, b2, tag = "Region 2 (u2±) [Regime 2]")
    }
    
    u2m <- s2[1]; u2p <- s2[2]
  }
  
  
  # ------------------------ Region 3: solve u3^± -------------------
  if (!regime2) {
    # (your original Region 3 block stays EXACTLY the same)
    
    # Row 1: H3(xU) = l
    m31m <- g3m(xU)
    m31p <- g3p(xU)
    b31  <- p$l - (mu/a3c$a1)*gamma - q_eval(q3c, xU)
    
    # Row 2: continuity at xλ (use anchors xa3 for region 3, xa2 for region 2)
    m32m <- g3m(xl)
    m32p <- g3p(xl)
    b32  <- u2m * g2m(xl) + u2p * g2p(xl) +
      mu*gamma*(1/a2c$a1 - 1/a3c$a1) + q_eval(q2c, xl) - q_eval(q3c, xl)
    
    M3_raw <- rbind(c(m31m, m31p), c(m32m, m32p))
    b3_raw <- c(b31, b32)
    
    # Condition number before scaling (to decide whether to scale)
    k3_raw <- tryCatch(kappa(M3_raw), error = function(e) Inf)
    
    if (is.finite(k3_raw) && k3_raw > kappa3_scale_thresh) {
      # Simple and effective: row scaling by row norms (or max-abs)
      s1 <- max(1, sqrt(sum(M3_raw[1,]^2)))  # guard with max(1, ·) to avoid over-scaling tiny rows
      s2 <- max(1, sqrt(sum(M3_raw[2,]^2)))
      M3  <- rbind(M3_raw[1,] / s1, M3_raw[2,] / s2)
      b3  <- c(b3_raw[1] / s1, b3_raw[2] / s2)
      s3  <- solve_2x2_LU(M3, b3, tag = "Region 3 (u3±) [row-scaled]")
      scaled_row_factors <- c(s1, s2)
    } else {
      M3  <- M3_raw
      b3  <- b3_raw
      s3  <- solve_2x2_LU(M3, b3, tag = "Region 3 (u3±)")
      scaled_row_factors <- c(1, 1)
    }
    
    u3m <- s3[1]; u3p <- s3[2]
  
  } else {
    # Regime 2: no Region 3 inside [xL, xU]. We keep placeholders for compatibility.
    u3m <- 0; u3p <- 0
    M3  <- matrix(NA_real_, 2, 2)
  }
  
  # ----------------------------- Build H and H' -----------------------------------
  H1_fun  <- function(x) u1m*g1m(x) + u1p*g1p(x) + q_eval(q1c, x) + q1c$gamma_fac*gamma
  H2_fun  <- function(x) u2m*g2m(x) + u2p*g2p(x) + q_eval(q2c, x) + q2c$gamma_fac*gamma
  
  if (!regime2) {
    H3_fun  <- function(x) u3m*g3m(x) + u3p*g3p(x) + q_eval(q3c, x) + q3c$gamma_fac*gamma
    Hp3_fun <- function(x) gp3m(x)*u3m + gp3p(x)*u3p + qprime_eval(q3c, x)
  } else {
    H3_fun  <- function(x) rep(p$l, length(x))
    Hp3_fun <- function(x) rep(0, length(x))
  }
  
  Hp1_fun <- function(x) gp1m(x)*u1m + gp1p(x)*u1p + qprime_eval(q1c, x)
  Hp2_fun <- function(x) gp2m(x)*u2m + gp2p(x)*u2p + qprime_eval(q2c, x)
  
  
  H <- function(x) {
    x <- as.numeric(x); out <- numeric(length(x))
    left  <- which(x <= xL)
    mid1  <- which(x >  xL & x <  xk)
    mid2  <- which(x >= xk & x <  xlam_eff)
    mid3  <- which(x >= xlam_eff & x <  xU)
    right <- which(x >= xU)
    if (length(left))  out[left]  <- -p$u
    if (length(mid1))  out[mid1]  <- H1_fun(x[mid1])
    if (length(mid2))  out[mid2]  <- H2_fun(x[mid2])
    if (length(mid3))  out[mid3]  <- H3_fun(x[mid3])
    if (length(right)) out[right] <-  p$l
    out
  }
  
  Hp <- function(x) {
    x <- as.numeric(x); out <- numeric(length(x))
    left  <- which(x <= xL)
    mid1  <- which(x >  xL & x <  xk)
    mid2  <- which(x >= xk & x <  xlam_eff)
    mid3  <- which(x >= xlam_eff & x <  xU)
    right <- which(x >= xU)
    if (length(left))  out[left]  <- 0
    if (length(mid1))  out[mid1]  <- Hp1_fun(x[mid1])
    if (length(mid2))  out[mid2]  <- Hp2_fun(x[mid2])
    if (length(mid3))  out[mid3]  <- Hp3_fun(x[mid3])
    if (length(right)) out[right] <- 0
    out
  }
  
  list(
    H=H, Hp=Hp, gamma=gamma,
    regime2 = regime2,
    pieces=list(
      region1=list(u_minus=u1m, u_plus=u1p, a=a1c, lam=lam1, q=q1c,
                   anchor_minus=anc1$minus, anchor_plus=anc1$plus, M=M1),
      region2=list(u_minus=u2m, u_plus=u2p, a=a2c, lam=lam2, q=q2c,
                   anchor_minus=anc2$minus, anchor_plus=anc2$plus, M=M2),
      region3=list(
        u_minus=u3m, u_plus=u3p, a=a3c, lam=lam3, q=q3c, 
        anchor_minus=anc3$minus, anchor_plus=anc3$plus, M=M3
      )
    ),
    params=p, x=list(xL=xL, xk=xk, xl=xl, xU=xU)
  )
}

# ------------------------ IH(xλ) residual (stable, eq. (49)) --------------------
# ------------------------ IH(xλ) residual (stable, generalized) ----------------
# NEW: works both when xλ<=xU (Regime 1) and when xλ>xU (Regime 2).
# In Regime 2, it correctly adds the constant-tail contribution from [xU, xλ].
IH_at_xlambda_stable <- function(sol, xla = NULL) {
  if (is.null(xla)) xla <- sol$x$xl
  
  p  <- sol$params
  mu <- p$mu
  xs <- sol$x
  xL <- xs$xL; xk <- xs$xk; xU <- xs$xU
  
  # last point that is inside the band (integration upper limit for the interior pieces)
  B <- min(xla, xU)
  
  R1 <- sol$pieces$region1
  R2 <- sol$pieces$region2
  
  q1c <- R1$q; q2c <- R2$q
  a1c <- R1$a; a2c <- R2$a
  lam1 <- R1$lam; lam2 <- R2$lam
  xa1m <- R1$anchor_minus; xa1p <- R1$anchor_plus
  xa2m <- R2$anchor_minus; xa2p <- R2$anchor_plus
  u1m <- R1$u_minus; u1p <- R1$u_plus
  u2m <- R2$u_minus; u2p <- R2$u_plus
  gamma <- sol$gamma
  
  # ----- Region 1 homogeneous contribution over [xL, xk], evaluated at xla -----
  if (abs(mu - lam1$lam_minus) < 1e-6) {
    term_c1_minus <-
      mu * u1m * exp(-mu*(xla - xk)) * exp(-lam1$lam_minus*(xk - xa1m)) *
      phi_1mexp(mu - lam1$lam_minus, (xk - xL))
  } else {
    term_c1_minus <-
      mu * u1m * (
        exp(-mu*(xla - xk)) * exp(-lam1$lam_minus*(xk - xa1m)) -
          exp(-mu*(xla - xL)) * exp(-lam1$lam_minus*(xL - xa1m))
      ) / (mu - lam1$lam_minus)
  }
  
  if (abs(mu - lam1$lam_plus) < 1e-6) {
    term_c1_plus <-
      mu * u1p * exp(-mu*(xla - xk)) * exp(-lam1$lam_plus*(xk - xa1p)) *
      phi_1mexp(mu - lam1$lam_plus, (xk - xL))
  } else {
    term_c1_plus <-
      mu * u1p * (
        exp(-mu*(xla - xk)) * exp(-lam1$lam_plus*(xk - xa1p)) -
          exp(-mu*(xla - xL)) * exp(-lam1$lam_plus*(xL - xa1p))
      ) / (mu - lam1$lam_plus)
  }
  
  term_c1 <- term_c1_minus + term_c1_plus
  
  # ----- Region 2 homogeneous contribution over [xk, B], evaluated at xla -----
  # NOTE: extra factor exp(-mu*(xla-B)) appears if xla>B (i.e. Regime 2)
  shift_B <- exp(-mu * (xla - B))
  
  term_c2 <-
    shift_B * (
      mu * u2m * exp(-lam2$lam_minus*(B - xa2m)) * phi_1mexp(mu - lam2$lam_minus, (B - xk)) +
        mu * u2p * exp(-lam2$lam_plus *(B - xa2p)) * phi_1mexp(mu - lam2$lam_plus , (B - xk))
    )
  
  # ----- Polynomial contributions over [xL, xk] and [xk, B] -----
  Q1_xk <- Q_poly(xk, q1c, a1c$a1, gamma, mu)
  Q1_xL <- Q_poly(xL, q1c, a1c$a1, gamma, mu)
  Q2_B  <- Q_poly(B , q2c, a2c$a1, gamma, mu)
  Q2_xk <- Q_poly(xk, q2c, a2c$a1, gamma, mu)
  
  term_poly <-
    mu * exp(-mu*xla) *
    (exp(mu*xk)*Q1_xk - exp(mu*xL)*Q1_xL + exp(mu*B)*Q2_B - exp(mu*xk)*Q2_xk)
  
  # ----- Constant tail when xla > xU: contribution from [xU, xla] where H=l -----
  tail <- if (xla > xU) p$l * (1 - exp(mu*(xU - xla))) else 0
  
  # ----- Full residual (should be 0 at the true xλ) -----
  -p$u * exp(mu*(xL - xla)) + term_c1 + term_c2 + term_poly + tail
}

# ================================= Diagnostics ==================================
diagnose <- function(sol, tol=1e-8) {
  p  <- sol$params; xs <- sol$x
  xL <- xs$xL; xk <- xs$xk; xl <- xs$xl; xU <- xs$xU
  R1 <- sol$pieces$region1; R2 <- sol$pieces$region2; R3 <- sol$pieces$region3
  gamma <- sol$gamma
  
  # Closures with PER-ROOT anchors
  g  <- function(lambda, x, xa)  exp(-lambda * (x - xa))
  gp <- function(lambda, x, xa) -lambda * exp(-lambda * (x - xa))
  
  H1_fun <- function(x) R1$u_minus*g(R1$lam$lam_minus,x,R1$anchor_minus) + R1$u_plus*g(R1$lam$lam_plus,x,R1$anchor_plus) +
    q_eval(R1$q, x) + R1$q$gamma_fac*gamma
  H2_fun <- function(x) R2$u_minus*g(R2$lam$lam_minus,x,R2$anchor_minus) + R2$u_plus*g(R2$lam$lam_plus,x,R2$anchor_plus) +
    q_eval(R2$q, x) + R2$q$gamma_fac*gamma
  H3_fun <- function(x) R3$u_minus*g(R3$lam$lam_minus,x,R3$anchor_minus) + R3$u_plus*g(R3$lam$lam_plus,x,R3$anchor_plus) +
    q_eval(R3$q, x) + R3$q$gamma_fac*gamma
  
  H1p_at <- function(x) gp(R1$lam$lam_minus,x,R1$anchor_minus)*R1$u_minus + gp(R1$lam$lam_plus,x,R1$anchor_plus)*R1$u_plus + qprime_eval(R1$q, x)
  H2p_at <- function(x) gp(R2$lam$lam_minus,x,R2$anchor_minus)*R2$u_minus + gp(R2$lam$lam_plus,x,R2$anchor_plus)*R2$u_plus + qprime_eval(R2$q, x)
  H3p_at <- function(x) gp(R3$lam$lam_minus,x,R3$anchor_minus)*R3$u_minus + gp(R3$lam$lam_plus,x,R3$anchor_plus)*R3$u_plus + qprime_eval(R3$q, x)
  
  cont_xL <- H1_fun(xL) + p$u
  cont_xk <- H1_fun(xk) - H2_fun(xk)
  cont_xl <- H2_fun(xl) - H3_fun(xl)
  cont_xU <- H3_fun(xU) - p$l
  
  root_k <- H1_fun(xk)               # H(xκ)=0
  IH_xl  <- IH_at_xlambda_stable(sol) # IH(xλ)=0
  
  dcont_xL <- H1p_at(xL) - 0
  dcont_xk <- H2p_at(xk) - H1p_at(xk)
  dcont_xl <- H3p_at(xl) - H2p_at(xl)
  dcont_xU <- 0 - H3p_at(xU)
  
  checks <- data.frame(
    check = c("continuity@xL", "continuity@xk", "continuity@xl", "continuity@xU",
              "H(xk)=0", "IH(xl)=0",
              "H'(jump)@xL", "H'(jump)@xk", "H'(jump)@xl", "H'(jump)@xU"),
    value = c(cont_xL, cont_xk, cont_xl, cont_xU, root_k, IH_xl,
              dcont_xL, dcont_xk, dcont_xl, dcont_xU),
    type  = c(rep("continuity", 4), rep("ambiguity_condition", 2), rep("derivative_jump", 4)),
    pass  = abs(c(cont_xL, cont_xk, cont_xl, cont_xU, root_k, IH_xl,
                  dcont_xL, dcont_xk, dcont_xl, dcont_xU)) <= tol
  )
  rownames(checks) <- NULL
  
  kappa2 <- function(M) {
    tryCatch(kappa(M), error=function(e) Inf)
  }
  kappas <- data.frame(
    block   = c("M1 (Region1)", "M2 (Region2)", "M3 (Region3)"),
    kappa2  = c(kappa2(R1$M), kappa2(R2$M), kappa2(R3$M))
  )
  list(checks=checks, condition_numbers=kappas)
}

# =========================== Optimal barriers (stable) ===========================
# NEW: allow xl and xU to be on either side, as long as both are > xk.
z_from_x <- function(xL, xk, xl, xU) {
  stopifnot(xL < xk, xk < xl, xk < xU)
  c(
    xL,
    log(xk - xL),   # > 0
    log(xl - xk),   # > 0
    log(xU - xk)    # > 0
  )
}

x_from_z <- function(z) {
  stopifnot(length(z) == 4)
  xL <- z[1]
  d1 <- exp(z[2])     # xk - xL
  dL <- exp(z[3])     # xl - xk
  dU <- exp(z[4])     # xU - xk
  xk <- xL + d1
  xl <- xk + dL
  xU <- xk + dU
  c(xL = xL, xk = xk, xl = xl, xU = xU)
}


# Residuals of optimality equations with Regime-2 support (xl may exceed xU)
opt_conditions_residuals <- function(
    z, p, tol_build=1e-10, tol_regime_switch = 1e-3,
    return_sol = FALSE,
    recorder = NULL               # <-- NEW
) {
  xs  <- x_from_z(z)
  sol <- build_suboptimal_H(
    p,
    as.numeric(xs["xL"]), as.numeric(xs["xk"]),
    as.numeric(xs["xl"]), as.numeric(xs["xU"]),
    tol_regime_switch = tol_regime_switch
  )
  
  xL <- xs["xL"]; xk <- xs["xk"]; xl <- xs["xl"]; xU <- xs["xU"]
  R1 <- sol$pieces$region1; R2 <- sol$pieces$region2
  
  H1p <- function(x) gp_eval(R1$lam$lam_minus,x,R1$anchor_minus)*R1$u_minus +
    gp_eval(R1$lam$lam_plus ,x,R1$anchor_plus )*R1$u_plus  +
    qprime_eval(R1$q, x)
  
  H2p <- function(x) gp_eval(R2$lam$lam_minus,x,R2$anchor_minus)*R2$u_minus +
    gp_eval(R2$lam$lam_plus ,x,R2$anchor_plus )*R2$u_plus  +
    qprime_eval(R2$q, x)
  
  if (!isTRUE(sol$regime2)) {
    R3 <- sol$pieces$region3
    H3p <- function(x) gp_eval(R3$lam$lam_minus,x,R3$anchor_minus)*R3$u_minus +
      gp_eval(R3$lam$lam_plus ,x,R3$anchor_plus )*R3$u_plus  +
      qprime_eval(R3$q, x)
    
    r <- c(
      H1p(xL),
      H3p(xU),
      H1p(xk) - H2p(xk),
      H2p(xl) - H3p(xl)
    )
  } else {
    r <- c(
      H1p(xL),
      H2p(xU),
      H1p(xk) - H2p(xk),
      IH_at_xlambda_stable(sol, xla = xl)
    )
  }
  
  # ---------- NEW: record here (side-effect), sol is guaranteed available ----------
  if (is.function(recorder)) {
    # never let recording crash the solver
    try(recorder(z, r, sol), silent = TRUE)
  }
  
  if (isTRUE(return_sol)) attr(r, "sol") <- sol
  r
}

# ------------------------ Outer solver: Broyden → Newton ------------------------
# Switch to a pure Newton method (with FD Jacobian) when close to a solution.
# Nearness is detected by max|f| <= switch_ftol (and after at least switch_iter_min iterations).
# Uses the stable inner builder from Secs. 5.1–5.3.

solve_optimal_barriers <- function(
    p, xL0, xk0, xl0, xU0,
    method = "Broyden",
    control = list(ftol=1e-11, xtol=1e-11, maxit=1e6, stepmax=1, trace=0),
    tol_build = 1e-10,
    verbose = FALSE,
    tol_regime_switch = 1e-3,
    # --- hybrid controls ---
    switch_to_newton = FALSE,
    switch_ftol = 5e-2,          # threshold to trigger Newton (on max|f|)
    switch_iter_min = 2L,        # don't switch on the very first step
    newton_control = list(ftol=1e-12, xtol=1e-12, maxit=1e5, stepmax=1, trace=0),
    fd_step = 1e-6,               # relative step for finite-difference Jacobian
    record_iterates = FALSE,
    record_xtol = 1e-12
) {
  if (!requireNamespace("nleqslv", quietly=TRUE))
    stop("Please install.packages('nleqslv')")
  if (!(xL0 < xk0 && xk0 < xl0 && xk0 < xU0))
    stop("Initial guess must satisfy xL0 < xk0, and xk0 < xl0, xk0 < xU0 (xl0 may be > xU0).")
  
  # ---- reparam z <-> x (keeps barriers ordered) ----
  z0 <- z_from_x(xL0, xk0, xl0, xU0)
  
  # -------------------- iterate recorder (NEW) --------------------
  rec_env <- new.env(parent = emptyenv())
  rec_env$hist <- list()
  rec_env$last_z <- NULL
  rec_env$record_on <- TRUE
  rec_env$stage <- "broyden"
  
  # .c1_jumps_from_sol <- function(sol) {
  #   xs <- sol$x
  #   xL <- xs$xL; xk <- xs$xk; xl <- xs$xl; xU <- xs$xU
  #   
  #   # epsilon chosen relative to smallest gap
  #   gaps <- c(xk - xL, xl - xk, xU - min(xl, xU))
  #   gmin <- max(1e-8, min(gaps[gaps > 0], na.rm=TRUE))
  #   eps  <- min(1e-6 * (1 + max(abs(c(xL,xk,xl,xU)))), 0.25*gmin)
  #   
  #   Hp <- sol$Hp
  #   c(
  #     Hp_xL  = as.numeric(Hp(xL + eps)),
  #     dHp_xk = as.numeric(Hp(xk - eps) - Hp(xk + eps)),
  #     dHp_xl = as.numeric(Hp(xl - eps) - Hp(xl + eps)),
  #     Hp_xU  = as.numeric(Hp(xU - eps))
  #   )
  # }
  
  .record_step <- function(z, r, sol) {
    if (!isTRUE(record_iterates)) return(invisible())
    
    # FORCE a copy (do NOT use as.numeric(z) alone)
    znum <- c(unname(as.double(z)))   # c(...) forces a new vector
    
    if (!is.null(rec_env$last_z)) {
      dz <- suppressWarnings(max(abs(znum - rec_env$last_z)))
      if (is.finite(dz) && dz <= record_xtol) return(invisible())
    }
    
    rec_env$last_z <- znum   # already a fresh copy
    
    # old
    # jumps <- .c1_jumps_from_sol(sol)
    
    # new: store the true residuals that define fmax
    jumps <- c(
      Hp_xL  = as.numeric(r[1]),
      dHp_xk = as.numeric(r[3]),
      dHp_xl = as.numeric(r[4]),  # note: in regime2 this is NOT a derivative jump
      Hp_xU  = as.numeric(r[2])
    )
    
    xs <- sol$x
    xL <- as.numeric(xs$xL); xk <- as.numeric(xs$xk)
    xl <- as.numeric(xs$xl); xU <- as.numeric(xs$xU)
    
    rec_env$hist[[length(rec_env$hist) + 1L]] <- list(
      stage = rec_env$stage,
      z = znum,
      xL = xL, xk = xk, xl = xl, xU = xU,
      gamma = as.numeric(sol$gamma),
      fmax  = max(abs(as.numeric(r))),
      Hp_xL = jumps["Hp_xL"], dHp_xk = jumps["dHp_xk"],
      dHp_xl = jumps["dHp_xl"], Hp_xU = jumps["Hp_xU"]
    )
    invisible()
  }
  
  # Residuals (Eqs. (52)–(55) in stable form)
  .fn_eval <- function(z, do_record = TRUE) {
    
    rec_cb <- if (isTRUE(record_iterates) && isTRUE(do_record)) {
      function(z_loc, r_loc, sol_loc) {
        .record_step(z_loc, as.numeric(r_loc), sol_loc)
      }
    } else NULL
    
    r <- opt_conditions_residuals(
      z, p,
      tol_build = tol_build,
      tol_regime_switch = tol_regime_switch,
      return_sol = FALSE,         # <-- no longer needed
      recorder = rec_cb           # <-- key change
    )
    
    r_num <- as.numeric(r)
    
    if (isTRUE(verbose)) {
      xs <- x_from_z(z)
      cat(sprintf(
        "xL=%.6g  xk=%.6g  xl=%.6g  xU=%.6g  |  r=[% .3e, % .3e, % .3e, % .3e]\n",
        xs["xL"], xs["xk"], xs["xl"], xs["xU"], r_num[1], r_num[2], r_num[3], r_num[4]
      ))
    }
    r_num
  }
  
  fn       <- function(z) .fn_eval(z, do_record = TRUE)
  fn_norec <- function(z) .fn_eval(z, do_record = FALSE)
  
  # Stable, central-difference Jacobian (O(h^2))
  jac_fd <- function(z, fz = NULL) {
    n <- length(z); J <- matrix(0.0, n, n)
    if (is.null(fz)) fz <- fn_norec(z)
    
    for (j in 1:n) {
      h <- fd_step * (1 + abs(z[j]))
      zp <- z; zm <- z
      zp[j] <- z[j] + h;  zm[j] <- z[j] - h
      fp <- fn_norec(zp); fm <- fn_norec(zm)
      J[, j] <- (fp - fm) / (2*h)
    }
    J
  }
  
  # -------------------- Phase 1: robust Broyden --------------------
  ans_broyden <- nleqslv::nleqslv(
    x = z0, fn = fn, method = method,
    global = "dbldog", xscalm = "auto", control = control
  )
  
  if (isTRUE(record_iterates)) {
    # record the final point returned by nleqslv (in case it wasn't evaluated last)
    .fn_eval(ans_broyden$x, do_record = TRUE)
  }
  
  if (isTRUE(record_iterates)) {
    z_end  <- as.numeric(ans_broyden$x)
    xs_end <- x_from_z(z_end)
    sol_end <- build_suboptimal_H(
      p,
      as.numeric(xs_end["xL"]), as.numeric(xs_end["xk"]),
      as.numeric(xs_end["xl"]), as.numeric(xs_end["xU"]),
      tol_regime_switch = tol_regime_switch
    )
    .record_step(z_end, as.numeric(ans_broyden$fvec), sol_end)
  }
  
  rec_env$stage <- "broyden"
  
  # Decide whether to switch to Newton
  do_switch <- isTRUE(switch_to_newton)
  if (do_switch) {
    fmax   <- max(abs(ans_broyden$fvec %||% fn(ans_broyden$x)))
    iters  <- ans_broyden$iter %||% 0L
    do_switch <- is.finite(fmax) && (fmax <= switch_ftol) && (iters >= switch_iter_min)
  }
  
  # -------------------- Phase 2: pure Newton (optional) --------------------
  if (do_switch) {
    rec_env$stage <- "newton"
    z_start <- ans_broyden$x
    jf <- function(z) jac_fd(z)  # supply Jacobian explicitly
    ans_newton <- tryCatch(
      nleqslv::nleqslv(
        x = z_start, fn = fn, jac = jf, method = "Newton",
        global = "none", xscalm = "auto", control = newton_control
      ),
      error = function(e) e
    )
    
    # If Newton succeeds and improves max|f|, take it; otherwise keep Broyden result
    if (!inherits(ans_newton, "error")) {
      f_b <- max(abs(ans_broyden$fvec %||% fn(ans_broyden$x)))
      f_n <- max(abs(ans_newton$fvec  %||% fn(ans_newton$x)))
      final_ans <- if (is.finite(f_n) && (f_n <= f_b)) ans_newton else ans_broyden
      chosen <- if (identical(final_ans, ans_newton)) "newton" else "broyden"
    } else {
      final_ans <- ans_broyden
      chosen <- "broyden"
      warning(sprintf("Newton phase skipped (error: %s). Keeping Broyden solution.", conditionMessage(ans_newton)))
    }
  } else {
    final_ans <- ans_broyden
    chosen <- "broyden"
  }
  
  # -------------------- Build final suboptimal H and return --------------------
  xs_hat <- x_from_z(final_ans$x)
  sol_subopt <- build_suboptimal_H(
    p, xs_hat[["xL"]], xs_hat[["xk"]], xs_hat[["xl"]], xs_hat[["xU"]],
    tol_regime_switch = tol_regime_switch
  )
  
  hist_df <- NULL
  if (isTRUE(record_iterates)) {
    if (length(rec_env$hist)) {
      hist_list <- lapply(seq_along(rec_env$hist), function(k) {
        h <- rec_env$hist[[k]]
        data.frame(
          step = k,
          stage = as.character(h$stage),
          xL = as.numeric(h$xL),
          xk = as.numeric(h$xk),
          xl = as.numeric(h$xl),
          xU = as.numeric(h$xU),
          gamma = as.numeric(h$gamma),
          fmax  = as.numeric(h$fmax),
          Hp_xL  = as.numeric(h$Hp_xL),
          dHp_xk = as.numeric(h$dHp_xk),
          dHp_xl = as.numeric(h$dHp_xl),
          Hp_xU  = as.numeric(h$Hp_xU),
          stringsAsFactors = FALSE
        )
      })
      
      hist_df <- do.call(rbind.data.frame, hist_list)
      rownames(hist_df) <- NULL
    } else {
      # keep as a 0-row DF WITH columns so downstream checks are meaningful
      hist_df <- data.frame(
        step=integer(), stage=character(),
        xL=numeric(), xk=numeric(), xl=numeric(), xU=numeric(),
        gamma=numeric(), fmax=numeric(),
        Hp_xL=numeric(), dHp_xk=numeric(), dHp_xl=numeric(), Hp_xU=numeric(),
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    H      = sol_subopt$H,
    Hp     = sol_subopt$Hp,
    gamma  = sol_subopt$gamma,
    pieces = sol_subopt$pieces,
    params = sol_subopt$params,
    x      = sol_subopt$x,
    nleqslv_primary = ans_broyden,
    nleqslv_refined = if (exists("ans_newton")) ans_newton else NULL,
    chosen_solver   = chosen,
    iter_history = hist_df
  )
}


# ================== Ambiguity thresholds with fixed reflecting barriers =========
# Fix xL and xU (the controller's reflecting barriers), and solve only for the
# ambiguity thresholds xk and xl under the current model.
z_from_x_fixed_barriers <- function(xL, xk, xl, xU) {
  if (!(xL < xk && xk < xU && xk < xl)) {
    stop("Initial guess must satisfy xL < xk < xU and xk < xl (xl may be > xU).")
  }
  span <- xU - xL
  s <- (xk - xL) / span
  c(
    qlogis(s),
    log(xl - xk)
  )
}

x_from_z_fixed_barriers <- function(z, xL, xU) {
  stopifnot(length(z) == 2, xL < xU)
  xk <- xL + (xU - xL) * plogis(z[1])
  xl <- xk + exp(z[2])
  c(xL = xL, xk = xk, xl = xl, xU = xU)
}

ambiguity_threshold_residuals_fixed_barriers <- function(
    z, p, xL, xU, tol_build = 1e-10, tol_regime_switch = 1e-3
) {
  xs  <- x_from_z_fixed_barriers(z, xL, xU)
  sol <- build_suboptimal_H(
    p, xs["xL"], xs["xk"], xs["xl"], xs["xU"],
    tol_regime_switch = tol_regime_switch
  )

  xk <- xs["xk"]; xl <- xs["xl"]
  R1 <- sol$pieces$region1
  R2 <- sol$pieces$region2

  H1p <- function(x) gp_eval(R1$lam$lam_minus, x, R1$anchor_minus) * R1$u_minus +
    gp_eval(R1$lam$lam_plus , x, R1$anchor_plus ) * R1$u_plus  +
    qprime_eval(R1$q, x)

  H2p <- function(x) gp_eval(R2$lam$lam_minus, x, R2$anchor_minus) * R2$u_minus +
    gp_eval(R2$lam$lam_plus , x, R2$anchor_plus ) * R2$u_plus  +
    qprime_eval(R2$q, x)

  r1 <- H1p(xk) - H2p(xk)

  if (!isTRUE(sol$regime2)) {
    R3 <- sol$pieces$region3
    H3p <- function(x) gp_eval(R3$lam$lam_minus, x, R3$anchor_minus) * R3$u_minus +
      gp_eval(R3$lam$lam_plus , x, R3$anchor_plus ) * R3$u_plus  +
      qprime_eval(R3$q, x)

    r2 <- H2p(xl) - H3p(xl)
  } else {
    r2 <- IH_at_xlambda_stable(sol, xla = xl)
  }

  c(r1, r2)
}

solve_ambiguity_thresholds_fixed_barriers <- function(
    p, xL, xU, xk0, xl0,
    method = "Broyden",
    control = list(ftol=1e-11, xtol=1e-11, maxit=1e6, stepmax=1, trace=0),
    tol_build = 1e-10,
    verbose = FALSE,
    tol_regime_switch = 1e-3,
    switch_to_newton = FALSE,
    switch_ftol = 5e-2,
    switch_iter_min = 2L,
    newton_control = list(ftol=1e-12, xtol=1e-12, maxit=1e5, stepmax=1, trace=0),
    fd_step = 1e-6
) {
  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Please install.packages('nleqslv')")
  }
  if (!(xL < xk0 && xk0 < xU && xk0 < xl0)) {
    stop("Initial guess must satisfy xL < xk0 < xU and xk0 < xl0 (xl0 may be > xU).")
  }

  z0 <- z_from_x_fixed_barriers(xL, xk0, xl0, xU)

  fn <- function(z) {
    r <- ambiguity_threshold_residuals_fixed_barriers(
      z, p, xL, xU,
      tol_build = tol_build,
      tol_regime_switch = tol_regime_switch
    )
    if (isTRUE(verbose)) {
      xs <- x_from_z_fixed_barriers(z, xL, xU)
      cat(sprintf(
        "xL=%.6g  xk=%.6g  xl=%.6g  xU=%.6g  |  r=[% .3e, % .3e]\n",
        xs["xL"], xs["xk"], xs["xl"], xs["xU"], r[1], r[2]
      ))
    }
    r
  }

  jac_fd <- function(z, fz = NULL) {
    n <- length(z); J <- matrix(0.0, n, n)
    if (is.null(fz)) fz <- fn(z)
    for (j in 1:n) {
      h <- fd_step * (1 + abs(z[j]))
      zp <- z; zm <- z
      zp[j] <- z[j] + h;  zm[j] <- z[j] - h
      fp <- fn(zp);       fm <- fn(zm)
      J[, j] <- (fp - fm) / (2*h)
    }
    J
  }

  ans_broyden <- nleqslv::nleqslv(
    x = z0, fn = fn, method = method,
    global = "dbldog", xscalm = "auto", control = control
  )

  do_switch <- isTRUE(switch_to_newton)
  if (do_switch) {
    fmax  <- max(abs(ans_broyden$fvec %||% fn(ans_broyden$x)))
    iters <- ans_broyden$iter %||% 0L
    do_switch <- is.finite(fmax) && (fmax <= switch_ftol) && (iters >= switch_iter_min)
  }

  if (do_switch) {
    z_start <- ans_broyden$x
    jf <- function(z) jac_fd(z)
    ans_newton <- tryCatch(
      nleqslv::nleqslv(
        x = z_start, fn = fn, jac = jf, method = "Newton",
        global = "none", xscalm = "auto", control = newton_control
      ),
      error = function(e) e
    )

    if (!inherits(ans_newton, "error")) {
      f_b <- max(abs(ans_broyden$fvec %||% fn(ans_broyden$x)))
      f_n <- max(abs(ans_newton$fvec  %||% fn(ans_newton$x)))
      final_ans <- if (is.finite(f_n) && (f_n <= f_b)) ans_newton else ans_broyden
      chosen <- if (identical(final_ans, ans_newton)) "newton" else "broyden"
    } else {
      final_ans <- ans_broyden
      chosen <- "broyden"
      warning(sprintf("Newton phase skipped (error: %s). Keeping Broyden solution.", conditionMessage(ans_newton)))
    }
  } else {
    final_ans <- ans_broyden
    chosen <- "broyden"
  }

  xs_hat <- x_from_z_fixed_barriers(final_ans$x, xL, xU)
  sol_subopt <- build_suboptimal_H(
    p, xs_hat[["xL"]], xs_hat[["xk"]], xs_hat[["xl"]], xs_hat[["xU"]],
    tol_regime_switch = tol_regime_switch
  )

  list(
    H      = sol_subopt$H,
    Hp     = sol_subopt$Hp,
    gamma  = sol_subopt$gamma,
    pieces = sol_subopt$pieces,
    params = sol_subopt$params,
    x      = sol_subopt$x,
    nleqslv_primary = ans_broyden,
    nleqslv_refined = if (exists("ans_newton")) ans_newton else NULL,
    chosen_solver   = chosen
  )
}



`%||%` <- function(a, b) if (!is.null(a)) a else b
