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
opt_conditions_residuals <- function(z, p, tol_build=1e-10, tol_regime_switch = 1e-3) {
  xs  <- x_from_z(z)
  sol <- build_suboptimal_H(p, xs["xL"], xs["xk"], xs["xl"], xs["xU"],
                            tol_regime_switch = tol_regime_switch)
  
  xL <- xs["xL"]; xk <- xs["xk"]; xl <- xs["xl"]; xU <- xs["xU"]
  R1 <- sol$pieces$region1; R2 <- sol$pieces$region2
  
  H1p <- function(x) gp_eval(R1$lam$lam_minus,x,R1$anchor_minus)*R1$u_minus +
    gp_eval(R1$lam$lam_plus ,x,R1$anchor_plus )*R1$u_plus  +
    qprime_eval(R1$q, x)
  
  H2p <- function(x) gp_eval(R2$lam$lam_minus,x,R2$anchor_minus)*R2$u_minus +
    gp_eval(R2$lam$lam_plus ,x,R2$anchor_plus )*R2$u_plus  +
    qprime_eval(R2$q, x)
  
  if (!isTRUE(sol$regime2)) {
    # Regime 1: your original 4 smooth-fit equations (needs Region 3)
    R3 <- sol$pieces$region3
    H3p <- function(x) gp_eval(R3$lam$lam_minus,x,R3$anchor_minus)*R3$u_minus +
      gp_eval(R3$lam$lam_plus ,x,R3$anchor_plus )*R3$u_plus  +
      qprime_eval(R3$q, x)
    
    r1 <- H1p(xL)                 # (52)
    r2 <- H3p(xU)                 # (53)
    r3 <- H1p(xk) - H2p(xk)       # (54)
    r4 <- H2p(xl) - H3p(xl)       # (55)
    c(r1, r2, r3, r4)
  } else {
    # Regime 2: no Region 3 inside band. Use:
    #   (i) smooth-fit at xL, xU, xk
    #   (ii) definition of xλ via IH(xλ)=0 (now xλ may be > xU)
    r1 <- H1p(xL)
    r2 <- H2p(xU)
    r3 <- H1p(xk) - H2p(xk)
    r4 <- IH_at_xlambda_stable(sol, xla = xl)
    c(r1, r2, r3, r4)
  }
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
    fd_step = 1e-6               # relative step for finite-difference Jacobian
) {
  if (!requireNamespace("nleqslv", quietly=TRUE))
    stop("Please install.packages('nleqslv')")
  if (!(xL0 < xk0 && xk0 < xl0 && xk0 < xU0))
    stop("Initial guess must satisfy xL0 < xk0, and xk0 < xl0, xk0 < xU0 (xl0 may be > xU0).")
  
  # ---- reparam z <-> x (keeps barriers ordered) ----
  z0 <- z_from_x(xL0, xk0, xl0, xU0)
  
  # Residuals (Eqs. (52)–(55) in stable form)
  fn <- function(z) {
    r <- opt_conditions_residuals(z, p, tol_build=tol_build, tol_regime_switch = tol_regime_switch)
    if (isTRUE(verbose)) {
      xs <- x_from_z(z)
      cat(sprintf(
        "xL=%.6g  xk=%.6g  xl=%.6g  xU=%.6g  |  r=[% .3e, % .3e, % .3e, % .3e]\n",
        xs["xL"], xs["xk"], xs["xl"], xs["xU"], r[1], r[2], r[3], r[4]
      ))
    }
    r
  }
  
  # Stable, central-difference Jacobian (O(h^2))
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
  
  # -------------------- Phase 1: robust Broyden --------------------
  ans_broyden <- nleqslv::nleqslv(
    x = z0, fn = fn, method = method,
    global = "dbldog", xscalm = "auto", control = control
  )
  
  # Decide whether to switch to Newton
  do_switch <- isTRUE(switch_to_newton)
  if (do_switch) {
    fmax   <- max(abs(ans_broyden$fvec %||% fn(ans_broyden$x)))
    iters  <- ans_broyden$iter %||% 0L
    do_switch <- is.finite(fmax) && (fmax <= switch_ftol) && (iters >= switch_iter_min)
  }
  
  # -------------------- Phase 2: pure Newton (optional) --------------------
  if (do_switch) {
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
    p, xs_hat[["xL"]], xs_hat[["xk"]], xs_hat[["xl"]], xs_hat[["xU"]]
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


# ================================ Figures =======================================
.with_plot_margins <- function(expr) { 
  op <- par(mar=c(2,2,0.5,0.5)+0.2)
  on.exit(par(op), add=TRUE); force(expr)
}

render_and_save <- function(
    fname_base, plotfun, show=interactive(), save=TRUE,
    dir="figures",
    width_in = 6.67, height_in = 4.67, dpi = 300,
    pointsize = NULL, family = NULL,
    match_current = FALSE, use_cairo_pdf = TRUE
) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  
  # Derive device settings (size, ps, family) from current device when available
  if (match_current && dev.cur() != 1L) {
    sz <- dev.size("in")
    width_in  <- sz[1]; height_in <- sz[2]
    if (is.null(pointsize)) pointsize <- par("ps")
    if (is.null(family))    family    <- par("family")
  }
  if (is.null(pointsize)) pointsize <- 12
  fam <- if (is.null(family) || !nzchar(family)) "sans" else family
  
  # SHOW: draw on the current device but force the same family for consistency
  if (isTRUE(show)) {
    op <- par(family = fam)               # temporary — reverts on exit
    on.exit(par(op), add = TRUE)
    .with_plot_margins(plotfun())
  }
  
  if (!isTRUE(save)) return(invisible())
  
  # PNG: exact physical size + pointsize
  png(file.path(dir, paste0(fname_base, ".png")),
      width = width_in, height = height_in, units = "in",
      res = dpi, pointsize = pointsize,
      type = getOption("bitmapType", "cairo"), antialias = "subpixel")
  .with_plot_margins(plotfun()); dev.off()
  
  # PDF: same physical size + pointsize + family
  if (isTRUE(use_cairo_pdf) && capabilities("cairo")) {
    cairo_pdf(file.path(dir, paste0(fname_base, ".pdf")),
              width = width_in, height = height_in,
              pointsize = pointsize, family = fam)
  } else {
    pdf(file.path(dir, paste0(fname_base, ".pdf")),
        width = width_in, height = height_in,
        pointsize = pointsize, family = fam, useDingbats = FALSE)
  }
  .with_plot_margins(plotfun()); dev.off()
}

plot_H <- function(sol, params, show=interactive(), save=TRUE, n = 600,
                   pad_x_frac = 0.10, top_blank = 0.10, bottom_blank = 0.05,
                   out_base="H_thresholds", dir="figures",
                   width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                   margins = c(3.2, 3.25, 1.05, 0.9),
                   axis_mgp = c(2.2, 0.7, 0),
                   tcl = -0.25,
                   show_x_axis_title = FALSE,
                   show_y_axis_title = TRUE,
                   x_axis_title = expression(x),
                   y_axis_title = expression(H(x) == V*minute*(x)),
                   show_tick_labels = TRUE,
                   tick_cex = 1.7,
                   axis_title_cex = 1,
                   base_cex = 1, mex = 1) {
  
  xs   <- unlist(sol$x)[c("xL","xk","xl","xU")]
  xL <- xs[1]; xk <- xs[2]; xl <- xs[3]; xU <- xs[4]
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22"); ltys <- c(1,2,2,1)
  
  x_max <- max(xU, xl)
  W <- x_max - xL
  x_from <- xL - pad_x_frac*W
  x_to   <- x_max + pad_x_frac*W
  
  frac_band <- 1 - top_blank - bottom_blank
  ylim_y <- range(sol$H(seq(x_from, x_to, l = n)))
  ylim_y <- ylim_y + diff(ylim_y)/frac_band * c(-bottom_blank, top_blank)
  
  plotfun <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add=TRUE)
    
    par(xaxs="i",
        oma = c(0,0,0,0),
        cex = base_cex,
        mex = mex,
        mar = margins,
        mgp = axis_mgp,
        tcl = tcl,
        cex.axis = tick_cex,
        cex.lab = axis_title_cex)
    
    curve(sol$H(x), from=x_from, to=x_to, n=n,
          xlab=if (show_x_axis_title) x_axis_title else "",
          ylab=if (show_y_axis_title) y_axis_title else "",
          col="black", lwd=1.4, ylim=ylim_y,
          xaxt="n", yaxt="n", bty="n")
    
    axis(1, labels = show_tick_labels)
    axis(2, labels = show_tick_labels)
    
    segments(x0=xs, y0=rep(-params$u, length(xs)), x1=xs, y1=rep(params$l, length(xs)),
             lty=ltys, col=cols, lwd=2)
    abline(h=-params$u, lty=3, col="#777777")
    abline(h=0,       lty=3, col="#BBBBBB")
    abline(h=params$l, lty=4, col="#777777")
    
    inch_to_user <- diff(grconvertX(c(0,1), from="in", to="user"))
    col_w <- 0.23 * inch_to_user
    tw <- c(0.36 * inch_to_user, rep(col_w, 5), 0.38 * inch_to_user)
    
    legend("topleft",
           legend=expression(H(x), underline(x), x^kappa, x^lambda, bar(x), -c[U], c[D]),
           lty=c(1,1,2,2,1,3,4),
           col=c("black", cols, "#777777", "#777777"),
           lwd=c(1.4, rep(2,6)),
           bty="n", horiz=TRUE, seg.len=1.5, x.intersp=0.5,
           text.width=tw)
    
    box()
  }
  
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}


plot_H_prime <- function(sol, params, show=interactive(), save=TRUE, n = 800,
                         pad_x_frac = 0.10, top_blank = 0.10, bottom_blank = 0.05,
                         out_base="Hprime_thresholds", dir="figures",
                         width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                         show_x_axis_title = FALSE,
                         show_y_axis_title = TRUE,
                         x_axis_title = expression(x),
                         y_axis_title = expression(H*minute*(x) == W(x)),
                         axis_title_cex = 1) {
  stopifnot(!is.null(sol$Hp))
  
  xs   <- unlist(sol$x)[c("xL","xk","xl","xU")]
  xL <- xs[1]; xk <- xs[2]; xl <- xs[3]; xU <- xs[4]
  
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22")  # xL,xk,xl,xU
  ltys <- c(1,2,2,1)
  
  x_max <- max(xU, xl)
  W <- x_max - xL
  x_from <- xL - pad_x_frac*W
  x_to   <- x_max + pad_x_frac*W
  
  # sample grid and break the line at thresholds to avoid artificial connections
  xg <- seq(x_from, x_to, length.out = n)
  yg <- sol$Hp(xg)
  
  # Force line breaks around the 4 special points (xL,xk,xl,xU)
  for (b in xs) {
    j <- which.min(abs(xg - b))
    for (k in c(j-1, j, j+1)) {
      if (k >= 1 && k <= length(yg)) yg[k] <- NA_real_
    }
  }
  
  # y-limits with controlled blank margins (like plot_H)
  frac_band <- 1 - top_blank - bottom_blank
  ylim_y <- range(yg, finite = TRUE)
  if (!all(is.finite(ylim_y))) ylim_y <- c(-1, 1)
  if (diff(ylim_y) == 0) {
    bump <- max(1, abs(ylim_y[1]))
    ylim_y <- ylim_y + c(-1, 1) * 0.05 * bump
  }
  ylim_y <- ylim_y + diff(ylim_y)/frac_band * c(-bottom_blank, top_blank)
  
  plotfun <- function() {
    op <- par(xaxs="i", cex.lab = axis_title_cex); on.exit(par(op), add=TRUE)
    
    plot(xg, yg, type="l",
         xlab=if (show_x_axis_title) x_axis_title else "",
         ylab=if (show_y_axis_title) y_axis_title else "",
         col="black", lwd=1.4, ylim=ylim_y)
    
    # thresholds as vertical segments spanning the panel
    segments(x0=xs, y0=rep(ylim_y[1], length(xs)),
             x1=xs, y1=rep(ylim_y[2], length(xs)),
             lty=ltys, col=cols, lwd=2)
    
    # reference line at 0
    abline(h=0, lty=3, col="#BBBBBB")
    
    # legend (same style as plot_H)
    inch_to_user <- diff(grconvertX(c(0,1), from="in", to="user"))
    col_w <- 0.23 * inch_to_user
    tw <- c(0.36 * inch_to_user,
            rep(col_w, 5),
            0.20 * inch_to_user)
    
    legend("topleft",
           legend = expression(H*minute*(x), underline(x), x^kappa, x^lambda, bar(x), 0),
           lty    = c(1, 1, 2, 2, 1, 3),
           col    = c("black", cols, "#BBBBBB"),
           lwd    = c(1.4, rep(2,4), 1.2),
           bty="n", horiz=TRUE, seg.len=1.5, x.intersp=0.5,
           text.width=tw)
    
    box()
  }
  
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

# --------------------------- Reflected process sim ------------------------------
simulate_reflected_jd <- function(T=8, dt=0.001, x0=NULL, params, thresholds,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- as.integer(round(T/dt)); tt <- seq(0, T, length.out=n+1)
  xL <- thresholds$xL; xU <- thresholds$xU
  if (is.null(x0)) x0 <- (thresholds$xL + thresholds$xU)/2
  EY <- -1/params$mu; drift <- params$b - params$r * EY
  X <- numeric(n+1); X[1] <- x0; U_proc <- L_proc <- numeric(n+1)
  jumped <- logical(n+1); dJ <- numeric(n+1)
  for (i in 1:n) {
    dW <- sqrt(dt)*rnorm(1); Nj <- rpois(1, params$r*dt)
    dJ[i+1] <- if (Nj>0) -sum(rexp(Nj, rate=params$mu)) else 0
    x_star <- X[i] + drift*dt + params$sigma*dW + dJ[i+1]
    U_inc <- L_inc <- 0
    if (x_star < xL) { U_inc <- xL - x_star; X[i+1] <- xL }
    else if (x_star > xU) { L_inc <- x_star - xU; X[i+1] <- xU }
    else X[i+1] <- x_star
    U_proc[i+1] <- U_proc[i] + U_inc; L_proc[i+1] <- L_proc[i] + L_inc
    jumped[i+1] <- Nj > 0
  }
  list(time=tt, X=X, U=U_proc, L=L_proc, jumped=jumped, dJ=dJ)
}

# --------------------------- Small reusable drawer for the middle panel --------
# Draw the reflected trajectory with thresholds; returns the ylim used (so panels
# can be tuned consistently if needed).
.draw_reflected_panel <- function(t, sim, xL, xU, xk, xl,
                                  top_blank=0.075, bottom_blank=0.05,
                                  draw_legend=TRUE) {
  Hcorr <- xU - xL
  frac  <- 1 - top_blank - bottom_blank
  Htot  <- Hcorr/frac
  ylim  <- c(xL - Htot*bottom_blank, xU + Htot*top_blank)
  
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22") # xL, xk, xl, xU
  ltys <- c(1,2,2,1)
  
  plot(t, sim, type="n", xlab="", ylab="", ylim=ylim, xaxt="n")
  usr <- par("usr")
  rect(usr[1], xL, usr[2], xU, col=adjustcolor("gray85", 0.6), border=NA)
  
  lines(t, sim, lwd=1.2)
  abline(h=c(xL, xU), col=c(cols[1], cols[4]), lty=c(ltys[1], ltys[4]), lwd=2)
  abline(h=c(xk, xl), col=c(cols[2], cols[3]), lty=c(ltys[2], ltys[3]), lwd=2)
  
  if (draw_legend) {
    inch_to_user <- diff(grconvertX(c(0,1), from="in", to="user"))
    col_w <- 0.23 * inch_to_user
    tw <- c(0.36 * inch_to_user, rep(col_w, 4))  # 5 entries total
    
    legend("topleft",
           legend = expression(bar(X)[t], underline(x), x^kappa, x^lambda, bar(x)),
           lty    = c(1, 1, 2, 2, 1),                # 5 entries
           col    = c("black", cols),                # 1 + 4 entries
           lwd    = c(1.2, rep(2,4)),
           bty    = "n", horiz = TRUE, x.intersp = 0.5, seg.len = 2,
           text.width = tw
    )
  }
  box()
  invisible(ylim)
}

# --------------------------- Original single plots (kept, just call helper) ----
plot_reflected_jd <- function(sol, params, sim=NULL, seed=123, show=interactive(), 
                              save=TRUE, out_base="reflected_jd", dir="figures",
                              width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                              top_blank=0.075, bottom_blank=0.05,
                              show_x_axis_title = TRUE,
                              show_y_axis_title = FALSE,
                              x_axis_title = "time",
                              y_axis_title = expression(bar(X)[t]),
                              axis_title_cex = 1) {
  xs <- sol$x
  if (is.null(sim)) { set.seed(seed); sim <- simulate_reflected_jd(params=params, thresholds=xs) }
  t <- sim$time; X <- sim$X; xL <- xs$xL; xU <- xs$xU; xk <- xs$xk; xl <- xs$xl
  plotfun <- function() {
    .draw_reflected_panel(t, X, xL, xU, xk, xl, top_blank=top_blank, bottom_blank=bottom_blank)
    axis(1)
    if (show_x_axis_title || show_y_axis_title) {
      title(xlab = if (show_x_axis_title) x_axis_title else "",
            ylab = if (show_y_axis_title) y_axis_title else "",
            cex.lab = axis_title_cex)
    }
  }
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_controls <- function(sim, show=interactive(), save=TRUE,
                          out_base="singular_controls", dir="figures",
                          width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                          axis_title_cex = 1) {
  plotfun <- function() {
    rng <- range(sim$U, sim$L)
    plot(sim$time, sim$L, type="s", xlab="time", ylab="cumulative push",
         lwd=1.6, ylim=rng, lty=2, cex.lab = axis_title_cex)
    lines(sim$time, sim$U, type="s", lwd=1.6)
    legend("topleft",
           legend=c(expression(D[t]~"(pushes down)"), expression(U[t]~"(pushes up)")),
           lty=c(2,1), lwd=2, bty="n")
    box()
  }
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

# --------------------------- NEW: merged 3-panel figure -------------------------
# Top:   L_t (pushes down)
# Middle: reflected process + thresholds
# Bottom: U_t (pushes up)
plot_reflected_with_controls <- function(sol, params, sim=NULL, seed=123,
                                         show=interactive(), save=TRUE,
                                         out_base="reflected_with_controls",
                                         dir="figures", draw_legend = TRUE,
                                         width_in = 6.67, height_in = 4.67, dpi = 300, match_current = FALSE,
                                         heights = c(0.7, 1.7, 0.7),
                                         top_blank=0.075, bottom_blank=0.05,
                                         margins = c(3.2, 3.25, 1.05, 0.9),
                                         axis_mgp = c(2.2, 0.7, 0),
                                         tcl = -0.25,
                                         show_x_axis_title = TRUE,
                                         show_y_axis_title = TRUE,
                                         x_axis_title = "time",
                                         y_axis_title = list(top = expression(L[t]),
                                                             bottom = expression(U[t])),
                                         show_tick_labels = TRUE,
                                         tick_cex = 1.7,
                                         axis_title_cex = 1,
                                         base_cex = 1, mex = 1) { 
  xs <- sol$x
  if (is.null(sim)) { set.seed(seed); sim <- simulate_reflected_jd(params=params, thresholds=xs) }
  t <- sim$time; X <- sim$X; xL <- xs$xL; xU <- xs$xU; xk <- xs$xk; xl <- xs$xl
  
  if (is.list(y_axis_title)) {
    y_axis_title_top <- y_axis_title[["top"]]
    y_axis_title_bottom <- y_axis_title[["bottom"]]
    if (is.null(y_axis_title_top) && length(y_axis_title) >= 1) {
      y_axis_title_top <- y_axis_title[[1]]
    }
    if (is.null(y_axis_title_bottom)) {
      y_axis_title_bottom <- y_axis_title_top
    }
  } else {
    y_axis_title_top <- y_axis_title
    y_axis_title_bottom <- y_axis_title
  }
  
  bottom <- margins[1]; left <- margins[2]; top <- margins[3]; right <- margins[4]
  
  plotfun <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add=TRUE)
    
    par(oma = c(0,0,0,0), cex = base_cex, mex = mex)
    
    layout(matrix(1:3, nrow=3), heights = heights)
    
    # --- TOP: D_t
    par(mar = c(0.4, left, top, right),
        mgp=axis_mgp, tcl=tcl,
        cex = base_cex, mex = mex,
        cex.axis=tick_cex, cex.lab = axis_title_cex)
    plot(t, sim$L, type="s",
         xlab="", ylab=if (show_y_axis_title) y_axis_title_top else "",
         xaxt="n", yaxt="n",
         lwd=2, col = "#228B22")
    axis(2, labels = show_tick_labels)
    legend("topleft", legend = expression(D[t]), cex = 1,
           lty = 1, lwd = 2, col = "#228B22", bty = "n")
    box()
    
    # --- MIDDLE
    par(mar = c(0.2, left, 0.2, right),
        mgp=axis_mgp, tcl=tcl,
        cex = base_cex, mex = mex,
        cex.axis=tick_cex, cex.lab = axis_title_cex)
    .draw_reflected_panel(t, X, xL, xU, xk, xl,
                          top_blank=top_blank, bottom_blank=bottom_blank,
                          draw_legend=draw_legend)
    
    # --- BOTTOM: U_t
    par(mar = c(bottom, left, 0.4, right),
        mgp=axis_mgp, tcl=tcl,
        cex = base_cex, mex = mex,
        cex.axis=tick_cex, cex.lab = axis_title_cex)
    plot(t, sim$U, type="s",
         xlab=if (show_x_axis_title) x_axis_title else "",
         ylab=if (show_y_axis_title) y_axis_title_bottom else "",
         xaxt="n", yaxt="n",
         lwd=2, col = "#B22222")
    axis(1, labels = show_tick_labels)
    axis(2, labels = show_tick_labels)
    legend("topleft", legend = expression(U[t]), cex = 1,
           lty = 1, lwd = 2, col = "#B22222", bty = "n")
    box()
  }
  
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}


# ============================= Comparative sweeps ================================
comparative_sweeper <- function(
    sweep_param = "b",
    sweep_values = seq(-10, 10, by = 0.1),
    b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1.0, l = 1.0,
    xL0 = -0.5, xk0 = -0.1, xl0 = 0.3, xU0 = 1.0,
    out_dir  = file.path("figures", "sweeps"),
    out_name = NULL,
    save = FALSE
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  sweep_param <- as.character(sweep_param)
  # NEW: allow "inv_mu" as synthetic parameter = 1/mu
  known_params <- c("b","delta","r","eps","sigma","mu","inv_mu","u","l")
  if (!sweep_param %in% known_params) {
    stop(sprintf("sweep_param must be one of: %s", paste(known_params, collapse = ", ")))
  }
  fixed_params <- list(b = b, delta = delta, r = r, eps = eps, sigma = sigma, mu = mu, u = u, l = l)
  
  # Metrics: thresholds + ergodic value γ
  res <- data.frame(
    sweep_param = sweep_param,
    sweep_value = as.numeric(sweep_values),
    xL = NA_real_, xk = NA_real_, xl = NA_real_, xU = NA_real_,
    gamma = NA_real_,
    converged = FALSE, error_msg = NA_character_,
    stringsAsFactors = FALSE
  )
  last_sol <- NULL
  
  for (i in seq_along(sweep_values)) {
    val <- sweep_values[i]
    message(sprintf("[sweep %s] %2d/%d  %s = % .6g ...",
                    sweep_param, i, length(sweep_values), sweep_param, val)); flush.console()
    
    p <- fixed_params
    
    # NEW: map 1/mu ("inv_mu") to actual mu > 0
    if (sweep_param == "inv_mu") {
      if (val <= 0) {
        res$converged[i] <- FALSE
        res$error_msg[i] <- "inv_mu (1/mu) must be > 0"
        next
      }
      p$mu <- 1/val
    } else {
      p[[sweep_param]] <- val
    }
    
    params <- do.call(make_params, p)
    
    # warm-start from last solution if available
    xL_init <- xL0; xk_init <- xk0; xl_init <- xl0; xU_init <- xU0
    if (!is.null(last_sol)) {
      xL_init <- last_sol$x$xL; xk_init <- last_sol$x$xk
      xl_init <- last_sol$x$xl; xU_init <- last_sol$x$xU
    }
    
    sol_opt <- tryCatch(
      solve_optimal_barriers(params, xL_init, xk_init, xl_init, xU_init, verbose = FALSE),
      error = function(e) e
    )
    
    if (inherits(sol_opt, "error")) {
      res$converged[i] <- FALSE
      res$error_msg[i] <- conditionMessage(sol_opt)
    } else {
      th <- sol_opt$x
      res[i, c("xL","xk","xl","xU")] <- unlist(th[c("xL","xk","xl","xU")])
      res$gamma[i] <- sol_opt$gamma
      res$converged[i] <- TRUE
      res$error_msg[i] <- NA_character_
      last_sol <- sol_opt
    }
  }
  
  # Default name (used only for CSV; plotting can choose its own if desired)
  if (is.null(out_name)) {
    step_str <- if (length(sweep_values) > 1) {
      diffs <- unique(round(diff(sweep_values), 10))
      if (length(diffs) == 1) sprintf("by_%g", diffs) else "custom_steps"
    } else "single_value"
    out_name <- sprintf("thresholds_vs_%s_%g_to_%g_%s",
                        sweep_param, min(sweep_values), max(sweep_values), step_str)
  }
  
  if (save) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    # NEW: build a filename suffix with all NON-swept parameter values
    fixed_for_name <- fixed_params
    if (sweep_param %in% names(fixed_for_name)) fixed_for_name[[sweep_param]] <- NULL
    if (sweep_param == "inv_mu") fixed_for_name$mu <- NULL  # mu is being varied via inv_mu
    
    fixed_tag <- paste(
      sprintf("%s_%s",
              names(fixed_for_name),
              vapply(fixed_for_name, function(x) sprintf("%.6g", x), "")
              
      ),
      collapse = "__"
    )
    
    csv_path <- file.path(out_dir, paste0(out_name, "__", fixed_tag, ".csv"))
    write.csv(res, csv_path, row.names = FALSE)
    message(sprintf("Saved sweep metrics to: %s", normalizePath(csv_path)))
  }
  
  # Return metrics object (no plotting)
  list(
    results     = res,
    out_dir     = out_dir,
    out_name    = out_name,
    fixed_params = fixed_params
  )
}


plot_sweep <- function(
    sweep_obj,
    out_dir  = file.path("figures", "sweeps"),
    out_name = NULL,
    show = interactive(),
    save = FALSE,
    width_in = 6.67,
    height_in = 4.67,
    dpi = 300,
    match_current = FALSE,
    cols = c("#B22222", "#1E90FF", "#8A2BE2", "#228B22"),
    ltys = c(1, 1, 1, 1),
    title = TRUE,
    show_x_axis_title = TRUE,
    show_y_axis_title = TRUE,
    x_axis_title = NULL,
    y_axis_title = NULL,
    axis_title_cex = 1,
    margins = c(2.2, 2.2, 1.5, 1.2),
    axis_mgp = c(3, 1, 0),
    plot_gamma = FALSE, 
    gamma_layout = c("stacked","separate")
) {
  gamma_layout <- match.arg(gamma_layout)
  bottom <- margins[1]; left <- margins[2]; top <- margins[3]; right <- margins[4]
  
  # Accept either the list returned by comparative_sweeper() or a bare data.frame
  res <- if (is.data.frame(sweep_obj)) {
    sweep_obj
  } else {
    sweep_obj$results %||% sweep_obj
  }
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  ok <- res$converged
  if (!any(ok)) {
    warning("No successful solutions; nothing to plot.")
    return(invisible(NULL))
  }
  
  sweep_param <- as.character(unique(res$sweep_param))
  if (length(sweep_param) != 1L) {
    warning("Multiple sweep_param values detected; using the first.")
    sweep_param <- sweep_param[1L]
  }
  
  sweep_values <- res$sweep_value
  
  # Default fname base if not provided (use raw sweep_param, not pretty label)
  if (is.null(out_name)) {
    step_str <- if (length(sweep_values) > 1) {
      diffs <- unique(round(diff(sweep_values), 10))
      if (length(diffs) == 1) sprintf("by_%g", diffs) else "custom_steps"
    } else "single_value"
    out_name <- sprintf("thresholds_vs_%s_%g_to_%g_%s",
                        sweep_param, min(sweep_values), max(sweep_values), step_str)
  }
  
  # append non-swept parameter values to filename (only if available)
  fixed_for_name <- NULL
  if (!is.data.frame(sweep_obj) && !is.null(sweep_obj$fixed_params)) {
    fixed_for_name <- sweep_obj$fixed_params
  }
  
  if (!is.null(fixed_for_name)) {
    if (sweep_param %in% names(fixed_for_name)) fixed_for_name[[sweep_param]] <- NULL
    if (sweep_param == "inv_mu") fixed_for_name$mu <- NULL
    
    fixed_tag <- paste(
      sprintf("%s_%s",
              names(fixed_for_name),
              vapply(fixed_for_name, function(x) sprintf("%.6g", x), "")
      ),
      collapse = "__"
    )
    fixed_tag <- gsub("\\s+", "", fixed_tag)
    out_name <- paste0(out_name, "__", fixed_tag)
  }
  
  # nice label for titles (e.g. "1/mu" instead of "inv_mu")
  param_label_str <- if (sweep_param == "inv_mu") "1/mu" else sweep_param
  default_x_axis_title <- switch(sweep_param,
                                 "b"      = expression(b),
                                 "delta"  = expression(delta),
                                 "r"      = expression(r),
                                 "eps"    = expression(epsilon),
                                 "sigma"  = expression(sigma),
                                 "mu"     = expression(mu),
                                 "inv_mu" = expression(1/mu),
                                 "u"      = expression(u),
                                 "l"      = expression(l),
                                 sweep_param)
  if (is.null(x_axis_title)) {
    x_axis_title <- default_x_axis_title
  }
  
  default_y_axis_title_thresholds <- "threshold value"
  default_y_axis_title_gamma <- expression(gamma)
  if (is.null(y_axis_title)) {
    y_axis_title_thresholds <- default_y_axis_title_thresholds
    y_axis_title_gamma <- default_y_axis_title_gamma
  } else if (is.list(y_axis_title)) {
    if (is.null(names(y_axis_title))) {
      y_axis_title_thresholds <- y_axis_title[[1]]
      y_axis_title_gamma <- if (length(y_axis_title) >= 2L) {
        y_axis_title[[2]]
      } else {
        y_axis_title[[1]]
      }
    } else {
      y_axis_title_thresholds <- y_axis_title[["thresholds"]]
      y_axis_title_gamma <- y_axis_title[["gamma"]]
      if (is.null(y_axis_title_thresholds)) {
        y_axis_title_thresholds <- default_y_axis_title_thresholds
      }
      if (is.null(y_axis_title_gamma)) {
        y_axis_title_gamma <- default_y_axis_title_gamma
      }
    }
  } else {
    y_axis_title_thresholds <- y_axis_title
    y_axis_title_gamma <- y_axis_title
  }
  
  # y-range for thresholds
  ylm_th <- range(as.numeric(unlist(res[ok, c("xL","xk","xl","xU")])), na.rm = TRUE)
  
  # y-range for gamma (if present)
  has_gamma <- "gamma" %in% names(res) && any(is.finite(res$gamma[ok]))
  if (!has_gamma) plot_gamma <- FALSE
  if (plot_gamma) {
    ylm_g <- range(res$gamma[ok], na.rm = TRUE)
  }
  
  title_str   <- if (title)
    sprintf("Ambiguity thresholds & barriers vs %s", param_label_str) else ""
  title_gamma <- if (title)
    sprintf("Ergodic value %s vs %s", "\u03b3", param_label_str) else ""
  
  ## ---------- Helper to draw thresholds-only panel ----------
  thresholds_panel <- function(show_x_title = TRUE, mar_override = NULL) {
    mar <- if (is.null(mar_override)) margins else mar_override
    op <- par(xaxs = "i", mar = mar, mgp = axis_mgp,
              cex.axis = 1.4, cex.lab = axis_title_cex)
    on.exit(par(op), add = TRUE)
    
    xlab <- if (show_x_title && show_x_axis_title) x_axis_title else ""
    plot(res$sweep_value, res$xL, type = "n", xlab = xlab,
         ylab = if (show_y_axis_title) y_axis_title_thresholds else "",
         ylim = ylm_th)
    grid()
    lines(res$sweep_value, res$xL, lwd = 2, lty = ltys[1], col = cols[1])
    lines(res$sweep_value, res$xk, lwd = 2, lty = ltys[2], col = cols[2])
    lines(res$sweep_value, res$xl, lwd = 2, lty = ltys[3], col = cols[3])
    lines(res$sweep_value, res$xU, lwd = 2, lty = ltys[4], col = cols[4])
    legend("topleft",
           legend = expression(underline(x), x^kappa, x^lambda, bar(x)),
           col = cols, lty = ltys, lwd = 2, bty = "n", horiz = TRUE,
           x.intersp = 0.5, seg.len = 2)
    title(title_str)
  }
  
  ## ---------- Helper to draw gamma-only panel ----------
  gamma_panel <- function(show_x_title = FALSE, show_x_axis_ticks = TRUE,
                          mar_override = NULL) {
    mar  <- if (is.null(mar_override)) margins else mar_override
    op   <- par(xaxs = "i", mar = mar, mgp = axis_mgp,
                cex.axis = 1.4, cex.lab = axis_title_cex)
    on.exit(par(op), add = TRUE)
    
    xlab <- if (show_x_title && show_x_axis_title) x_axis_title else ""
    xaxt <- if (show_x_axis_ticks) "s" else "n"
    
    plot(res$sweep_value, res$gamma, type = "l", lwd = 2,
         xlab = xlab,
         ylab = if (show_y_axis_title) y_axis_title_gamma else "",
         ylim = ylm_g, xaxt = xaxt)
    grid()
    legend("topleft", legend = expression(gamma), lwd = 2, bty = "n")
    title(title_gamma)
  }
  
  ## ======================= CASE 1: thresholds only =======================
  if (!plot_gamma) {
    plotfun <- function() {
      thresholds_panel(show_x_title = TRUE)
    }
    
    if (exists("render_and_save")) {
      render_and_save(fname_base = out_name, plotfun = plotfun,
                      show = show, save = save, dir = out_dir,
                      width_in = width_in, height_in = height_in, dpi = dpi,
                      match_current = match_current)
    } else {
      if (isTRUE(show)) plotfun()
      if (isTRUE(save)) {
        png(file.path(out_dir, paste0(out_name, ".png")),
            width = 1200, height = 800, res = 150); plotfun(); dev.off()
        pdf(file.path(out_dir, paste0(out_name, ".pdf")),
            width = 9, height = 6); plotfun(); dev.off()
      }
    }
    
    ## ======================= CASE 2: gamma + thresholds, separate =======================
  } else if (gamma_layout == "separate") {
    plotfun_th <- function() thresholds_panel(show_x_title = TRUE)
    plotfun_g  <- function() gamma_panel(show_x_title = FALSE,
                                         show_x_axis_ticks = TRUE)
    
    if (exists("render_and_save")) {
      render_and_save(fname_base = out_name,
                      plotfun = plotfun_th, show = show, save = save, dir = out_dir,
                      width_in = width_in, height_in = height_in, dpi = dpi,
                      match_current = match_current)
      render_and_save(fname_base = paste0(out_name, "_gamma"),
                      plotfun = plotfun_g, show = show, save = save, dir = out_dir,
                      width_in = width_in, height_in = height_in, dpi = dpi,
                      match_current = match_current)
    } else {
      if (isTRUE(show)) {
        plotfun_th()
        plotfun_g()
      }
      if (isTRUE(save)) {
        png(file.path(out_dir, paste0(out_name, ".png")),
            width = 1200, height = 800, res = 150); plotfun_th(); dev.off()
        pdf(file.path(out_dir, paste0(out_name, ".pdf")),
            width = 9, height = 6); plotfun_th(); dev.off()
        png(file.path(out_dir, paste0(out_name, "_gamma.png")),
            width = 1200, height = 800, res = 150); plotfun_g(); dev.off()
        pdf(file.path(out_dir, paste0(out_name, "_gamma.pdf")),
            width = 9, height = 6); plotfun_g(); dev.off()
      }
    }
    
    ## ======================= CASE 3: gamma + thresholds, stacked 2×1 =======================
  } else {  # gamma_layout == "stacked"
    plotfun <- function() {
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      layout(matrix(1:2, nrow = 2), heights = c(0.35, 0.65))  # (% gamma, % thresholds)
      
      # Top: gamma (shared x, no axis)
      gamma_panel(show_x_title = FALSE, show_x_axis_ticks = FALSE,
                  mar_override = c(0.4, left, top, right))
      
      # Bottom: thresholds (x-axis shown)
      thresholds_panel(show_x_title = TRUE,
                       mar_override = c(bottom, left, 0.4, right))
    }
    
    if (exists("render_and_save")) {
      render_and_save(fname_base = out_name, plotfun = plotfun,
                      show = show, save = save, dir = out_dir,
                      width_in = width_in, height_in = height_in, dpi = dpi,
                      match_current = match_current)
    } else {
      if (isTRUE(show)) plotfun()
      if (isTRUE(save)) {
        png(file.path(out_dir, paste0(out_name, ".png")),
            width = 1600, height = 1600, res = 150); plotfun(); dev.off()
        pdf(file.path(out_dir, paste0(out_name, ".pdf")),
            width = 12, height = 12); plotfun(); dev.off()
      }
    }
  }
  
  invisible(NULL)
}


# ======================= Cost of misspecification surfaces =======================
# For each (delta, eps), compute:
#   gamma_with_ambiguity:
#     optimal robust gamma under that ambiguous model
#   gamma_zero_ambiguity_policy:
#     gamma obtained by fixing the certainty-optimal reflecting barriers,
#     re-solving the ambiguity thresholds under the same ambiguous model,
#     and evaluating the resulting ergodic value
#   relative_cost:
#     unscaled ratio kept for internal diagnostics/backward compatibility
#   relative_cost_pct:
#     public-facing RMC, defined as 100 * relative_cost
.safe_relative_change <- function(num, den, tol = 1e-12) {
  out <- rep(NA_real_, length(num))
  ok <- is.finite(num) & is.finite(den) & (abs(den) > tol)
  out[ok] <- (num[ok] - den[ok]) / den[ok]
  out
}

# Public-facing reporting now defines RMC in percentage terms. Keep the
# unscaled relative-cost ratio internally, but silently map public reporting
# requests onto the percentage-scaled surface.
.misspecification_reporting_surface <- function(surface) {
  if (identical(surface, "relative_cost")) "relative_cost_pct" else surface
}

.ensure_misspecification_surface_has_finite_values <- function(grid_obj, surface) {
  zmat <- grid_obj[[surface]]
  if (is.null(zmat) || !is.matrix(zmat)) {
    stop(sprintf("Surface '%s' is not available in grid_obj.", surface))
  }
  if (any(is.finite(zmat))) {
    return(invisible(zmat))
  }
  
  msg <- sprintf("Surface '%s' does not contain any finite values.", surface)
  certainty_error_msg <- grid_obj$certainty_solution_error
  certainty_failed <- is.null(grid_obj$certainty_solution) &&
    is.character(certainty_error_msg) &&
    length(certainty_error_msg) >= 1L &&
    !is.na(certainty_error_msg[1]) &&
    nzchar(certainty_error_msg[1])
  
  if (certainty_failed && surface %in% c(
    "relative_cost_pct",
    "relative_cost",
    "gamma_difference",
    "gamma_zero_ambiguity_policy"
  )) {
    msg <- paste0(
      msg,
      " The certainty baseline could not be solved, so the RMC comparison surface is undefined. ",
      "certainty_solution_error: ", certainty_error_msg[1]
    )
    
    cp <- grid_obj$certainty_params
    if (!is.null(cp$b) && isTRUE(all.equal(as.numeric(cp$b), 0))) {
      msg <- paste0(
        msg,
        " At certainty (delta = 0, epsilon = 0), this means a_{i,1} = mu * b = 0, ",
        "which is exactly the unimplemented S1 cubic branch."
      )
    }
  }
  
  stop(msg)
}

.misspecification_rdata_path <- function(out_dir, out_name) {
  file.path(out_dir, paste0(out_name, ".RData"))
}

save_misspecification_cost_rdata <- function(
    grid_obj,
    out_dir = grid_obj$out_dir %||% file.path("figures", "misspecification"),
    out_name = NULL,
    object_name = "grid_obj"
) {
  if (is.null(out_name)) {
    out_name <- grid_obj$out_name %||% "misspecification"
  }
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_path <- .misspecification_rdata_path(out_dir, out_name)
  grid_obj$out_dir <- out_dir
  grid_obj$out_name <- out_name
  grid_obj$rdata_path <- out_path
  
  save_env <- new.env(parent = emptyenv())
  assign(object_name, grid_obj, envir = save_env)
  save(list = object_name, file = out_path, envir = save_env)
  
  message(sprintf("Saved misspecification grid to: %s", normalizePath(out_path)))
  invisible(out_path)
}

load_misspecification_cost_rdata <- function(
    path = NULL,
    out_dir = file.path("figures", "misspecification"),
    out_name = NULL,
    object_name = NULL
) {
  if (is.null(path)) {
    if (is.null(out_name)) {
      stop("Provide either 'path' or 'out_name' to load a misspecification grid.")
    }
    path <- .misspecification_rdata_path(out_dir, out_name)
  }
  
  if (!file.exists(path)) {
    stop(sprintf("Misspecification grid file not found: %s", path))
  }
  
  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = load_env)
  
  if (is.null(object_name)) {
    if (length(loaded_names) != 1L) {
      stop("The .RData file contains multiple objects; please specify 'object_name'.")
    }
    object_name <- loaded_names[[1]]
  }
  
  if (!exists(object_name, envir = load_env, inherits = FALSE)) {
    stop(sprintf("Object '%s' was not found in %s", object_name, path))
  }
  
  grid_obj <- get(object_name, envir = load_env, inherits = FALSE)
  if (!is.list(grid_obj)) {
    stop(sprintf("Object '%s' in %s is not a misspecification grid object.", object_name, path))
  }
  
  if (is.null(grid_obj$out_dir) || !nzchar(grid_obj$out_dir)) {
    grid_obj$out_dir <- dirname(path)
  }
  if (is.null(grid_obj$out_name) || !nzchar(grid_obj$out_name)) {
    grid_obj$out_name <- tools::file_path_sans_ext(basename(path))
  }
  grid_obj$rdata_path <- path
  class(grid_obj) <- unique(c("misspecification_cost_grid", class(grid_obj)))
  
  grid_obj
}

misspecification_cost_grid <- function(
    delta_values,
    eps_values,
    b = 0, r = 1, sigma = 1, mu = 1, u = 1.0, l = 1.0,
    certainty_delta = 0, certainty_eps = 0,
    xL0 = -0.5, xk0 = -0.1, xl0 = 0.3, xU0 = 1.0,
    method = "Broyden",
    control = list(ftol=1e-11, xtol=1e-11, maxit=1e6, stepmax=1, trace=0),
    tol_build = 1e-10,
    solver_verbose = FALSE,
    tol_regime_switch = 1e-3,
    switch_to_newton = FALSE,
    switch_ftol = 5e-2,
    switch_iter_min = 2L,
    newton_control = list(ftol=1e-12, xtol=1e-12, maxit=1e5, stepmax=1, trace=0),
    fd_step = 1e-6,
    relative_tol = 1e-12,
    verbose = TRUE,
    out_dir = file.path("figures", "misspecification"),
    out_name = NULL,
    save = FALSE
) {
  delta_values <- as.numeric(delta_values)
  eps_values <- as.numeric(eps_values)
  
  if (!length(delta_values) || !length(eps_values)) {
    stop("delta_values and eps_values must both be non-empty.")
  }
  if (any(!is.finite(delta_values)) || any(delta_values < 0)) {
    stop("delta_values must be finite and non-negative.")
  }
  if (any(!is.finite(eps_values)) || any(eps_values < 0 | eps_values > 1)) {
    stop("eps_values must be finite and lie in [0, 1].")
  }
  
  certainty_params <- make_params(
    b = b, delta = certainty_delta, r = r, eps = certainty_eps,
    sigma = sigma, mu = mu, u = u, l = l
  )
  
  certainty_sol_try <- tryCatch(
    solve_optimal_barriers(
      certainty_params, xL0, xk0, xl0, xU0,
      method = method,
      control = control,
      tol_build = tol_build,
      verbose = solver_verbose,
      tol_regime_switch = tol_regime_switch,
      switch_to_newton = switch_to_newton,
      switch_ftol = switch_ftol,
      switch_iter_min = switch_iter_min,
      newton_control = newton_control,
      fd_step = fd_step
    ),
    error = function(e) e
  )
  certainty_sol <- if (inherits(certainty_sol_try, "error")) NULL else certainty_sol_try
  certainty_error_msg <- if (inherits(certainty_sol_try, "error")) {
    conditionMessage(certainty_sol_try)
  } else {
    NA_character_
  }
  
  dim_labs <- function(x) format(signif(x, 8), trim = TRUE, scientific = FALSE)
  dims <- list(delta = dim_labs(delta_values), eps = dim_labs(eps_values))
  
  gamma_with_ambiguity <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                                 dimnames = dims)
  gamma_zero_ambiguity_policy <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                                        dimnames = dims)
  gamma_difference <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                             dimnames = dims)
  relative_cost <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                          dimnames = dims)
  relative_cost_pct <- matrix(NA_real_, nrow = length(delta_values), ncol = length(eps_values),
                              dimnames = dims)
  robust_converged <- matrix(FALSE, nrow = length(delta_values), ncol = length(eps_values),
                             dimnames = dims)
  zero_policy_converged <- matrix(FALSE, nrow = length(delta_values), ncol = length(eps_values),
                                  dimnames = dims)
  error_msg <- matrix("", nrow = length(delta_values), ncol = length(eps_values),
                      dimnames = dims)
  
  append_error <- function(existing, msg) {
    if (!nzchar(existing)) msg else paste(existing, msg, sep = " | ")
  }
  start_vals <- function(sol) {
    if (!is.null(sol) && !is.null(sol$x)) {
      return(sol$x[c("xL", "xk", "xl", "xU")])
    }
    c(xL = xL0, xk = xk0, xl = xl0, xU = xU0)
  }
  
  total <- length(delta_values) * length(eps_values)
  counter <- 0L
  
  for (i in seq_along(delta_values)) {
    last_robust_sol <- certainty_sol
    last_zero_sol <- certainty_sol
    
    for (j in seq_along(eps_values)) {
      counter <- counter + 1L
      delta_i <- delta_values[i]
      eps_j <- eps_values[j]
      
      if (isTRUE(verbose)) {
        message(sprintf("[misspecification] %3d/%d  delta = %.6g, eps = %.6g",
                        counter, total, delta_i, eps_j))
        flush.console()
      }
      
      amb_params <- make_params(
        b = b, delta = delta_i, r = r, eps = eps_j,
        sigma = sigma, mu = mu, u = u, l = l
      )
      
      same_as_certainty <-
        isTRUE(all.equal(delta_i, certainty_delta)) &&
        isTRUE(all.equal(eps_j, certainty_eps))
      robust_start <- start_vals(last_robust_sol)
      
      robust_try <- tryCatch(
        if (same_as_certainty && !is.null(certainty_sol)) {
          certainty_sol
        } else {
          solve_optimal_barriers(
            amb_params,
            robust_start[["xL"]],
            robust_start[["xk"]],
            robust_start[["xl"]],
            robust_start[["xU"]],
            method = method,
            control = control,
            tol_build = tol_build,
            verbose = solver_verbose,
            tol_regime_switch = tol_regime_switch,
            switch_to_newton = switch_to_newton,
            switch_ftol = switch_ftol,
            switch_iter_min = switch_iter_min,
            newton_control = newton_control,
            fd_step = fd_step
          )
        },
        error = function(e) e
      )
      
      if (inherits(robust_try, "error")) {
        error_msg[i, j] <- append_error(
          error_msg[i, j],
          paste0("robust policy: ", conditionMessage(robust_try))
        )
      } else {
        gamma_with_ambiguity[i, j] <- robust_try$gamma
        robust_converged[i, j] <- TRUE
        last_robust_sol <- robust_try
      }
      
      if (is.null(certainty_sol)) {
        zero_try <- structure(
          list(message = paste0(
            "certainty baseline unavailable: ", certainty_error_msg
          )),
          class = c("simpleError", "error", "condition")
        )
      } else {
        zero_start <- start_vals(last_zero_sol)
        zero_try <- tryCatch(
          if (same_as_certainty) {
            certainty_sol
          } else {
            solve_ambiguity_thresholds_fixed_barriers(
              amb_params,
              certainty_sol$x$xL,
              certainty_sol$x$xU,
              zero_start[["xk"]],
              zero_start[["xl"]],
              method = method,
              control = control,
              tol_build = tol_build,
              verbose = solver_verbose,
              tol_regime_switch = tol_regime_switch,
              switch_to_newton = switch_to_newton,
              switch_ftol = switch_ftol,
              switch_iter_min = switch_iter_min,
              newton_control = newton_control,
              fd_step = fd_step
            )
          },
          error = function(e) e
        )
      }
      
      if (inherits(zero_try, "error")) {
        error_msg[i, j] <- append_error(
          error_msg[i, j],
          paste0("certainty-barriers policy: ", conditionMessage(zero_try))
        )
      } else {
        gamma_zero_ambiguity_policy[i, j] <- zero_try$gamma
        zero_policy_converged[i, j] <- TRUE
        last_zero_sol <- zero_try
      }
      
      if (robust_converged[i, j] && zero_policy_converged[i, j]) {
        gamma_difference[i, j] <-
          gamma_zero_ambiguity_policy[i, j] - gamma_with_ambiguity[i, j]
        relative_cost[i, j] <-
          .safe_relative_change(gamma_zero_ambiguity_policy[i, j],
                                gamma_with_ambiguity[i, j],
                                tol = relative_tol)
        relative_cost_pct[i, j] <- 100 * relative_cost[i, j]
      }
    }
  }
  
  long_results <- data.frame(
    delta = rep(delta_values, times = length(eps_values)),
    eps = rep(eps_values, each = length(delta_values)),
    gamma_robust = as.vector(gamma_with_ambiguity),
    gamma_nonrobust_worst_case = as.vector(gamma_zero_ambiguity_policy),
    gamma_with_ambiguity = as.vector(gamma_with_ambiguity),
    gamma_zero_ambiguity_policy = as.vector(gamma_zero_ambiguity_policy),
    gamma_difference = as.vector(gamma_difference),
    relative_cost = as.vector(relative_cost),
    relative_cost_pct = as.vector(relative_cost_pct),
    robust_converged = as.vector(robust_converged),
    zero_policy_converged = as.vector(zero_policy_converged),
    converged = as.vector(robust_converged & zero_policy_converged),
    error_msg = as.vector(error_msg),
    stringsAsFactors = FALSE
  )
  long_results$error_msg[!nzchar(long_results$error_msg)] <- NA_character_
  
  if (is.null(out_name)) {
    out_name <- sprintf("misspecification_delta_%g_to_%g__eps_%g_to_%g",
                        min(delta_values), max(delta_values),
                        min(eps_values), max(eps_values))
  }
  
  out <- list(
    delta_values = delta_values,
    eps_values = eps_values,
    certainty_params = certainty_params,
    certainty_solution = certainty_sol,
    certainty_solution_error = certainty_error_msg,
    gamma_robust = gamma_with_ambiguity,
    gamma_nonrobust_worst_case = gamma_zero_ambiguity_policy,
    gamma_with_ambiguity = gamma_with_ambiguity,
    gamma_zero_ambiguity_policy = gamma_zero_ambiguity_policy,
    gamma_difference = gamma_difference,
    relative_cost = relative_cost,
    relative_cost_pct = relative_cost_pct,
    robust_converged = robust_converged,
    zero_policy_converged = zero_policy_converged,
    error_msg = error_msg,
    long_results = long_results,
    out_dir = out_dir,
    out_name = out_name,
    rdata_path = NULL
  )
  class(out) <- c("misspecification_cost_grid", class(out))
  
  if (save) {
    out$rdata_path <- save_misspecification_cost_rdata(
      out,
      out_dir = out_dir,
      out_name = out_name
    )
  }
  
  out
}

export_misspecification_cost_latex <- function(
    grid_obj,
    surface = c("relative_cost_pct", "relative_cost"),
    out_dir = file.path("figures", "misspecification"),
    out_name = NULL,
    save = TRUE,
    digits = NULL,
    format = "f",
    na_string = "--",
    include_table_env = TRUE,
    placement = "htbp",
    centering = TRUE,
    booktabs = TRUE,
    caption = NULL,
    label = NULL,
    row_header = "\\diagbox[width=1.6cm,height=0.7cm]{$\\delta$}{$\\epsilon$}",
    column_align = NULL
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  
  if (is.null(grid_obj[[surface]]) || !is.matrix(grid_obj[[surface]])) {
    stop(sprintf("Surface '%s' is not available in grid_obj.", surface))
  }
  
  zmat <- grid_obj[[surface]]
  if (is.null(digits)) {
    digits <- 2L
  }
  
  if (is.null(caption)) {
    caption <- "RMC matrix"
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_table")
  }
  
  if (is.null(label) && isTRUE(include_table_env)) {
    label <- paste0("tab:", gsub("[^A-Za-z0-9]+", "_", out_name))
  }
  
  row_labels <- rownames(zmat)
  if (is.null(row_labels)) {
    row_labels <- format(signif(seq_len(nrow(zmat)), 8), trim = TRUE, scientific = FALSE)
  }
  col_labels <- colnames(zmat)
  if (is.null(col_labels)) {
    col_labels <- format(signif(seq_len(ncol(zmat)), 8), trim = TRUE, scientific = FALSE)
  }
  
  if (is.null(column_align)) {
    column_align <- paste0("r|", paste(rep("r", ncol(zmat)), collapse = ""))
  }
  
  format_cell <- function(x) {
    out <- rep(na_string, length(x))
    ok <- is.finite(x)
    if (any(ok)) {
      xx <- x[ok]
      xx[abs(xx) < 0.5 * 10^(-digits)] <- 0
      out[ok] <- formatC(xx, format = format, digits = digits)
    }
    out
  }
  
  header_row <- paste0(
    paste(c(row_header, col_labels), collapse = " & "),
    " \\\\"
  )
  body_rows <- vapply(
    seq_len(nrow(zmat)),
    function(i) {
      paste0(
        paste(c(row_labels[i], format_cell(zmat[i, ])), collapse = " & "),
        " \\\\"
      )
    },
    character(1)
  )
  
  lines <- character(0)
  if (isTRUE(include_table_env)) {
    lines <- c(lines, paste0("\\begin{table}[", placement, "]"))
    if (isTRUE(centering)) lines <- c(lines, "\\centering")
    if (!is.null(caption)) lines <- c(lines, paste0("\\caption{", caption, "}"))
    if (!is.null(label))   lines <- c(lines, paste0("\\label{", label, "}"))
  }
  
  lines <- c(lines, paste0("\\begin{tabular}{", column_align, "}"))
  lines <- c(lines, if (isTRUE(booktabs)) "\\toprule" else "\\hline")
  lines <- c(lines, header_row)
  lines <- c(lines, if (isTRUE(booktabs)) "\\midrule" else "\\hline")
  lines <- c(lines, body_rows)
  lines <- c(lines, if (isTRUE(booktabs)) "\\bottomrule" else "\\hline")
  lines <- c(lines, "\\end{tabular}")
  
  if (isTRUE(include_table_env)) {
    lines <- c(lines, "\\end{table}")
  }
  
  latex <- paste(lines, collapse = "\n")
  
  out_path <- NULL
  if (isTRUE(save)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_path <- file.path(out_dir, paste0(out_name, ".tex"))
    writeLines(lines, out_path)
    message(sprintf("Saved LaTeX table to: %s", normalizePath(out_path)))
  }
  
  invisible(list(
    latex = latex,
    path = out_path,
    matrix = zmat,
    surface = surface
  ))
}

.validate_margin_spec <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 4L || any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("'%s' must be a numeric vector of length 4 with non-negative entries.", name))
  }
  x
}

.validate_fig_spec <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 4L || any(!is.finite(x))) {
    stop(sprintf("'%s' must be a numeric vector of length 4.", name))
  }
  if (any(x < 0 | x > 1)) {
    stop(sprintf("'%s' entries must lie between 0 and 1.", name))
  }
  if (x[1] >= x[2] || x[3] >= x[4]) {
    stop(sprintf("'%s' must satisfy left < right and bottom < top.", name))
  }
  x
}

.misspecification_surface_plot_meta <- function(surface) {
  list(
    z_axis_title = switch(
      surface,
      relative_cost = expression(plain(RMC) / 100 == (gamma[NR*","*W] - gamma[R]) / gamma[R]),
      relative_cost_pct = expression(plain(RMC)),
      gamma_difference = expression(gamma[NR*","*W] - gamma[R]),
      gamma_with_ambiguity = expression(gamma[R]),
      gamma_zero_ambiguity_policy = expression(gamma[NR*","*W])
    ),
    color_scale_title = switch(
      surface,
      relative_cost = expression(plain(RMC) / 100),
      relative_cost_pct = expression(plain(RMC)),
      gamma_difference = expression(gamma[NR*","*W] - gamma[R]),
      gamma_with_ambiguity = expression(gamma[R]),
      gamma_zero_ambiguity_policy = expression(gamma[NR*","*W])
    ),
    main = switch(
      surface,
      relative_cost = "RMC / 100",
      relative_cost_pct = "RMC",
      gamma_difference = "Excess ergodic cost of the non-robust policy",
      gamma_with_ambiguity = "Robust ergodic cost",
      gamma_zero_ambiguity_policy = "Non-robust policy under the worst-case model"
    )
  )
}

.validate_unit_interval_scalar <- function(x, name) {
  x <- as.numeric(x)[1]
  if (!is.finite(x) || x < 0 || x > 1) {
    stop(sprintf("'%s' must be a single number between 0 and 1.", name))
  }
  x
}

.validate_nonnegative_scalar <- function(x, name) {
  x <- as.numeric(x)[1]
  if (!is.finite(x) || x < 0) {
    stop(sprintf("'%s' must be a single non-negative number.", name))
  }
  x
}

.validate_xy_nudge <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 2L || any(!is.finite(x))) {
    stop(sprintf("'%s' must be a numeric vector of length 2.", name))
  }
  x
}

.coerce_plot_label <- function(label) {
  if (is.null(label) || !length(label)) return("")
  if (is.language(label) && !is.expression(label)) {
    return(as.expression(label))
  }
  label
}

.draw_parametrizer_legend_key <- function(
    label,
    corner = c("topleft", "bottomleft"),
    inset = c(0.04, 0.06),
    line_length = 0.18,
    cex = 1,
    col = "#333333",
    lty = 1,
    lwd = 1.5
) {
  label <- .coerce_plot_label(label)
  corner <- match.arg(corner)
  inset <- as.numeric(inset)
  if (length(inset) != 2L || any(!is.finite(inset)) ||
      any(inset < 0 | inset > 1)) {
    stop("'inset' must be a numeric vector of length 2 with entries in [0, 1].")
  }
  line_length <- .validate_unit_interval_scalar(line_length, "line_length")
  cex <- .validate_nonnegative_scalar(cex, "cex")
  
  usr <- par("usr")
  xspan <- diff(usr[1:2])
  yspan <- diff(usr[3:4])
  if (!is.finite(xspan) || !is.finite(yspan) || xspan <= 0 || yspan <= 0) {
    return(invisible(NULL))
  }
  
  x0 <- usr[1] + inset[1] * xspan
  x1 <- min(usr[2], x0 + line_length * xspan)
  y0 <- if (identical(corner, "topleft")) {
    usr[4] - inset[2] * yspan
  } else {
    usr[3] + inset[2] * yspan
  }
  
  total_width <- x1 - x0
  if (!is.finite(total_width) || total_width <= 0) {
    return(invisible(NULL))
  }
  
  gap_width <- max(
    1.25 * strwidth(label, cex = cex, units = "user"),
    0.025 * xspan
  )
  gap_width <- min(gap_width, 0.8 * total_width)
  seg_width <- 0.5 * (total_width - gap_width)
  
  if (seg_width > 0) {
    segments(x0, y0, x0 + seg_width, y0, col = col, lty = lty, lwd = lwd)
    segments(x1 - seg_width, y0, x1, y0, col = col, lty = lty, lwd = lwd)
  }
  
  text(x0 + total_width/2, y0, labels = label, cex = cex, col = col)
  
  invisible(NULL)
}

.persp_axis_edge_candidates <- function(axis, xlim, ylim, zlim) {
  axis <- match.arg(axis, c("x", "y", "z"))
  x0 <- xlim[1]; x1 <- xlim[2]
  y0 <- ylim[1]; y1 <- ylim[2]
  z0 <- zlim[1]; z1 <- zlim[2]
  
  switch(
    axis,
    x = list(
      rbind(c(x0, y0, z0), c(x1, y0, z0)),
      rbind(c(x0, y1, z0), c(x1, y1, z0)),
      rbind(c(x0, y0, z1), c(x1, y0, z1)),
      rbind(c(x0, y1, z1), c(x1, y1, z1))
    ),
    y = list(
      rbind(c(x0, y0, z0), c(x0, y1, z0)),
      rbind(c(x1, y0, z0), c(x1, y1, z0)),
      rbind(c(x0, y0, z1), c(x0, y1, z1)),
      rbind(c(x1, y0, z1), c(x1, y1, z1))
    ),
    z = list(
      rbind(c(x0, y0, z0), c(x0, y0, z1)),
      rbind(c(x1, y0, z0), c(x1, y0, z1)),
      rbind(c(x0, y1, z0), c(x0, y1, z1)),
      rbind(c(x1, y1, z0), c(x1, y1, z1))
    )
  )
}

.resolve_persp_axis_edge <- function(axis, pmat, xlim, ylim, zlim) {
  candidates <- .persp_axis_edge_candidates(axis, xlim, ylim, zlim)
  center_proj <- trans3d(mean(xlim), mean(ylim), mean(zlim), pmat = pmat)
  center_xy <- c(center_proj$x, center_proj$y)
  
  mids <- do.call(
    rbind,
    lapply(candidates, function(edge) {
      proj <- trans3d(edge[, 1], edge[, 2], edge[, 3], pmat = pmat)
      c(mean(proj$x), mean(proj$y))
    })
  )
  
  pick <- switch(
    axis,
    x = order(mids[, 2], -abs(mids[, 1] - center_xy[1]))[1],
    y = order(mids[, 2], -abs(mids[, 1] - center_xy[1]))[1],
    z = order(-abs(mids[, 1] - center_xy[1]), mids[, 2])[1]
  )
  
  candidates[[pick]]
}

.persp_axis_geometry <- function(axis, pmat, xlim, ylim, zlim) {
  edge <- .resolve_persp_axis_edge(axis, pmat = pmat, xlim = xlim, ylim = ylim, zlim = zlim)
  proj <- trans3d(edge[, 1], edge[, 2], edge[, 3], pmat = pmat)
  p0 <- c(proj$x[1], proj$y[1])
  p1 <- c(proj$x[2], proj$y[2])
  v <- p1 - p0
  len <- sqrt(sum(v^2))
  if (!is.finite(len) || len <= 0) {
    return(NULL)
  }
  
  dir <- v / len
  normal <- c(-dir[2], dir[1])
  center_proj <- trans3d(mean(xlim), mean(ylim), mean(zlim), pmat = pmat)
  center_xy <- c(center_proj$x, center_proj$y)
  mid <- 0.5 * (p0 + p1)
  if (sum(normal * (mid - center_xy)) < 0) {
    normal <- -normal
  }
  
  angle <- atan2(dir[2], dir[1]) * 180 / pi
  if (angle > 90) angle <- angle - 180
  if (angle < -90) angle <- angle + 180
  
  list(
    edge = edge,
    p0 = p0,
    p1 = p1,
    v = v,
    dir = dir,
    normal = normal,
    angle = angle
  )
}

.validate_positive_scalar_or_null <- function(x, name) {
  if (is.null(x)) return(NULL)
  x <- as.numeric(x)[1]
  if (!is.finite(x) || x <= 0) {
    stop(sprintf("'%s' must be NULL or a single positive number.", name))
  }
  x
}

.map_axis_values_to_persp <- function(values, from, stretch = 1, scale = TRUE) {
  dims <- dim(values)
  dim_names <- dimnames(values)
  values <- as.numeric(values)
  from <- as.numeric(from)
  span <- diff(from)
  if (!is.finite(span) || span <= 0) {
    out <- rep(0.5 * stretch, length(values))
  } else if (isTRUE(scale)) {
    out <- (values - from[1]) / span * stretch
  } else {
    out <- (values - from[1]) * stretch
  }
  
  if (!is.null(dims)) {
    dim(out) <- dims
    dimnames(out) <- dim_names
  }
  out
}

.pretty_ticks_within <- function(lim, n = 5) {
  lim <- sort(as.numeric(lim))
  at <- pretty(lim, n = n)
  tol <- sqrt(.Machine$double.eps) * max(1, diff(lim))
  at <- at[is.finite(at) & at >= (lim[1] - tol) & at <= (lim[2] + tol)]
  if (!length(at)) {
    at <- lim
  }
  sort(unique(at))
}

.format_projected_axis_labels <- function(values) {
  values <- as.numeric(values)
  tol <- sqrt(.Machine$double.eps) * max(1, max(abs(values), na.rm = TRUE))
  values[abs(values) < tol] <- 0
  format(signif(values, 8), trim = TRUE, scientific = FALSE)
}

.draw_manual_persp_axis_ticks <- function(
    axis,
    pmat,
    xlim,
    ylim,
    zlim,
    axis_limits,
    at = NULL,
    n = 5,
    cex = 1,
    col = par("col.axis"),
    font = par("font.axis"),
    lwd = 1,
    tcl = par("tcl"),
    ticktype = c("detailed", "simple"),
    axis_mgp = c(2.2, 0.7, 0),
    draw_axis_line = TRUE
) {
  ticktype <- match.arg(ticktype)
  geom <- .persp_axis_geometry(axis, pmat = pmat, xlim = xlim, ylim = ylim, zlim = zlim)
  if (is.null(geom)) {
    return(invisible(NULL))
  }
  
  axis_limits <- sort(as.numeric(axis_limits))
  plot_limits <- switch(axis, x = xlim, y = ylim, z = zlim)
  if (is.null(at)) {
    at <- .pretty_ticks_within(axis_limits, n = n)
  } else {
    at <- sort(unique(as.numeric(at)))
    at <- at[is.finite(at)]
  }
  tol <- sqrt(.Machine$double.eps) * max(1, diff(axis_limits))
  at <- at[at >= (axis_limits[1] - tol) & at <= (axis_limits[2] + tol)]
  if (!length(at)) {
    return(invisible(NULL))
  }
  
  at_plot <- .map_axis_values_to_persp(
    values = at,
    from = axis_limits,
    stretch = diff(plot_limits),
    scale = TRUE
  ) + plot_limits[1]
  
  coords <- switch(
    axis,
    x = cbind(at_plot, rep(geom$edge[1, 2], length(at_plot)), rep(geom$edge[1, 3], length(at_plot))),
    y = cbind(rep(geom$edge[1, 1], length(at_plot)), at_plot, rep(geom$edge[1, 3], length(at_plot))),
    z = cbind(rep(geom$edge[1, 1], length(at_plot)), rep(geom$edge[1, 2], length(at_plot)), at_plot)
  )
  proj_ticks <- trans3d(coords[, 1], coords[, 2], coords[, 3], pmat = pmat)
  tick_xy <- cbind(proj_ticks$x, proj_ticks$y)
  
  char_size <- max(
    strwidth("0", units = "user", cex = cex),
    strheight("0", units = "user", cex = cex)
  )
  tick_scale <- abs(as.numeric(tcl)[1])
  if (!is.finite(tick_scale)) tick_scale <- 0.25
  tick_length <- if (tick_scale == 0) 0 else {
    (0.35 + if (identical(ticktype, "detailed")) 0.12 else 0) *
      char_size * (tick_scale / 0.25)
  }
  label_gap <- tick_length + (0.35 + max(axis_mgp[2], 0)) * char_size
  
  if (isTRUE(draw_axis_line)) {
    segments(geom$p0[1], geom$p0[2], geom$p1[1], geom$p1[2], col = col, lwd = lwd, xpd = NA)
  }
  if (tick_length > 0) {
    tick_end <- tick_xy + matrix(geom$normal, nrow = nrow(tick_xy), ncol = 2, byrow = TRUE) * tick_length
    segments(tick_xy[, 1], tick_xy[, 2], tick_end[, 1], tick_end[, 2], col = col, lwd = lwd, xpd = NA)
  }
  
  label_xy <- tick_xy + matrix(geom$normal, nrow = nrow(tick_xy), ncol = 2, byrow = TRUE) * label_gap
  text(
    label_xy[, 1],
    label_xy[, 2],
    labels = .format_projected_axis_labels(at),
    cex = cex,
    col = col,
    font = font,
    srt = geom$angle,
    adj = c(0.5, 0.5),
    xpd = NA
  )
  
  invisible(at)
}

.draw_manual_persp_axis_title <- function(
    label,
    axis,
    pmat,
    xlim,
    ylim,
    zlim,
    at = 0.5,
    pad = 0.08,
    nudge = c(0, 0),
    cex = 1,
    col = par("col.lab"),
    font = par("font.lab")
) {
  label <- .coerce_plot_label(label)
  if (is.character(label) && (!length(label) || !nzchar(label[1]))) {
    return(invisible(NULL))
  }
  
  geom <- .persp_axis_geometry(axis, pmat = pmat, xlim = xlim, ylim = ylim, zlim = zlim)
  if (is.null(geom)) {
    return(invisible(NULL))
  }
  
  usr <- par("usr")
  pad_units <- pad * max(diff(usr[1:2]), diff(usr[3:4]))
  nudge_units <- c(
    nudge[1] * diff(usr[1:2]),
    nudge[2] * diff(usr[3:4])
  )
  pos <- geom$p0 + at * geom$v + pad_units * geom$normal + nudge_units
  
  text(
    pos[1],
    pos[2],
    labels = label,
    srt = geom$angle,
    adj = c(0.5, 0.5),
    xpd = NA,
    cex = cex,
    col = col,
    font = font
  )
  
  invisible(pos)
}

.resolve_color_scale_fig <- function(
    color_scale_position,
    color_scale_fig,
    color_scale_outer_fig
) {
  if (identical(color_scale_position, "figure")) {
    return(color_scale_fig)
  }
  
  omd <- par("omd")
  plt <- par("plt")
  if (!all(is.finite(omd)) || !all(is.finite(plt))) {
    stop("Could not determine the current plotting regions for the color scale.")
  }
  
  x_bounds <- c(omd[2], 1)
  if ((x_bounds[2] - x_bounds[1]) <= 0) {
    stop("Set outer_margins[4] > 0 to use color_scale_position = 'right_outer'.")
  }
  
  y_bounds <- c(
    omd[3] + (omd[4] - omd[3]) * plt[3],
    omd[3] + (omd[4] - omd[3]) * plt[4]
  )
  if ((y_bounds[2] - y_bounds[1]) <= 0) {
    stop("The plot region is too small to place the color scale in the right outer margin.")
  }
  
  c(
    x_bounds[1] + (x_bounds[2] - x_bounds[1]) * color_scale_outer_fig[1],
    x_bounds[1] + (x_bounds[2] - x_bounds[1]) * color_scale_outer_fig[2],
    y_bounds[1] + (y_bounds[2] - y_bounds[1]) * color_scale_outer_fig[3],
    y_bounds[1] + (y_bounds[2] - y_bounds[1]) * color_scale_outer_fig[4]
  )
}

.draw_vertical_color_scale_panel <- function(
    zlim,
    color_palette,
    color_scale_title,
    color_scale_axis_mgp,
    color_scale_nticks,
    color_scale_title_cex,
    color_scale_tick_cex,
    color_scale_border
) {
  op_scale <- par(
    mgp = color_scale_axis_mgp,
    xpd = NA,
    cex.axis = color_scale_tick_cex
  )
  on.exit(par(op_scale), add = TRUE)
  
  plot.new()
  plot.window(xlim = c(0, 1), ylim = zlim, xaxs = "i", yaxs = "i")
  
  y_breaks <- seq(zlim[1], zlim[2], length.out = length(color_palette) + 1L)
  for (k in seq_along(color_palette)) {
    rect(0, y_breaks[k], 1, y_breaks[k + 1L], col = color_palette[k], border = NA)
  }
  box(col = color_scale_border)
  tick_at <- pretty(zlim, n = color_scale_nticks)
  tick_tol <- sqrt(.Machine$double.eps) * max(1, diff(zlim))
  tick_at <- tick_at[tick_at >= (zlim[1] - tick_tol) & tick_at <= (zlim[2] + tick_tol)]
  if (!length(tick_at)) {
    tick_at <- zlim
  }
  axis(4, at = tick_at, las = 1)
  mtext(color_scale_title, side = 3, line = 0.15, cex = color_scale_title_cex)
  
  invisible(NULL)
}

.draw_vertical_color_scale <- function(
    zlim,
    color_palette,
    color_scale_title,
    legend_fig,
    color_scale_margins,
    color_scale_axis_mgp,
    color_scale_nticks,
    color_scale_title_cex,
    color_scale_tick_cex,
    color_scale_border
) {
  op_legend <- par(
    fig = legend_fig,
    mar = color_scale_margins,
    new = TRUE
  )
  on.exit(par(op_legend), add = TRUE)
  
  .draw_vertical_color_scale_panel(
    zlim = zlim,
    color_palette = color_palette,
    color_scale_title = color_scale_title,
    color_scale_axis_mgp = color_scale_axis_mgp,
    color_scale_nticks = color_scale_nticks,
    color_scale_title_cex = color_scale_title_cex,
    color_scale_tick_cex = color_scale_tick_cex,
    color_scale_border = color_scale_border
  )
}

plot_misspecification_cost_3d <- function(
    grid_obj,                                                   # output of misspecification_cost_grid()
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"), # which matrix in grid_obj to plot as the surface height
    out_dir = file.path("figures", "misspecification"),         # directory where the figure is saved
    out_name = NULL,                                            # base filename for saved outputs
    show = interactive(),                                       # if TRUE, draw on the current device
    save = FALSE,                                               # if TRUE, save PNG and PDF copies
    width_in = 6.67,                                            # saved figure width in inches
    height_in = 4.67,                                           # saved figure height in inches
    dpi = 300,                                                  # PNG resolution used on save
    match_current = FALSE,                                      # if TRUE, save using the current device size instead of width_in/height_in
    title = TRUE,                                               # if TRUE, use the default plot title
    main = NULL,                                                # custom main title; overrides the default when supplied
    show_x_axis_title = TRUE,                                   # if TRUE, show the x-axis title
    show_y_axis_title = TRUE,                                   # if TRUE, show the y-axis title
    show_z_axis_title = TRUE,                                   # if TRUE, show the z-axis title
    x_axis_title = expression(delta),                           # x-axis title shown along the projected x-axis
    y_axis_title = expression(epsilon),                         # y-axis title shown along the projected y-axis
    z_axis_title = NULL,                                        # z-axis title; if NULL, use the default for the chosen surface
    axis_title_cex = 1,                                         # size of the axis titles
    x_axis_title_at = 0.5,                                      # position of the x-axis title along its projected axis (0 = start, 1 = end)
    y_axis_title_at = 0.5,                                      # position of the y-axis title along its projected axis (0 = start, 1 = end)
    z_axis_title_at = 0.5,                                      # position of the z-axis title along its projected axis (0 = start, 1 = end)
    x_axis_title_pad = 0.08,                                    # outward distance of the x-axis title from its projected axis
    y_axis_title_pad = 0.08,                                    # outward distance of the y-axis title from its projected axis
    z_axis_title_pad = 0.10,                                    # outward distance of the z-axis title from its projected axis
    tick_cex = 1.2,                                             # size of the axis tick labels
    title_cex = 1,                                              # size of the main title
    margins = c(3.2, 3.2, 1.8, 0.8),                           # margins around the 3D panel, in par(mar) order: bottom, left, top, right
    axis_mgp = c(2.2, 0.7, 0),                                 # spacing control used by the custom projected tick labels
    tcl = -0.25,                                                # tick length
    theta = 35,                                                 # horizontal viewing angle
    phi = 25,                                                   # vertical viewing angle
    r = sqrt(3),                                                # distance of the eyepoint from the center of the 3D box
    d = 2,                                                      # strength of the perspective effect
    scale = TRUE,                                               # if TRUE, normalize x, y, and z before applying the internal stretch; if FALSE, preserve their original relative scales
    auto_stretch = TRUE,                                        # if TRUE, compute default internal x/y stretch factors from the actual panel aspect ratio
    x_stretch = NULL,                                           # manual stretch applied to the internal delta direction; NULL uses the automatic value when available
    y_stretch = NULL,                                           # manual stretch applied to the internal epsilon direction; NULL uses the automatic value when available
    expand = 0.7,                                               # vertical expansion factor of the surface
    shade = 0.35,                                               # amount of lighting/shading on the surface
    ltheta = -135,                                              # horizontal angle of the light source
    lphi = 0,                                                   # vertical angle of the light source
    use_z_colormap = TRUE,                                      # if TRUE, color the surface according to z-values
    color_min = "#1E90FF",                                      # color used for the minimum z-values
    color_max = "#B22222",                                      # color used for the maximum z-values
    color_palette = NULL,                                       # custom vector of colors; overrides color_min/color_max when supplied
    color_steps = 64,                                           # number of colors in the generated palette
    add_color_scale = TRUE,                                     # if TRUE, draw a vertical color scale in a narrow panel to the right
    color_scale_width = 0.14,                                   # fraction of the total figure width reserved for the right-hand color-scale panel
    color_scale_title = NULL,                                   # title shown above the color scale
    color_scale_nticks = 5,                                     # approximate number of tick marks on the color scale
    color_scale_title_cex = 0.9,                                # size of the color-scale title
    color_scale_tick_cex = 0.85,                                # size of the color-scale tick labels
    color_scale_border = "#444444",                             # border color of the color scale
    col = "#7FA7A2",                                            # fallback constant surface color when use_z_colormap = FALSE
    border = NA,                                                # border color of the surface facets
    box = TRUE,                                                 # if TRUE, draw the 3D bounding box
    axes = TRUE,                                                # if TRUE, draw custom projected axes and tick labels using the true, un-stretched values
    nticks = 5,                                                 # approximate number of tick marks per custom projected axis
    ticktype = c("detailed", "simple"),                         # style used for the custom projected tick marks
    zlim = NULL,                                                # custom z-range; if NULL, compute it from the plotted surface
    restore_par = TRUE                                          # if FALSE, do not restore par() on exit; useful inside a shared layout
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  ticktype <- match.arg(ticktype)
  surface_meta <- .misspecification_surface_plot_meta(surface)
  
  margins <- .validate_margin_spec(margins, "margins")
  x_axis_title_at <- .validate_unit_interval_scalar(x_axis_title_at, "x_axis_title_at")
  y_axis_title_at <- .validate_unit_interval_scalar(y_axis_title_at, "y_axis_title_at")
  z_axis_title_at <- .validate_unit_interval_scalar(z_axis_title_at, "z_axis_title_at")
  x_axis_title_pad <- .validate_nonnegative_scalar(x_axis_title_pad, "x_axis_title_pad")
  y_axis_title_pad <- .validate_nonnegative_scalar(y_axis_title_pad, "y_axis_title_pad")
  z_axis_title_pad <- .validate_nonnegative_scalar(z_axis_title_pad, "z_axis_title_pad")
  x_stretch <- .validate_positive_scalar_or_null(x_stretch, "x_stretch")
  y_stretch <- .validate_positive_scalar_or_null(y_stretch, "y_stretch")
  color_scale_width <- .validate_unit_interval_scalar(color_scale_width, "color_scale_width")
  if (color_scale_width <= 0 || color_scale_width >= 1) {
    stop("'color_scale_width' must be strictly between 0 and 1.")
  }
  color_steps <- as.integer(color_steps)[1]
  if (!is.finite(color_steps) || color_steps < 2L) {
    stop("'color_steps' must be an integer >= 2.")
  }
  
  if (is.null(grid_obj$delta_values) || is.null(grid_obj$eps_values)) {
    stop("grid_obj must be the output of misspecification_cost_grid().")
  }
  
  delta_values <- as.numeric(grid_obj$delta_values)
  eps_values <- as.numeric(grid_obj$eps_values)
  zmat <- grid_obj[[surface]]
  
  zmat <- .ensure_misspecification_surface_has_finite_values(grid_obj, surface)
  if (length(delta_values) < 2L || length(eps_values) < 2L) {
    stop("A 3D surface plot needs at least two delta values and two eps values.")
  }
  
  delta_ord <- order(delta_values)
  eps_ord <- order(eps_values)
  delta_values <- delta_values[delta_ord]
  eps_values <- eps_values[eps_ord]
  zmat <- zmat[delta_ord, eps_ord, drop = FALSE]
  zmat[!is.finite(zmat)] <- NA_real_
  
  if (is.null(zlim)) {
    zlim <- range(zmat, finite = TRUE)
  }
  if (diff(zlim) == 0) {
    bump <- max(1, abs(zlim[1]))
    zlim <- zlim + c(-1, 1) * 0.05 * bump
  }
  
  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(c(color_min, color_max))(color_steps)
  } else {
    color_palette <- as.character(color_palette)
  }
  if (!length(color_palette)) {
    stop("color_palette must contain at least one color.")
  }
  
  draw_color_scale <- isTRUE(use_z_colormap) && isTRUE(add_color_scale)
  
  if (is.null(z_axis_title)) {
    z_axis_title <- surface_meta$z_axis_title
  }
  if (is.null(color_scale_title)) {
    color_scale_title <- surface_meta$color_scale_title
  }
  if (is.null(main)) {
    main <- if (isTRUE(title)) surface_meta$main else ""
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_3d")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  if (isTRUE(use_z_colormap)) {
    n_delta <- nrow(zmat)
    n_eps <- ncol(zmat)
    z11 <- zmat[-n_delta, -n_eps, drop = FALSE]
    z21 <- zmat[-1,      -n_eps, drop = FALSE]
    z12 <- zmat[-n_delta, -1,     drop = FALSE]
    z22 <- zmat[-1,      -1,      drop = FALSE]
    facet_ok <- is.finite(z11) & is.finite(z21) & is.finite(z12) & is.finite(z22)
    facet_z <- matrix(NA_real_, nrow = n_delta - 1L, ncol = n_eps - 1L)
    facet_z[facet_ok] <- (z11[facet_ok] + z21[facet_ok] + z12[facet_ok] + z22[facet_ok]) / 4
    
    breaks <- seq(zlim[1], zlim[2], length.out = length(color_palette) + 1L)
    facet_cols <- matrix(grDevices::adjustcolor("white", alpha.f = 0),
                         nrow = n_delta - 1L, ncol = n_eps - 1L)
    bins <- cut(facet_z[facet_ok], breaks = breaks, include.lowest = TRUE, labels = FALSE)
    facet_cols[facet_ok] <- color_palette[bins]
  } else {
    facet_cols <- col
  }
  
  plotfun <- function() {
    if (isTRUE(restore_par)) {
      op <- par(c("mar", "mgp", "tcl", "cex.axis", "cex.lab", "cex.main"))
      on.exit(par(op), add = TRUE)
    }
    
    if (draw_color_scale) {
      layout(matrix(c(1, 2), nrow = 1), widths = c(1 - color_scale_width, color_scale_width))
      on.exit(layout(1), add = TRUE)
    }
    
    par(
      mar = margins,
      mgp = axis_mgp,
      tcl = tcl,
      cex.axis = tick_cex,
      cex.lab = axis_title_cex,
      cex.main = title_cex
    )
    
    panel_pin <- par("pin")
    panel_aspect <- if (all(is.finite(panel_pin)) && panel_pin[2] > 0) {
      panel_pin[1] / panel_pin[2]
    } else {
      1
    }
    auto_x_stretch <- max(1, panel_aspect)
    auto_y_stretch <- max(1, 1 / panel_aspect)
    x_stretch_use <- if (isTRUE(auto_stretch) && is.null(x_stretch)) auto_x_stretch else (x_stretch %||% 1)
    y_stretch_use <- if (isTRUE(auto_stretch) && is.null(y_stretch)) auto_y_stretch else (y_stretch %||% 1)
    
    delta_plot <- .map_axis_values_to_persp(
      values = delta_values,
      from = range(delta_values),
      stretch = x_stretch_use,
      scale = scale
    )
    eps_plot <- .map_axis_values_to_persp(
      values = eps_values,
      from = range(eps_values),
      stretch = y_stretch_use,
      scale = scale
    )
    z_plot <- .map_axis_values_to_persp(
      values = zmat,
      from = zlim,
      stretch = 1,
      scale = scale
    )
    zlim_plot <- .map_axis_values_to_persp(
      values = zlim,
      from = zlim,
      stretch = 1,
      scale = scale
    )
    
    pmat <- persp(
      x = delta_plot,
      y = eps_plot,
      z = z_plot,
      xlab = "",
      ylab = "",
      zlab = "",
      main = main,
      theta = theta,
      phi = phi,
      r = r,
      d = d,
      scale = FALSE,
      expand = expand,
      shade = shade,
      ltheta = ltheta,
      lphi = lphi,
      col = facet_cols,
      border = border,
      box = box,
      axes = FALSE,
      nticks = nticks,
      ticktype = ticktype,
      zlim = zlim_plot
    )
    
    if (isTRUE(axes)) {
      .draw_manual_persp_axis_ticks(
        axis = "x",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        axis_limits = range(delta_values),
        n = nticks,
        cex = tick_cex,
        lwd = 1,
        tcl = tcl,
        ticktype = ticktype,
        axis_mgp = axis_mgp,
        draw_axis_line = !isTRUE(box)
      )
      .draw_manual_persp_axis_ticks(
        axis = "y",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        axis_limits = range(eps_values),
        n = nticks,
        cex = tick_cex,
        lwd = 1,
        tcl = tcl,
        ticktype = ticktype,
        axis_mgp = axis_mgp,
        draw_axis_line = !isTRUE(box)
      )
      .draw_manual_persp_axis_ticks(
        axis = "z",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        axis_limits = zlim,
        n = nticks,
        cex = tick_cex,
        lwd = 1,
        tcl = tcl,
        ticktype = ticktype,
        axis_mgp = axis_mgp,
        draw_axis_line = !isTRUE(box)
      )
    }
    
    if (isTRUE(show_x_axis_title)) {
      .draw_manual_persp_axis_title(
        label = x_axis_title,
        axis = "x",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        at = x_axis_title_at,
        pad = x_axis_title_pad,
        cex = axis_title_cex
      )
    }
    if (isTRUE(show_y_axis_title)) {
      .draw_manual_persp_axis_title(
        label = y_axis_title,
        axis = "y",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        at = y_axis_title_at,
        pad = y_axis_title_pad,
        cex = axis_title_cex
      )
    }
    if (isTRUE(show_z_axis_title)) {
      .draw_manual_persp_axis_title(
        label = z_axis_title,
        axis = "z",
        pmat = pmat,
        xlim = range(delta_plot),
        ylim = range(eps_plot),
        zlim = zlim_plot,
        at = z_axis_title_at,
        pad = z_axis_title_pad,
        cex = axis_title_cex
      )
    }
    
    if (draw_color_scale) {
      par(
        mar = c(margins[1], 0.35, margins[3], 1.45),
        tcl = tcl,
        cex.lab = axis_title_cex,
        cex.main = title_cex
      )
      .draw_vertical_color_scale_panel(
        zlim = zlim,
        color_palette = color_palette,
        color_scale_title = color_scale_title,
        color_scale_axis_mgp = c(1.2, 0.25, 0),
        color_scale_nticks = color_scale_nticks,
        color_scale_title_cex = color_scale_title_cex,
        color_scale_tick_cex = color_scale_tick_cex,
        color_scale_border = color_scale_border
      )
    }
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_misspecification_cost_contour <- function(
    grid_obj,                                                   # output of misspecification_cost_grid()
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"), # which matrix in grid_obj to plot
    out_dir = file.path("figures", "misspecification"),         # directory where the figure is saved
    out_name = NULL,                                            # base filename for saved outputs
    show = interactive(),                                       # if TRUE, draw on the current device
    save = FALSE,                                               # if TRUE, save PNG and PDF copies
    width_in = 6.67,                                            # saved figure width in inches
    height_in = 4.67,                                           # saved figure height in inches
    dpi = 300,                                                  # PNG resolution used on save
    match_current = FALSE,                                      # if TRUE, save using the current device size instead
    title = TRUE,                                               # if TRUE, use the default plot title
    main = NULL,                                                # custom main title; overrides the default when supplied
    show_x_axis_title = TRUE,                                   # if TRUE, show the x-axis title
    show_y_axis_title = TRUE,                                   # if TRUE, show the y-axis title
    x_axis_title = expression(delta),                           # x-axis title
    y_axis_title = expression(epsilon),                         # y-axis title
    axis_title_cex = 1,                                         # size of the axis titles
    tick_cex = 1.2,                                             # size of the axis tick labels
    title_cex = 1,                                              # size of the main title
    base_cex = 1,                                               # global cex scaling inside the panel
    mex = 1,                                                    # margin expansion factor
    margins = c(3.2, 3.2, 1.8, 0.8),                           # par("mar")
    outer_margins = c(0, 0, 0, 0),                              # par("oma")
    axis_mgp = c(2.2, 0.7, 0),                                 # par("mgp")
    tcl = -0.25,                                                # tick length
    axes = TRUE,                                                # if TRUE, draw axes
    box = TRUE,                                                 # if TRUE, draw a box around the panel
    x_axis_at = NULL,                                           # custom x-axis tick positions
    y_axis_at = NULL,                                           # custom y-axis tick positions
    x_axis_labels = TRUE,                                       # x-axis labels (logical or custom vector)
    y_axis_labels = TRUE,                                       # y-axis labels (logical or custom vector)
    axis_las = 1,                                               # orientation of axis tick labels
    xlim = NULL,                                                # custom x-range
    ylim = NULL,                                                # custom y-range
    zlim = NULL,                                                # custom z-range for colors/contours
    asp = NA,                                                   # aspect ratio; NA leaves it unconstrained
    xaxs = "i",                                                 # x-axis style used by plot.window()
    yaxs = "i",                                                 # y-axis style used by plot.window()
    add_grid = FALSE,                                           # if TRUE, add a 2D reference grid
    grid_col = "#DDDDDD",                                       # grid color
    grid_lty = 3,                                               # grid line type
    grid_lwd = 0.8,                                             # grid line width
    use_z_colormap = FALSE,                                     # if TRUE, draw a filled color background; default leaves the panel white
    color_min = "#1E90FF",                                      # minimum-color anchor
    color_max = "#B22222",                                      # maximum-color anchor
    color_palette = NULL,                                       # custom vector of colors; overrides color_min/color_max
    color_steps = 64,                                           # number of colors in the generated palette
    useRaster = FALSE,                                          # passed to image()
    interpolate = FALSE,                                        # passed to image()
    add_contours = TRUE,                                        # if TRUE, overlay contour lines
    contour_levels = NULL,                                      # explicit contour levels; if NULL, use pretty()
    contour_nlevels = 10,                                       # target number of contour levels when contour_levels is NULL
    contour_drawlabels = TRUE,                                  # if TRUE, label contour lines
    contour_labcex = 0.8,                                       # contour-label size
    contour_method = c("flattest", "simple", "edge"),          # label-placement method for contour()
    contour_lwd = 1.3,                                          # contour-line width
    contour_lty = 1,                                            # contour-line type
    contour_col = "#333333",                                    # contour-line color
    contour_vfont = c("sans serif", "bold"),                    # contour-label vector font
    add_points = FALSE,                                         # if TRUE, show the evaluated (delta, eps) grid points
    point_pch = 16,                                             # point character for the grid points
    point_cex = 0.55,                                           # point size for the grid points
    point_col = "#111111",                                      # point color for the grid points
    add_parametrizer_legend = TRUE,                             # if TRUE and add_contours = TRUE, draw a small broken-line key for the contour level variable
    parametrizer_legend_label = NULL,                           # label shown in the contour-level key; NULL -> surface-specific default
    parametrizer_legend_inset = c(0.04, 0.08),                 # inset of the contour-level key from the lower-left corner
    parametrizer_legend_line_length = 0.18,                     # line length of the contour-level key as a fraction of panel width
    parametrizer_legend_cex = NULL,                             # text size of the contour-level key; NULL -> tick_cex
    parametrizer_legend_col = NULL,                             # color of the contour-level key; NULL -> contour_col
    parametrizer_legend_lty = NULL,                             # line type of the contour-level key; NULL -> contour_lty
    parametrizer_legend_lwd = NULL,                             # line width of the contour-level key; NULL -> contour_lwd
    add_color_scale = FALSE,                                    # if TRUE and use_z_colormap = TRUE, add the color-scale legend
    color_scale_position = c("figure", "right_outer"),          # legend placement mode
    color_scale_title = NULL,                                   # title shown above the color-scale legend
    color_scale_fig = c(0.86, 0.94, 0.26, 0.82),                # legend location in device coordinates when color_scale_position = "figure"
    color_scale_outer_fig = c(0.18, 0.82, 0.06, 0.94),          # legend location within the right outer-margin strip when color_scale_position = "right_outer"
    color_scale_outer_margin = 4.5,                             # minimum right outer margin reserved automatically for "right_outer"
    color_scale_margins = c(0.4, 0.1, 0.8, 1.8),                # margins used inside the legend panel
    color_scale_axis_mgp = c(1.2, 0.25, 0),                    # placement of ticks/title inside the legend panel
    color_scale_nticks = 5,                                     # number of tick marks on the color-scale legend
    color_scale_title_cex = 0.9,                                # size of the color-scale title
    color_scale_tick_cex = 0.85,                                # size of the color-scale tick labels
    color_scale_border = "#444444",                             # border color of the color-scale legend
    restore_par = TRUE                                          # if FALSE, do not restore par() on exit; useful inside a shared layout
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  contour_method <- match.arg(contour_method)
  color_scale_position <- match.arg(color_scale_position)
  surface_meta <- .misspecification_surface_plot_meta(surface)
  
  margins <- .validate_margin_spec(margins, "margins")
  outer_margins <- .validate_margin_spec(outer_margins, "outer_margins")
  color_scale_outer_margin <- as.numeric(color_scale_outer_margin)[1]
  if (!is.finite(color_scale_outer_margin) || color_scale_outer_margin < 0) {
    stop("'color_scale_outer_margin' must be a non-negative number.")
  }
  if (isTRUE(use_z_colormap) && isTRUE(add_color_scale) &&
      identical(color_scale_position, "right_outer")) {
    outer_margins[4] <- max(outer_margins[4], color_scale_outer_margin)
  }
  if (identical(color_scale_position, "figure")) {
    color_scale_fig <- .validate_fig_spec(color_scale_fig, "color_scale_fig")
  } else {
    color_scale_outer_fig <- .validate_fig_spec(color_scale_outer_fig, "color_scale_outer_fig")
  }
  
  color_steps <- as.integer(color_steps)[1]
  if (!is.finite(color_steps) || color_steps < 2L) {
    stop("'color_steps' must be an integer >= 2.")
  }
  contour_nlevels <- as.integer(contour_nlevels)[1]
  if (!is.finite(contour_nlevels) || contour_nlevels < 1L) {
    stop("'contour_nlevels' must be an integer >= 1.")
  }
  
  validate_range <- function(x, name) {
    x <- as.numeric(x)
    if (length(x) != 2L || any(!is.finite(x))) {
      stop(sprintf("'%s' must be a numeric vector of length 2.", name))
    }
    x <- sort(x)
    if (x[1] == x[2]) {
      stop(sprintf("'%s' must have distinct endpoints.", name))
    }
    x
  }
  
  if (is.null(grid_obj$delta_values) || is.null(grid_obj$eps_values)) {
    stop("grid_obj must be the output of misspecification_cost_grid().")
  }
  
  delta_values <- as.numeric(grid_obj$delta_values)
  eps_values <- as.numeric(grid_obj$eps_values)
  zmat <- grid_obj[[surface]]
  
  zmat <- .ensure_misspecification_surface_has_finite_values(grid_obj, surface)
  if (length(delta_values) < 2L || length(eps_values) < 2L) {
    stop("A contour plot needs at least two delta values and two eps values.")
  }
  if (!isTRUE(use_z_colormap) && !isTRUE(add_contours) && !isTRUE(add_points)) {
    stop("At least one of use_z_colormap, add_contours, or add_points must be TRUE.")
  }
  
  delta_ord <- order(delta_values)
  eps_ord <- order(eps_values)
  delta_values <- delta_values[delta_ord]
  eps_values <- eps_values[eps_ord]
  zmat <- zmat[delta_ord, eps_ord, drop = FALSE]
  zmat[!is.finite(zmat)] <- NA_real_
  
  if (is.null(xlim)) {
    xlim <- range(delta_values, finite = TRUE)
  } else {
    xlim <- validate_range(xlim, "xlim")
  }
  if (is.null(ylim)) {
    ylim <- range(eps_values, finite = TRUE)
  } else {
    ylim <- validate_range(ylim, "ylim")
  }
  
  if (is.null(zlim)) {
    zlim <- range(zmat, finite = TRUE)
  } else {
    zlim <- validate_range(zlim, "zlim")
  }
  if (diff(zlim) == 0) {
    bump <- max(1, abs(zlim[1]))
    zlim <- zlim + c(-1, 1) * 0.05 * bump
  }
  
  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(c(color_min, color_max))(color_steps)
  } else {
    color_palette <- as.character(color_palette)
  }
  if (!length(color_palette)) {
    stop("color_palette must contain at least one color.")
  }
  
  if (is.null(contour_levels)) {
    contour_levels <- pretty(zlim, n = contour_nlevels)
  }
  contour_levels <- sort(unique(as.numeric(contour_levels)))
  contour_levels <- contour_levels[is.finite(contour_levels)]
  contour_levels <- contour_levels[
    contour_levels >= zlim[1] & contour_levels <= zlim[2]
  ]
  if (isTRUE(add_contours) && !length(contour_levels)) {
    stop("No contour levels fall inside 'zlim'.")
  }
  
  if (is.null(color_scale_title)) {
    color_scale_title <- surface_meta$color_scale_title
  }
  if (is.null(parametrizer_legend_label)) {
    parametrizer_legend_label <- surface_meta$color_scale_title
  }
  if (is.null(parametrizer_legend_cex)) {
    parametrizer_legend_cex <- tick_cex
  }
  if (is.null(parametrizer_legend_col)) {
    parametrizer_legend_col <- contour_col
  }
  if (is.null(parametrizer_legend_lty)) {
    parametrizer_legend_lty <- contour_lty
  }
  if (is.null(parametrizer_legend_lwd)) {
    parametrizer_legend_lwd <- contour_lwd
  }
  if (is.null(main)) {
    main <- if (isTRUE(title)) surface_meta$main else ""
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_contour")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  plotfun <- function() {
    if (isTRUE(restore_par)) {
      op <- par(c("oma", "mar", "mgp", "tcl", "cex", "mex",
                  "cex.axis", "cex.lab", "cex.main"))
      on.exit(par(op), add = TRUE)
    }
    
    par_args <- list(
      mar = margins,
      mgp = axis_mgp,
      tcl = tcl,
      cex = base_cex,
      mex = mex,
      cex.axis = tick_cex,
      cex.lab = axis_title_cex,
      cex.main = title_cex
    )
    if (isTRUE(restore_par)) {
      par_args$oma <- outer_margins
    }
    do.call(par, par_args)
    
    plot.new()
    window_args <- list(xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
    if (!is.na(asp)) window_args$asp <- asp
    do.call(plot.window, window_args)
    
    if (isTRUE(use_z_colormap)) {
      image(
        x = delta_values,
        y = eps_values,
        z = zmat,
        zlim = zlim,
        col = color_palette,
        add = TRUE,
        xaxs = "i",
        yaxs = "i",
        useRaster = useRaster,
        interpolate = interpolate
      )
    }
    
    if (isTRUE(add_grid)) {
      grid(col = grid_col, lty = grid_lty, lwd = grid_lwd)
    }
    
    if (isTRUE(add_contours)) {
      contour(
        x = delta_values,
        y = eps_values,
        z = zmat,
        levels = contour_levels,
        add = TRUE,
        drawlabels = contour_drawlabels,
        method = contour_method,
        labcex = contour_labcex,
        lwd = contour_lwd,
        lty = contour_lty,
        col = contour_col,
        vfont = contour_vfont,
        axes = FALSE,
        frame.plot = FALSE
      )
    }
    
    if (isTRUE(add_points)) {
      pts <- expand.grid(delta = delta_values, eps = eps_values)
      points(pts$delta, pts$eps, pch = point_pch, cex = point_cex, col = point_col)
    }
    
    if (isTRUE(axes)) {
      axis(1, at = x_axis_at, labels = x_axis_labels, las = axis_las)
      axis(2, at = y_axis_at, labels = y_axis_labels, las = axis_las)
    }
    if (isTRUE(box)) {
      box()
    }
    
    title(
      main = main,
      xlab = if (show_x_axis_title) x_axis_title else "",
      ylab = if (show_y_axis_title) y_axis_title else ""
    )
    
    if (isTRUE(add_parametrizer_legend) && isTRUE(add_contours)) {
      .draw_parametrizer_legend_key(
        label = parametrizer_legend_label,
        corner = "bottomleft",
        inset = parametrizer_legend_inset,
        line_length = parametrizer_legend_line_length,
        cex = parametrizer_legend_cex,
        col = parametrizer_legend_col,
        lty = parametrizer_legend_lty,
        lwd = parametrizer_legend_lwd
      )
    }
    
    if (isTRUE(use_z_colormap) && isTRUE(add_color_scale)) {
      legend_fig <- .resolve_color_scale_fig(
        color_scale_position = color_scale_position,
        color_scale_fig = color_scale_fig,
        color_scale_outer_fig = color_scale_outer_fig
      )
      .draw_vertical_color_scale(
        zlim = zlim,
        color_palette = color_palette,
        color_scale_title = color_scale_title,
        legend_fig = legend_fig,
        color_scale_margins = color_scale_margins,
        color_scale_axis_mgp = color_scale_axis_mgp,
        color_scale_nticks = color_scale_nticks,
        color_scale_title_cex = color_scale_title_cex,
        color_scale_tick_cex = color_scale_tick_cex,
        color_scale_border = color_scale_border
      )
    }
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_misspecification_cost_vs_delta <- function(
    grid_obj,                                                   # output of misspecification_cost_grid()
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"), # which matrix in grid_obj to plot on the y-axis
    eps_values = NULL,                                          # optional subset of epsilon values to display
    out_dir = file.path("figures", "misspecification"),         # directory where the figure is saved
    out_name = NULL,                                            # base filename for saved outputs
    show = interactive(),                                       # if TRUE, draw on the current device
    save = FALSE,                                               # if TRUE, save PNG and PDF copies
    width_in = 6.67,                                            # saved figure width in inches
    height_in = 4.67,                                           # saved figure height in inches
    dpi = 300,                                                  # PNG resolution used on save
    match_current = FALSE,                                      # if TRUE, save using the current device size instead
    title = TRUE,                                               # if TRUE, use the default plot title
    main = NULL,                                                # custom main title; overrides the default when supplied
    show_x_axis_title = TRUE,                                   # if TRUE, show the x-axis title
    show_y_axis_title = TRUE,                                   # if TRUE, show the y-axis title
    x_axis_title = expression(delta),                           # x-axis title
    y_axis_title = NULL,                                        # y-axis title; if NULL, use the default for the chosen surface
    axis_title_cex = 1,                                         # size of the axis titles
    tick_cex = 1.2,                                             # size of the axis tick labels
    title_cex = 1,                                              # size of the main title
    base_cex = 1,                                               # global cex scaling inside the panel
    mex = 1,                                                    # margin expansion factor
    margins = c(3.2, 3.6, 1.8, 0.8),                           # par("mar")
    outer_margins = c(0, 0, 0, 0),                              # par("oma")
    axis_mgp = c(2.2, 0.7, 0),                                 # par("mgp")
    tcl = -0.25,                                                # tick length
    axes = TRUE,                                                # if TRUE, draw axes
    box = TRUE,                                                 # if TRUE, draw a box around the panel
    x_axis_at = NULL,                                           # custom x-axis tick positions
    y_axis_at = NULL,                                           # custom y-axis tick positions
    x_axis_labels = TRUE,                                       # x-axis labels (logical or custom vector)
    y_axis_labels = TRUE,                                       # y-axis labels (logical or custom vector)
    axis_las = 1,                                               # orientation of axis tick labels
    xlim = NULL,                                                # custom x-range
    ylim = NULL,                                                # custom y-range
    xaxs = "i",                                                 # x-axis style used by plot.window()
    yaxs = "i",                                                 # y-axis style used by plot.window()
    add_grid = FALSE,                                           # if TRUE, add a 2D reference grid
    grid_col = "#DDDDDD",                                       # grid color
    grid_lty = 3,                                               # grid line type
    grid_lwd = 0.8,                                             # grid line width
    type = c("l", "p", "b", "o"),                              # line/point style used for the epsilon curves
    cols = NULL,                                                # line colors; generated automatically when NULL
    ltys = 1,                                                   # line types for the epsilon curves
    lwds = 2,                                                   # line widths for the epsilon curves
    add_points = FALSE,                                         # if TRUE, add points even when type = "l"
    point_pch = 16,                                             # point character for the epsilon curves
    point_cex = 0.7,                                            # point size for the epsilon curves
    point_cols = NULL,                                          # point colors; defaults to cols
    label_eps_on_curve = TRUE,                                  # if TRUE, annotate each curve with epsilon directly on the line
    label_x = 0.9,                                              # target delta value where epsilon labels are placed
    label_gap = NULL,                                           # horizontal gap left in the line around the label; NULL -> automatic
    label_cex = 0.85,                                           # size of the inline epsilon labels
    label_col = NULL,                                           # label colors; defaults to the curve colors
    add_parametrizer_legend = TRUE,                             # if TRUE, draw a small broken-line key for the parametrizing variable
    add_legend = FALSE,                                         # if TRUE, show the epsilon legend
    legend_title = expression(epsilon),                         # legend title
    legend_labels = NULL,                                       # custom legend labels; if NULL, use the epsilon values
    legend_position = "topleft",                                # legend position keyword or c(x, y)
    legend_inset = 0.01,                                        # legend inset when using a keyword position
    legend_ncol = 1,                                            # number of legend columns
    legend_cex = 0.9,                                           # legend text size
    legend_bty = "n",                                           # legend box type
    legend_seg_len = 2.5,                                       # legend line-segment length
    restore_par = TRUE                                          # if FALSE, do not restore par() on exit; useful inside a shared layout
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  type <- match.arg(type)
  surface_meta <- .misspecification_surface_plot_meta(surface)
  
  margins <- .validate_margin_spec(margins, "margins")
  outer_margins <- .validate_margin_spec(outer_margins, "outer_margins")
  
  validate_range <- function(x, name) {
    x <- as.numeric(x)
    if (length(x) != 2L || any(!is.finite(x))) {
      stop(sprintf("'%s' must be a numeric vector of length 2.", name))
    }
    x <- sort(x)
    if (x[1] == x[2]) {
      stop(sprintf("'%s' must have distinct endpoints.", name))
    }
    x
  }
  
  if (is.null(grid_obj$delta_values) || is.null(grid_obj$eps_values)) {
    stop("grid_obj must be the output of misspecification_cost_grid().")
  }
  
  delta_values <- as.numeric(grid_obj$delta_values)
  eps_grid <- as.numeric(grid_obj$eps_values)
  zmat <- grid_obj[[surface]]
  
  zmat <- .ensure_misspecification_surface_has_finite_values(grid_obj, surface)
  if (!length(delta_values) || !length(eps_grid)) {
    stop("The selected grid must contain at least one delta value and one eps value.")
  }
  
  delta_ord <- order(delta_values)
  eps_ord <- order(eps_grid)
  delta_values <- delta_values[delta_ord]
  eps_grid <- eps_grid[eps_ord]
  zmat <- zmat[delta_ord, eps_ord, drop = FALSE]
  zmat[!is.finite(zmat)] <- NA_real_
  
  if (is.null(eps_values)) {
    eps_plot <- eps_grid
  } else {
    eps_values <- as.numeric(eps_values)
    if (!length(eps_values) || any(!is.finite(eps_values)) ||
        any(eps_values < 0 | eps_values > 1)) {
      stop("'eps_values' must be a non-empty numeric vector with entries in [0, 1].")
    }
    
    tol <- sqrt(.Machine$double.eps)
    idx <- vapply(
      eps_values,
      function(v) {
        hits <- which(abs(eps_grid - v) <= tol * max(1, abs(v)))
        if (length(hits)) hits[1] else NA_integer_
      },
      integer(1)
    )
    if (anyNA(idx)) {
      missing_eps <- format(signif(eps_values[is.na(idx)], 8), trim = TRUE, scientific = FALSE)
      stop(sprintf("The following eps_values were not found in grid_obj: %s",
                   paste(missing_eps, collapse = ", ")))
    }
    
    eps_plot <- eps_grid[idx]
    zmat <- zmat[, idx, drop = FALSE]
  }
  
  keep_curve <- apply(zmat, 2, function(y) any(is.finite(y)))
  if (!all(keep_curve)) {
    warning(sprintf("Dropping %d epsilon curve(s) with no finite values.", sum(!keep_curve)))
    eps_plot <- eps_plot[keep_curve]
    zmat <- zmat[, keep_curve, drop = FALSE]
  }
  if (!length(eps_plot) || !any(is.finite(zmat))) {
    stop("The selected surface does not contain any finite values to plot.")
  }
  
  if (is.null(xlim)) {
    xlim <- range(delta_values, finite = TRUE)
  } else {
    xlim <- validate_range(xlim, "xlim")
  }
  if (is.null(ylim)) {
    ylim <- range(zmat, finite = TRUE)
  } else {
    ylim <- validate_range(ylim, "ylim")
  }
  if (diff(ylim) == 0) {
    bump <- max(1, abs(ylim[1]))
    ylim <- ylim + c(-1, 1) * 0.05 * bump
  }
  
  label_x <- as.numeric(label_x)[1]
  if (!is.finite(label_x)) {
    stop("'label_x' must be a single finite number.")
  }
  if (!is.null(label_gap)) {
    label_gap <- as.numeric(label_gap)[1]
    if (!is.finite(label_gap) || label_gap < 0) {
      stop("'label_gap' must be NULL or a single non-negative number.")
    }
  }
  
  n_curves <- ncol(zmat)
  if (is.null(cols)) {
    cols <- tryCatch(
      grDevices::hcl.colors(n_curves, palette = "Dark 3"),
      error = function(e) grDevices::rainbow(n_curves, end = 0.85)
    )
  }
  cols <- rep_len(cols, n_curves)
  ltys <- rep_len(ltys, n_curves)
  lwds <- rep_len(lwds, n_curves)
  point_pch <- rep_len(point_pch, n_curves)
  point_cex <- rep_len(point_cex, n_curves)
  point_cols <- rep_len(point_cols %||% cols, n_curves)
  label_cols <- rep_len(label_col %||% cols, n_curves)
  
  if (is.null(legend_labels)) {
    legend_labels <- format(signif(eps_plot, 8), trim = TRUE, scientific = FALSE)
  } else if (length(legend_labels) != n_curves) {
    stop("'legend_labels' must have the same length as the number of plotted epsilon curves.")
  }
  
  if (is.null(y_axis_title)) {
    y_axis_title <- surface_meta$color_scale_title
  }
  if (is.null(main)) {
    main <- if (isTRUE(title)) paste(surface_meta$main, "vs delta") else ""
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_vs_delta")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  show_lines <- type %in% c("l", "b", "o")
  show_points <- isTRUE(add_points) || type %in% c("p", "b", "o")
  
  draw_curve_with_gap <- function(x, y, col, lty, lwd, x_gap_left, x_gap_right) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    
    if (!length(x)) return(invisible(NULL))
    if (length(x) == 1L) {
      return(invisible(NULL))
    }
    
    if (!is.finite(x_gap_left) || !is.finite(x_gap_right) || x_gap_left >= x_gap_right ||
        x_gap_right <= min(x) || x_gap_left >= max(x)) {
      lines(x, y, col = col, lty = lty, lwd = lwd)
      return(invisible(NULL))
    }
    
    x_gap_left <- max(x_gap_left, min(x))
    x_gap_right <- min(x_gap_right, max(x))
    if (x_gap_left >= x_gap_right) {
      lines(x, y, col = col, lty = lty, lwd = lwd)
      return(invisible(NULL))
    }
    
    y_gap_left <- approx(x, y, xout = x_gap_left, ties = mean, rule = 1)$y
    y_gap_right <- approx(x, y, xout = x_gap_right, ties = mean, rule = 1)$y
    
    x_left <- x[x < x_gap_left]
    y_left <- y[x < x_gap_left]
    if (is.finite(y_gap_left) && x_gap_left > min(x)) {
      x_left <- c(x_left, x_gap_left)
      y_left <- c(y_left, y_gap_left)
    }
    
    x_right <- x[x > x_gap_right]
    y_right <- y[x > x_gap_right]
    if (is.finite(y_gap_right) && x_gap_right < max(x)) {
      x_right <- c(x_gap_right, x_right)
      y_right <- c(y_gap_right, y_right)
    }
    
    if (length(x_left) >= 2L) {
      lines(x_left, y_left, col = col, lty = lty, lwd = lwd)
    }
    if (length(x_right) >= 2L) {
      lines(x_right, y_right, col = col, lty = lty, lwd = lwd)
    }
    
    invisible(NULL)
  }
  
  curve_label_info <- function(x, y, target_x, gap_width) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    
    if (length(x) < 2L) return(NULL)
    if (target_x < min(x) || target_x > max(x)) return(NULL)
    
    y_target <- approx(x, y, xout = target_x, ties = mean, rule = 1)$y
    if (!is.finite(y_target)) return(NULL)
    
    left_x <- max(min(x), target_x - gap_width/2)
    right_x <- min(max(x), target_x + gap_width/2)
    if (left_x >= right_x) return(NULL)
    
    list(x = target_x, y = y_target, left_x = left_x, right_x = right_x)
  }
  
  plotfun <- function() {
    if (isTRUE(restore_par)) {
      op <- par(c("oma", "mar", "mgp", "tcl", "cex", "mex",
                  "cex.axis", "cex.lab", "cex.main"))
      on.exit(par(op), add = TRUE)
    }
    
    par_args <- list(
      mar = margins,
      mgp = axis_mgp,
      tcl = tcl,
      cex = base_cex,
      mex = mex,
      cex.axis = tick_cex,
      cex.lab = axis_title_cex,
      cex.main = title_cex
    )
    if (isTRUE(restore_par)) {
      par_args$oma <- outer_margins
    }
    do.call(par, par_args)
    
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
    
    if (isTRUE(add_grid)) {
      grid(col = grid_col, lty = grid_lty, lwd = grid_lwd)
    }
    
    for (j in seq_len(n_curves)) {
      yj <- zmat[, j]
      okj <- is.finite(yj)
      if (!any(okj)) next
      
      xj <- delta_values[okj]
      yj <- yj[okj]
      label_info <- NULL
      if (isTRUE(label_eps_on_curve)) {
        gap_width_j <- label_gap
        if (is.null(gap_width_j)) {
          gap_width_j <- max(
            1.15 * strwidth(legend_labels[j], cex = label_cex, units = "user"),
            0.015 * diff(xlim)
          )
        }
        label_info <- curve_label_info(
          x = xj,
          y = yj,
          target_x = label_x,
          gap_width = gap_width_j
        )
      }
      
      if (show_lines) {
        if (is.null(label_info)) {
          lines(xj, yj, col = cols[j], lty = ltys[j], lwd = lwds[j])
        } else {
          draw_curve_with_gap(
            x = xj,
            y = yj,
            col = cols[j],
            lty = ltys[j],
            lwd = lwds[j],
            x_gap_left = label_info$left_x,
            x_gap_right = label_info$right_x
          )
        }
      }
      if (show_points) {
        if (!is.null(label_info)) {
          keep_pts <- (xj < label_info$left_x) | (xj > label_info$right_x)
          x_pts <- xj[keep_pts]
          y_pts <- yj[keep_pts]
        } else {
          x_pts <- xj
          y_pts <- yj
        }
        if (length(x_pts)) {
          points(x_pts, y_pts, pch = point_pch[j], cex = point_cex[j], col = point_cols[j])
        }
      }
      if (!is.null(label_info)) {
        text(
          x = label_info$x,
          y = label_info$y,
          labels = legend_labels[j],
          cex = label_cex,
          col = label_cols[j],
          xpd = NA
        )
      }
    }
    
    if (isTRUE(axes)) {
      axis(1, at = x_axis_at, labels = x_axis_labels, las = axis_las)
      axis(2, at = y_axis_at, labels = y_axis_labels, las = axis_las)
    }
    if (isTRUE(box)) {
      box()
    }
    
    title(
      main = main,
      xlab = if (show_x_axis_title) x_axis_title else "",
      ylab = if (show_y_axis_title) y_axis_title else ""
    )
    
    if (isTRUE(add_parametrizer_legend)) {
      .draw_parametrizer_legend_key(
        label = expression(epsilon),
        cex = tick_cex
      )
    }
    
    if (isTRUE(add_legend)) {
      legend_args <- list(
        legend = legend_labels,
        title = legend_title,
        col = cols,
        ncol = legend_ncol,
        cex = legend_cex,
        bty = legend_bty,
        seg.len = legend_seg_len,
        merge = TRUE
      )
      
      if (show_lines) {
        legend_args$lty <- ltys
        legend_args$lwd <- lwds
      }
      if (show_points) {
        legend_args$pch <- point_pch
        legend_args$pt.cex <- point_cex
      }
      
      if (is.character(legend_position) && length(legend_position) == 1L) {
        legend_args$x <- legend_position
        legend_args$inset <- legend_inset
      } else if (is.numeric(legend_position) && length(legend_position) == 2L &&
                 all(is.finite(legend_position))) {
        legend_args$x <- legend_position[1]
        legend_args$y <- legend_position[2]
      } else {
        stop("'legend_position' must be a single legend keyword or a numeric vector c(x, y).")
      }
      
      do.call(legend, legend_args)
    }
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_misspecification_cost_vs_epsilon <- function(
    grid_obj,                                                   # output of misspecification_cost_grid()
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"), # which matrix in grid_obj to plot on the y-axis
    delta_values = NULL,                                        # optional subset of delta values to display
    out_dir = file.path("figures", "misspecification"),         # directory where the figure is saved
    out_name = NULL,                                            # base filename for saved outputs
    show = interactive(),                                       # if TRUE, draw on the current device
    save = FALSE,                                               # if TRUE, save PNG and PDF copies
    width_in = 6.67,                                            # saved figure width in inches
    height_in = 4.67,                                           # saved figure height in inches
    dpi = 300,                                                  # PNG resolution used on save
    match_current = FALSE,                                      # if TRUE, save using the current device size instead
    title = TRUE,                                               # if TRUE, use the default plot title
    main = NULL,                                                # custom main title; overrides the default when supplied
    show_x_axis_title = TRUE,                                   # if TRUE, show the x-axis title
    show_y_axis_title = TRUE,                                   # if TRUE, show the y-axis title
    x_axis_title = expression(epsilon),                         # x-axis title
    y_axis_title = NULL,                                        # y-axis title; if NULL, use the default for the chosen surface
    axis_title_cex = 1,                                         # size of the axis titles
    tick_cex = 1.2,                                             # size of the axis tick labels
    title_cex = 1,                                              # size of the main title
    base_cex = 1,                                               # global cex scaling inside the panel
    mex = 1,                                                    # margin expansion factor
    margins = c(3.2, 3.6, 1.8, 0.8),                           # par("mar")
    outer_margins = c(0, 0, 0, 0),                              # par("oma")
    axis_mgp = c(2.2, 0.7, 0),                                 # par("mgp")
    tcl = -0.25,                                                # tick length
    axes = TRUE,                                                # if TRUE, draw axes
    box = TRUE,                                                 # if TRUE, draw a box around the panel
    x_axis_at = NULL,                                           # custom x-axis tick positions
    y_axis_at = NULL,                                           # custom y-axis tick positions
    x_axis_labels = TRUE,                                       # x-axis labels (logical or custom vector)
    y_axis_labels = TRUE,                                       # y-axis labels (logical or custom vector)
    axis_las = 1,                                               # orientation of axis tick labels
    xlim = NULL,                                                # custom x-range
    ylim = NULL,                                                # custom y-range
    xaxs = "i",                                                 # x-axis style used by plot.window()
    yaxs = "i",                                                 # y-axis style used by plot.window()
    add_grid = FALSE,                                           # if TRUE, add a 2D reference grid
    grid_col = "#DDDDDD",                                       # grid color
    grid_lty = 3,                                               # grid line type
    grid_lwd = 0.8,                                             # grid line width
    type = c("l", "p", "b", "o"),                              # line/point style used for the delta curves
    cols = NULL,                                                # line colors; generated automatically when NULL
    ltys = 1,                                                   # line types for the delta curves
    lwds = 2,                                                   # line widths for the delta curves
    add_points = FALSE,                                         # if TRUE, add points even when type = "l"
    point_pch = 16,                                             # point character for the delta curves
    point_cex = 0.7,                                            # point size for the delta curves
    point_cols = NULL,                                          # point colors; defaults to cols
    label_delta_on_curve = TRUE,                                # if TRUE, annotate each curve with delta directly on the line
    label_x = 0.9,                                              # target epsilon value where delta labels are placed
    label_gap = NULL,                                           # horizontal gap left in the line around the label; NULL -> automatic
    label_cex = 0.85,                                           # size of the inline delta labels
    label_col = NULL,                                           # label colors; defaults to the curve colors
    add_parametrizer_legend = TRUE,                             # if TRUE, draw a small broken-line key for the parametrizing variable
    add_legend = FALSE,                                         # if TRUE, show the delta legend
    legend_title = expression(delta),                           # legend title
    legend_labels = NULL,                                       # custom legend labels; if NULL, use the delta values
    legend_position = "topleft",                                # legend position keyword or c(x, y)
    legend_inset = 0.01,                                        # legend inset when using a keyword position
    legend_ncol = 1,                                            # number of legend columns
    legend_cex = 0.9,                                           # legend text size
    legend_bty = "n",                                           # legend box type
    legend_seg_len = 2.5,                                       # legend line-segment length
    restore_par = TRUE                                          # if FALSE, do not restore par() on exit; useful inside a shared layout
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  type <- match.arg(type)
  surface_meta <- .misspecification_surface_plot_meta(surface)
  
  margins <- .validate_margin_spec(margins, "margins")
  outer_margins <- .validate_margin_spec(outer_margins, "outer_margins")
  
  validate_range <- function(x, name) {
    x <- as.numeric(x)
    if (length(x) != 2L || any(!is.finite(x))) {
      stop(sprintf("'%s' must be a numeric vector of length 2.", name))
    }
    x <- sort(x)
    if (x[1] == x[2]) {
      stop(sprintf("'%s' must have distinct endpoints.", name))
    }
    x
  }
  
  if (is.null(grid_obj$delta_values) || is.null(grid_obj$eps_values)) {
    stop("grid_obj must be the output of misspecification_cost_grid().")
  }
  
  delta_grid <- as.numeric(grid_obj$delta_values)
  eps_values <- as.numeric(grid_obj$eps_values)
  zmat <- grid_obj[[surface]]
  
  zmat <- .ensure_misspecification_surface_has_finite_values(grid_obj, surface)
  if (!length(delta_grid) || !length(eps_values)) {
    stop("The selected grid must contain at least one delta value and one eps value.")
  }
  
  delta_ord <- order(delta_grid)
  eps_ord <- order(eps_values)
  delta_grid <- delta_grid[delta_ord]
  eps_values <- eps_values[eps_ord]
  zmat <- zmat[delta_ord, eps_ord, drop = FALSE]
  zmat[!is.finite(zmat)] <- NA_real_
  
  if (is.null(delta_values)) {
    delta_plot <- delta_grid
  } else {
    delta_values <- as.numeric(delta_values)
    if (!length(delta_values) || any(!is.finite(delta_values)) || any(delta_values < 0)) {
      stop("'delta_values' must be a non-empty numeric vector with non-negative entries.")
    }
    
    tol <- sqrt(.Machine$double.eps)
    idx <- vapply(
      delta_values,
      function(v) {
        hits <- which(abs(delta_grid - v) <= tol * max(1, abs(v)))
        if (length(hits)) hits[1] else NA_integer_
      },
      integer(1)
    )
    if (anyNA(idx)) {
      missing_delta <- format(signif(delta_values[is.na(idx)], 8), trim = TRUE, scientific = FALSE)
      stop(sprintf("The following delta_values were not found in grid_obj: %s",
                   paste(missing_delta, collapse = ", ")))
    }
    
    delta_plot <- delta_grid[idx]
    zmat <- zmat[idx, , drop = FALSE]
  }
  
  keep_curve <- apply(zmat, 1, function(y) any(is.finite(y)))
  if (!all(keep_curve)) {
    warning(sprintf("Dropping %d delta curve(s) with no finite values.", sum(!keep_curve)))
    delta_plot <- delta_plot[keep_curve]
    zmat <- zmat[keep_curve, , drop = FALSE]
  }
  if (!length(delta_plot) || !any(is.finite(zmat))) {
    stop("The selected surface does not contain any finite values to plot.")
  }
  
  if (is.null(xlim)) {
    xlim <- range(eps_values, finite = TRUE)
  } else {
    xlim <- validate_range(xlim, "xlim")
  }
  if (is.null(ylim)) {
    ylim <- range(zmat, finite = TRUE)
  } else {
    ylim <- validate_range(ylim, "ylim")
  }
  if (diff(ylim) == 0) {
    bump <- max(1, abs(ylim[1]))
    ylim <- ylim + c(-1, 1) * 0.05 * bump
  }
  
  label_x <- as.numeric(label_x)[1]
  if (!is.finite(label_x)) {
    stop("'label_x' must be a single finite number.")
  }
  if (!is.null(label_gap)) {
    label_gap <- as.numeric(label_gap)[1]
    if (!is.finite(label_gap) || label_gap < 0) {
      stop("'label_gap' must be NULL or a single non-negative number.")
    }
  }
  
  n_curves <- nrow(zmat)
  if (is.null(cols)) {
    cols <- tryCatch(
      grDevices::hcl.colors(n_curves, palette = "Dark 3"),
      error = function(e) grDevices::rainbow(n_curves, end = 0.85)
    )
  }
  cols <- rep_len(cols, n_curves)
  ltys <- rep_len(ltys, n_curves)
  lwds <- rep_len(lwds, n_curves)
  point_pch <- rep_len(point_pch, n_curves)
  point_cex <- rep_len(point_cex, n_curves)
  point_cols <- rep_len(point_cols %||% cols, n_curves)
  label_cols <- rep_len(label_col %||% cols, n_curves)
  
  if (is.null(legend_labels)) {
    legend_labels <- format(signif(delta_plot, 8), trim = TRUE, scientific = FALSE)
  } else if (length(legend_labels) != n_curves) {
    stop("'legend_labels' must have the same length as the number of plotted delta curves.")
  }
  
  if (is.null(y_axis_title)) {
    y_axis_title <- surface_meta$color_scale_title
  }
  if (is.null(main)) {
    main <- if (isTRUE(title)) paste(surface_meta$main, "vs epsilon") else ""
  }
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_vs_epsilon")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  show_lines <- type %in% c("l", "b", "o")
  show_points <- isTRUE(add_points) || type %in% c("p", "b", "o")
  
  draw_curve_with_gap <- function(x, y, col, lty, lwd, x_gap_left, x_gap_right) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    
    if (!length(x)) return(invisible(NULL))
    if (length(x) == 1L) {
      return(invisible(NULL))
    }
    
    if (!is.finite(x_gap_left) || !is.finite(x_gap_right) || x_gap_left >= x_gap_right ||
        x_gap_right <= min(x) || x_gap_left >= max(x)) {
      lines(x, y, col = col, lty = lty, lwd = lwd)
      return(invisible(NULL))
    }
    
    x_gap_left <- max(x_gap_left, min(x))
    x_gap_right <- min(x_gap_right, max(x))
    if (x_gap_left >= x_gap_right) {
      lines(x, y, col = col, lty = lty, lwd = lwd)
      return(invisible(NULL))
    }
    
    y_gap_left <- approx(x, y, xout = x_gap_left, ties = mean, rule = 1)$y
    y_gap_right <- approx(x, y, xout = x_gap_right, ties = mean, rule = 1)$y
    
    x_left <- x[x < x_gap_left]
    y_left <- y[x < x_gap_left]
    if (is.finite(y_gap_left) && x_gap_left > min(x)) {
      x_left <- c(x_left, x_gap_left)
      y_left <- c(y_left, y_gap_left)
    }
    
    x_right <- x[x > x_gap_right]
    y_right <- y[x > x_gap_right]
    if (is.finite(y_gap_right) && x_gap_right < max(x)) {
      x_right <- c(x_gap_right, x_right)
      y_right <- c(y_gap_right, y_right)
    }
    
    if (length(x_left) >= 2L) {
      lines(x_left, y_left, col = col, lty = lty, lwd = lwd)
    }
    if (length(x_right) >= 2L) {
      lines(x_right, y_right, col = col, lty = lty, lwd = lwd)
    }
    
    invisible(NULL)
  }
  
  curve_label_info <- function(x, y, target_x, gap_width) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    
    if (length(x) < 2L) return(NULL)
    if (target_x < min(x) || target_x > max(x)) return(NULL)
    
    y_target <- approx(x, y, xout = target_x, ties = mean, rule = 1)$y
    if (!is.finite(y_target)) return(NULL)
    
    left_x <- max(min(x), target_x - gap_width/2)
    right_x <- min(max(x), target_x + gap_width/2)
    if (left_x >= right_x) return(NULL)
    
    list(x = target_x, y = y_target, left_x = left_x, right_x = right_x)
  }
  
  plotfun <- function() {
    if (isTRUE(restore_par)) {
      op <- par(c("oma", "mar", "mgp", "tcl", "cex", "mex",
                  "cex.axis", "cex.lab", "cex.main"))
      on.exit(par(op), add = TRUE)
    }
    
    par_args <- list(
      mar = margins,
      mgp = axis_mgp,
      tcl = tcl,
      cex = base_cex,
      mex = mex,
      cex.axis = tick_cex,
      cex.lab = axis_title_cex,
      cex.main = title_cex
    )
    if (isTRUE(restore_par)) {
      par_args$oma <- outer_margins
    }
    do.call(par, par_args)
    
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
    
    if (isTRUE(add_grid)) {
      grid(col = grid_col, lty = grid_lty, lwd = grid_lwd)
    }
    
    for (j in seq_len(n_curves)) {
      yj <- zmat[j, ]
      okj <- is.finite(yj)
      if (!any(okj)) next
      
      xj <- eps_values[okj]
      yj <- yj[okj]
      label_info <- NULL
      if (isTRUE(label_delta_on_curve)) {
        gap_width_j <- label_gap
        if (is.null(gap_width_j)) {
          gap_width_j <- max(
            1.15 * strwidth(legend_labels[j], cex = label_cex, units = "user"),
            0.015 * diff(xlim)
          )
        }
        label_info <- curve_label_info(
          x = xj,
          y = yj,
          target_x = label_x,
          gap_width = gap_width_j
        )
      }
      
      if (show_lines) {
        if (is.null(label_info)) {
          lines(xj, yj, col = cols[j], lty = ltys[j], lwd = lwds[j])
        } else {
          draw_curve_with_gap(
            x = xj,
            y = yj,
            col = cols[j],
            lty = ltys[j],
            lwd = lwds[j],
            x_gap_left = label_info$left_x,
            x_gap_right = label_info$right_x
          )
        }
      }
      if (show_points) {
        if (!is.null(label_info)) {
          keep_pts <- (xj < label_info$left_x) | (xj > label_info$right_x)
          x_pts <- xj[keep_pts]
          y_pts <- yj[keep_pts]
        } else {
          x_pts <- xj
          y_pts <- yj
        }
        if (length(x_pts)) {
          points(x_pts, y_pts, pch = point_pch[j], cex = point_cex[j], col = point_cols[j])
        }
      }
      if (!is.null(label_info)) {
        text(
          x = label_info$x,
          y = label_info$y,
          labels = legend_labels[j],
          cex = label_cex,
          col = label_cols[j],
          xpd = NA
        )
      }
    }
    
    if (isTRUE(axes)) {
      axis(1, at = x_axis_at, labels = x_axis_labels, las = axis_las)
      axis(2, at = y_axis_at, labels = y_axis_labels, las = axis_las)
    }
    if (isTRUE(box)) {
      box()
    }
    
    title(
      main = main,
      xlab = if (show_x_axis_title) x_axis_title else "",
      ylab = if (show_y_axis_title) y_axis_title else ""
    )
    
    if (isTRUE(add_parametrizer_legend)) {
      .draw_parametrizer_legend_key(
        label = expression(delta),
        cex = tick_cex
      )
    }
    
    if (isTRUE(add_legend)) {
      legend_args <- list(
        legend = legend_labels,
        title = legend_title,
        col = cols,
        ncol = legend_ncol,
        cex = legend_cex,
        bty = legend_bty,
        seg.len = legend_seg_len,
        merge = TRUE
      )
      
      if (show_lines) {
        legend_args$lty <- ltys
        legend_args$lwd <- lwds
      }
      if (show_points) {
        legend_args$pch <- point_pch
        legend_args$pt.cex <- point_cex
      }
      
      if (is.character(legend_position) && length(legend_position) == 1L) {
        legend_args$x <- legend_position
        legend_args$inset <- legend_inset
      } else if (is.numeric(legend_position) && length(legend_position) == 2L &&
                 all(is.finite(legend_position))) {
        legend_args$x <- legend_position[1]
        legend_args$y <- legend_position[2]
      } else {
        stop("'legend_position' must be a single legend keyword or a numeric vector c(x, y).")
      }
      
      do.call(legend, legend_args)
    }
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

plot_misspecification_cost_canvas <- function(
    grid_obj,
    grid_obj_vs_delta = grid_obj,
    grid_obj_vs_epsilon = grid_obj,
    surface = c("relative_cost_pct", "relative_cost", "gamma_difference",
                "gamma_with_ambiguity", "gamma_zero_ambiguity_policy"),
    out_dir = file.path("figures", "misspecification"),
    out_name = NULL,
    show = interactive(),
    save = FALSE,
    width_in = 12,
    height_in = 8.5,
    dpi = 300,
    match_current = FALSE,
    panel_layout = c("1x3", "3x1"),
    column_widths = NULL,
    row_heights = NULL,
    contour_args = list(),
    vs_delta_args = list(),
    vs_epsilon_args = list()
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  panel_layout <- match.arg(panel_layout)
  
  validate_positive_vector <- function(x, name, n) {
    x <- as.numeric(x)
    if (length(x) != n || any(!is.finite(x)) || any(x <= 0)) {
      stop(sprintf("'%s' must be a numeric vector of length %d with positive entries.", name, n))
    }
    x
  }
  validate_list <- function(x, name) {
    if (!is.list(x)) {
      stop(sprintf("'%s' must be a list.", name))
    }
    x
  }
  
  if (identical(panel_layout, "1x3")) {
    if (is.null(column_widths)) column_widths <- c(1, 1, 1)
    if (is.null(row_heights)) row_heights <- 1
    column_widths <- validate_positive_vector(column_widths, "column_widths", 3L)
    row_heights <- validate_positive_vector(row_heights, "row_heights", 1L)
  } else {
    if (is.null(column_widths)) column_widths <- 1
    if (is.null(row_heights)) row_heights <- c(1, 1, 1)
    column_widths <- validate_positive_vector(column_widths, "column_widths", 1L)
    row_heights <- validate_positive_vector(row_heights, "row_heights", 3L)
  }
  contour_args <- validate_list(contour_args, "contour_args")
  vs_delta_args <- validate_list(vs_delta_args, "vs_delta_args")
  vs_epsilon_args <- validate_list(vs_epsilon_args, "vs_epsilon_args")
  
  if (is.null(out_name)) {
    out_name <- paste0((grid_obj$out_name %||% "misspecification"), "_", surface, "_canvas")
  }
  if (missing(out_dir) && !is.null(grid_obj$out_dir)) {
    out_dir <- grid_obj$out_dir
  }
  
  build_panel_args <- function(defaults, extras, grid) {
    args <- utils::modifyList(defaults, extras)
    args$grid_obj <- grid
    args$surface <- surface
    args$out_dir <- out_dir
    args$show <- TRUE
    args$save <- FALSE
    args$title <- FALSE
    args$restore_par <- FALSE
    args
  }
  
  contour_call <- build_panel_args(
    defaults = list(add_color_scale = FALSE),
    extras = contour_args,
    grid = grid_obj
  )
  contour_call$add_color_scale <- FALSE
  
  vs_delta_call <- build_panel_args(
    defaults = list(),
    extras = vs_delta_args,
    grid = grid_obj_vs_delta
  )
  vs_epsilon_call <- build_panel_args(
    defaults = list(),
    extras = vs_epsilon_args,
    grid = grid_obj_vs_epsilon
  )
  
  plotfun <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    layout_mat <- if (identical(panel_layout, "1x3")) {
      matrix(1:3, nrow = 1, byrow = TRUE)
    } else {
      matrix(1:3, ncol = 1, byrow = TRUE)
    }
    layout(layout_mat, widths = column_widths, heights = row_heights)
    on.exit(layout(1), add = TRUE)
    
    do.call(plot_misspecification_cost_vs_delta, vs_delta_call)
    do.call(plot_misspecification_cost_contour, contour_call)
    do.call(plot_misspecification_cost_vs_epsilon, vs_epsilon_call)
  }
  
  render_and_save(out_name, plotfun, show = show, save = save, dir = out_dir,
                  width_in = width_in, height_in = height_in, dpi = dpi,
                  match_current = match_current)
}

fit_best_family <- function(
    x, y,
    families  = c("monomial", "exponential"),
    criterion = c("AIC", "RSS"),
    verbose   = TRUE,
    max_tries = 6              # how many jittered tries for each nonlinear family
) {
  criterion <- match.arg(criterion)
  families  <- intersect(families, c("monomial", "exponential"))
  if (length(families) == 0L) stop("No valid families given.")
  
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (!all(ok)) {
    if (verbose) message("Dropping ", sum(!ok), " non-finite (x,y) pairs.")
    x <- x[ok]; y <- y[ok]
  }
  if (length(x) < 3L) stop("Need at least 3 data points.")
  
  df <- data.frame(x = x, y = y)
  xrange <- range(df$x)
  xspan  <- diff(xrange)
  if (xspan <= 0) stop("x must not be constant.")
  
  # -------- safe building blocks to avoid NaN/Inf -----------------
  safe_pow <- function(x, b, c, eps = 1e-8, max_power = 50) {
    # base >= eps so no negative or zero to weird powers
    base <- pmax(x - b, eps)
    val  <- base^c
    val[!is.finite(val)] <- max_power
    pmin(val, max_power)
  }
  
  safe_exp_pow <- function(x, b, c, eps = 1e-8, max_exponent = 700) {
    base <- pmax(x - b, eps)
    power <- base^c
    power[!is.finite(power)] <- max_exponent
    exp(pmin(power, max_exponent))
  }
  
  fits <- list()
  crit <- data.frame(
    family = character(),
    RSS    = numeric(),
    AIC    = numeric(),
    stringsAsFactors = FALSE
  )
  
  use_nlsLM <- requireNamespace("minpack.lm", quietly = TRUE)
  if (verbose && use_nlsLM) {
    message("Using minpack.lm::nlsLM for nonlinear fits.")
  }
  
  # ---------- helper: robust nonlinear fitter with jittered starts ----------
  fit_nonlinear_family <- function(formula, start, lower, upper, family_label) {
    best_fit <- NULL
    best_RSS <- Inf
    last_err <- NULL
    
    for (k in seq_len(max_tries)) {
      st <- start
      
      # Jitter after the first attempt
      if (k > 1) {
        for (nm in names(st)) {
          v <- st[[nm]]
          if (!is.finite(v)) next
          j  <- rnorm(1, mean = 0, sd = 0.3)
          st[[nm]] <- v * (1 + j)
        }
        # enforce bounds
        for (nm in names(st)) {
          if (!is.null(lower) && nm %in% names(lower)) {
            st[[nm]] <- max(st[[nm]], lower[[nm]])
          }
          if (!is.null(upper) && nm %in% names(upper)) {
            st[[nm]] <- min(st[[nm]], upper[[nm]])
          }
        }
      }
      
      fit_try <- try(
        if (use_nlsLM) {
          minpack.lm::nlsLM(
            formula, data = df, start = st,
            lower = lower, upper = upper,
            control = minpack.lm::nls.lm.control(maxiter = 1000)
          )
        } else {
          nls(
            formula, data = df, start = st,
            algorithm = "port",
            lower = lower, upper = upper,
            control = list(maxiter = 500)
          )
        },
        silent = TRUE
      )
      
      if (!inherits(fit_try, "try-error")) {
        RSS <- sum(resid(fit_try)^2)
        if (RSS < best_RSS) {
          best_RSS <- RSS
          best_fit <- fit_try
        }
      } else {
        last_err <- fit_try
      }
    }
    
    if (is.null(best_fit) && verbose && !is.null(last_err)) {
      cond <- attr(last_err, "condition")
      msg  <- if (inherits(cond, "condition")) conditionMessage(cond) else as.character(last_err)
      message("  ", family_label, " fit failed after ", max_tries, " tries: ", msg)
    }
    best_fit
  }
  
  ## ------------- 2) Monomial: y = d + a * safe_pow(x, b, c) -------------
  if ("monomial" %in% families) {
    if (verbose) message("Fitting monomial model y = d + a*safe_pow(x, b, c) ...")
    
    # base starting values
    b0 <- xrange[1] - 0.1 * xspan      # somewhere to the left
    t0 <- safe_pow(df$x, b0, 1)        # essentially x-b0
    d0 <- min(df$y)
    
    # crude a0, c0 via log–log if y-d0>0 and t0>0
    if (all(t0 > 0)) {
      y_shift <- df$y - d0
      if (all(y_shift > 0)) {
        llm <- try(lm(log(y_shift) ~ log(t0), data = df), silent = TRUE)
        if (!inherits(llm, "try-error")) {
          c0  <- coef(llm)[2]
          a0  <- exp(coef(llm)[1])
        } else {
          a0 <- mean(y_shift)
          c0 <- 1
        }
      } else {
        a0 <- mean(df$y)
        c0 <- 1
      }
    } else {
      a0 <- mean(df$y)
      c0 <- 1
    }
    
    # bounds: c in [0, 5] to avoid wild curvature; b within a generous x-range; d free
    lower_m <- c(a = -Inf, b = xrange[1] - 5*xspan, c = -5,   d = -Inf)
    upper_m <- c(a =  Inf, b = xrange[2] + 5*xspan, c = 5, d =  Inf)
    start_m <- list(a = a0, b = b0, c = c0, d = d0)
    
    m_mono <- fit_nonlinear_family(
      y ~ d + a * safe_pow(x, b, c),
      start = start_m,
      lower = lower_m,
      upper = upper_m,
      family_label = "Monomial"
    )
    
    if (!is.null(m_mono)) {
      RSS_m <- sum(resid(m_mono)^2)
      AIC_m <- AIC(m_mono)
      fits$monomial <- m_mono
      crit <- rbind(crit,
                    data.frame(family = "monomial", RSS = RSS_m, AIC = AIC_m))
    }
  }
  
  ## --- 3) Exponential: y = d + a * safe_exp_pow(x, b, c) (or decay version) ---
  if ("exponential" %in% families) {
    if (verbose) message("Fitting exponential model y = d + a*safe_exp_pow(x, b, c) ...")
    
    b0 <- xrange[1] - 0.1 * xspan
    d0 <- min(df$y)
    a0 <- max(df$y) - d0
    if (!is.finite(a0) || a0 == 0) a0 <- mean(df$y) - d0
    c0 <- 1
    
    lower_e <- c(a = -Inf, b = xrange[1] - 5*xspan, c = -Inf,   d = -Inf)
    upper_e <- c(a =  Inf, b = xrange[2] + 5*xspan, c = 5.0, d =  Inf)
    start_e <- list(a = a0, b = b0, c = c0, d = d0)
    
    m_exp <- fit_nonlinear_family(
      # If you want explicit decay, switch the model to:
      # y ~ d + a * (1 / safe_exp_pow(x, b, c))
      y ~ d + a * safe_exp_pow(x, b, c),
      start = start_e,
      lower = lower_e,
      upper = upper_e,
      family_label = "Exponential"
    )
    
    if (!is.null(m_exp)) {
      RSS_e <- sum(resid(m_exp)^2)
      AIC_e <- AIC(m_exp)
      fits$exponential <- m_exp
      crit <- rbind(crit,
                    data.frame(family = "exponential", RSS = RSS_e, AIC = AIC_e))
    }
  }
  
  if (nrow(crit) == 0L) stop("No model could be fitted.")
  
  ## ------------------------- Choose best model -------------------------
  best_row <- if (criterion == "AIC") which.min(crit$AIC) else which.min(crit$RSS)
  best_family <- crit$family[best_row]
  best_fit    <- fits[[best_family]]
  
  if (verbose) {
    message(sprintf("Best family: %s (criterion = %s)", best_family, criterion))
    print(crit)
  }
  
  list(
    best_family = best_family,
    best_fit    = best_fit,
    all_fits    = fits,
    criteria    = crit
  )
}











# # ---------------------- EXAMPLE -----------------------------------------------
# Parameters
p <- make_params(b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1, l = 1)
xL0 <- -0.5; xk0 <- -0.1; xl0 <- 0.5; xU0 <- 1

# Suboptimal H (inner solver)
subopt_H <- build_suboptimal_H(p, xL0, xk0, xl0, xU0)
diagnose(subopt_H)
# Plot H
plot_H(subopt_H, p, show = TRUE, save = FALSE)
plot_H_prime(subopt_H, p, show = TRUE, save = FALSE)

# Optimal H (outer solver)
sol_opt <- solve_optimal_barriers(p, xL0, xk0, xl0, xU0, verbose=TRUE, 
                                  switch_to_newton = FALSE, tol_regime_switch = 1e-3)
diagnose(sol_opt)
# Plot H and simulate the reflected process
plot_H(sol_opt, p, show = TRUE, save = TRUE, 
       top_blank = 0.075, bottom_blank = 0.05, 
       margins = c(3, 2, 1, 1),
       show_x_axis_title = TRUE, show_y_axis_title = FALSE, 
       axis_mgp = c(2.2, 0.7, 0), axis_title_cex = 1.4,
       tick_cex = 1)
plot_H_prime(sol_opt, p, show = TRUE, save = FALSE)
sim <- simulate_reflected_jd(params=p, thresholds=sol_opt$x, seed = 123)
plot_reflected_jd(sol_opt, p, sim=sim, show=TRUE, save=TRUE)
plot_controls(sim, show=TRUE, save=TRUE)
plot_reflected_with_controls(sol_opt, p, sim, top_blank = 0, bottom_blank = 0, 
                             heights = c(0.7, 1.7, 0.85), draw_legend = FALSE,
                             margins = c(3, 2, 1, 1), axis_mgp = c(2.2, 0.7, 0),
                             show_x_axis_title = TRUE, x_axis_title = "t",
                             show_y_axis_title = FALSE, axis_title_cex = 1.4,
                             tick_cex = 1,
                             show=TRUE, save=TRUE)

# # ---------------------- Cost of misspecification on a (delta, eps) grid ---------------------
# Cost of misspecification on a (delta, eps) grid
.format_misspecification_param_value <- function(x) {
  gsub("\\s+", "", format(signif(as.numeric(x), 8), trim = TRUE, scientific = FALSE))
}

.format_misspecification_tag_value <- function(x) {
  x <- .format_misspecification_param_value(x)
  x <- gsub("-", "minus_", x, fixed = TRUE)
  x <- gsub("\\.", "p", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

.misspecification_example_tag <- function(b_value, r) {
  paste0(
    "b_", .format_misspecification_tag_value(b_value),
    "__r_", .format_misspecification_tag_value(r)
  )
}

.resolve_misspecification_example_tag <- function(tag = NULL, b_value = NULL, r = NULL) {
  if (!is.null(tag)) {
    tag <- as.character(tag)[1]
    if (nzchar(tag)) {
      return(tag)
    }
  }
  if (is.null(b_value) || is.null(r)) {
    stop("Provide either 'tag' or both 'b_value' and 'r' to identify the misspecification example.")
  }
  .misspecification_example_tag(b_value = b_value, r = r)
}

.misspecification_param_suffix <- function(sigma, mu, u, l) {
  paste0(
    "__sigma_", .format_misspecification_param_value(sigma),
    "__mu_", .format_misspecification_param_value(mu),
    "__u_", .format_misspecification_param_value(u),
    "__l_", .format_misspecification_param_value(l)
  )
}

.misspecification_cost_example_spec <- function(
    tag = NULL, b_value = NULL,
    r = 1, sigma = 1, mu = 1, u = 1, l = 1
) {
  tag <- .resolve_misspecification_example_tag(tag = tag, b_value = b_value, r = r)
  base_name <- paste0("misspec_", tag, .misspecification_param_suffix(sigma, mu, u, l))
  base_out_dir <- file.path("figures", "misspecification", tag)
  
  list(
    out_dir = base_out_dir,
    base_name = base_name,
    main_grid = list(
      delta_values = seq(0, 1, by = 0.01),
      eps_values = seq(0, 1, by = 0.01),
      out_name = paste0(base_name, "_surface")
    ),
    grid_vs_delta = list(
      delta_values = seq(0, 1, by = 0.025),
      eps_values = seq(0, 1, by = 0.2),
      out_name = paste0(base_name, "_vs_delta_data")
    ),
    grid_vs_epsilon = list(
      delta_values = seq(0, 1, by = 0.2),
      eps_values = seq(0, 1, by = 0.025),
      out_name = paste0(base_name, "_vs_epsilon_data")
    ),
    text_grid = list(
      delta_values = seq(0, 1, by = 0.2),
      eps_values = seq(0, 1, by = 0.2),
      out_name = paste0(base_name, "_table_data")
    )
  )
}

generate_misspecification_cost_example_data <- function(
    b_value, tag = NULL,
    r = 1, sigma = 1, mu = 1, u = 1, l = 1,
    force = FALSE
) {
  tag <- .resolve_misspecification_example_tag(tag = tag, b_value = b_value, r = r)
  spec <- .misspecification_cost_example_spec(
    tag = tag,
    b_value = b_value,
    r = r, sigma = sigma, mu = mu, u = u, l = l
  )
  grid_names <- c("main_grid", "grid_vs_delta", "grid_vs_epsilon", "text_grid")
  
  grid_common_args <- list(
    b = b_value, r = r, sigma = sigma, mu = mu, u = u, l = l,
    xL0 = -0.5, xk0 = -0.1, xl0 = 0.5, xU0 = 1,
    verbose = FALSE,
    out_dir = spec$out_dir,
    save = TRUE
  )
  
  out_paths <- setNames(vector("list", length(grid_names)), grid_names)
  
  for (grid_name in grid_names) {
    grid_args <- spec[[grid_name]]
    out_path <- .misspecification_rdata_path(spec$out_dir, grid_args$out_name)
    
    if (file.exists(out_path) && !isTRUE(force)) {
      message(sprintf("Using cached misspecification grid: %s", normalizePath(out_path)))
      out_paths[[grid_name]] <- out_path
      next
    }
    
    grid_obj <- do.call(
      misspecification_cost_grid,
      c(grid_args, grid_common_args)
    )
    out_paths[[grid_name]] <- grid_obj$rdata_path
  }
  
  invisible(out_paths)
}

load_misspecification_cost_example_data <- function(
    tag = NULL,
    b_value = NULL,
    r = 1, sigma = 1, mu = 1, u = 1, l = 1
) {
  tag <- .resolve_misspecification_example_tag(tag = tag, b_value = b_value, r = r)
  spec <- .misspecification_cost_example_spec(
    tag = tag,
    b_value = b_value,
    r = r, sigma = sigma, mu = mu, u = u, l = l
  )
  load_one <- function(grid_args) {
    load_misspecification_cost_rdata(
      out_dir = spec$out_dir,
      out_name = grid_args$out_name
    )
  }
  
  list(
    main_grid = load_one(spec$main_grid),
    grid_vs_delta = load_one(spec$grid_vs_delta),
    grid_vs_epsilon = load_one(spec$grid_vs_epsilon),
    text_grid = load_one(spec$text_grid)
  )
}

.misspecification_grid_model_params <- function(grid_obj) {
  cp <- grid_obj$certainty_params
  needed <- c("b", "r", "sigma", "mu", "u", "l")
  if (is.null(cp) || !all(needed %in% names(cp))) {
    return(NULL)
  }
  as.list(cp[needed])
}

plot_misspecification_cost_example_set <- function(
    tag = NULL,
    b_value = NULL,
    grids = NULL,
    r = 1,
    sigma = 1,
    mu = 1,
    u = 1,
    l = 1,
    surface = c("relative_cost_pct", "relative_cost"),
    table_surface = NULL,
    vs_delta_label_x = 0.9,
    canvas_panel_layout = "1x3",
    canvas_column_widths = NULL,
    canvas_row_heights = NULL,
    canvas_width_in = NULL,
    canvas_height_in = NULL
) {
  surface <- .misspecification_reporting_surface(match.arg(surface))
  if (is.null(table_surface)) {
    table_surface <- surface
  }
  table_surface <- .misspecification_reporting_surface(
    match.arg(table_surface, c("relative_cost_pct", "relative_cost"))
  )
  
  if (!is.null(grids)) {
    grid_params <- .misspecification_grid_model_params(grids$main_grid)
    if (!is.null(grid_params)) {
      b_value <- grid_params$b
      r <- grid_params$r
      sigma <- grid_params$sigma
      mu <- grid_params$mu
      u <- grid_params$u
      l <- grid_params$l
    }
  }
  
  tag <- .resolve_misspecification_example_tag(tag = tag, b_value = b_value, r = r)
  
  if (is.null(grids)) {
    grids <- load_misspecification_cost_example_data(
      tag = tag,
      b_value = b_value,
      r = r, sigma = sigma, mu = mu, u = u, l = l
    )
  }
  
  spec <- .misspecification_cost_example_spec(
    tag = tag,
    b_value = b_value,
    r = r, sigma = sigma, mu = mu, u = u, l = l
  )
  base_out_dir <- spec$out_dir
  base_name <- spec$base_name
  persp_r <- 0.5
  
  misspec_grid <- grids$main_grid
  misspec_grid_cost_vs_delta <- grids$grid_vs_delta
  misspec_grid_cost_vs_epsilon <- grids$grid_vs_epsilon
  misspec_grid_text <- grids$text_grid
  
  if (is.null(canvas_width_in)) {
    canvas_width_in <- if (identical(canvas_panel_layout, "3x1")) 4.67 else 12
  }
  if (is.null(canvas_height_in)) {
    canvas_height_in <- if (identical(canvas_panel_layout, "3x1")) 12 else 4.67
  }
  if (is.null(canvas_column_widths)) {
    canvas_column_widths <- if (identical(canvas_panel_layout, "3x1")) 1 else c(1, 1, 1)
  }
  if (is.null(canvas_row_heights)) {
    canvas_row_heights <- if (identical(canvas_panel_layout, "3x1")) c(1, 1, 1) else 1
  }
  
  plot_misspecification_cost_3d(
    misspec_grid,
    surface = surface,
    show = TRUE,
    save = TRUE,
    add_color_scale = FALSE,
    height_in = 4.67, width_in = 6.67,
    margins = c(0.9, 2.1, 0.0, 0.0),
    x_axis_title_pad = 0.12,
    y_axis_title_pad = 0.12,
    z_axis_title_pad = 0.1,
    x_axis_title_at = 0.5,
    y_axis_title_at = 0.5,
    z_axis_title_at = 0.5,
    r = persp_r,
    expand = 0.65,
    title = FALSE,
    tick_cex = 1,
    axis_title_cex = 1.4
  )
  plot_misspecification_cost_contour(
    misspec_grid,
    surface = surface,
    show = TRUE,
    save = TRUE,
    title = FALSE,
    show_x_axis_title = TRUE,
    show_y_axis_title = TRUE,
    axis_title_cex = 1.4,
    tick_cex = 1,
    margins = c(3.5, 3.8, 1.5, 1),
    axis_mgp = c(2.6, 0.8, 0),
    outer_margins = c(0, 0, 0, 0),
    add_grid = TRUE,
    contour_nlevels = 12,
    use_z_colormap = FALSE,
    add_color_scale = FALSE,
    contour_lwd = 1.5, contour_labcex = 1
  )
  plot_misspecification_cost_vs_delta(
    misspec_grid_cost_vs_delta,
    surface = surface,
    show = TRUE,
    save = TRUE,
    title = FALSE,
    show_x_axis_title = TRUE,
    show_y_axis_title = TRUE,
    axis_title_cex = 1.4,
    tick_cex = 1,
    margins = c(3.8, 3.8, 1.5, 1),
    axis_mgp = c(2.6, 0.8, 0),
    add_grid = TRUE,
    label_eps_on_curve = TRUE,
    label_x = vs_delta_label_x,
    add_legend = FALSE, cols = "black"
  )
  plot_misspecification_cost_vs_epsilon(
    misspec_grid_cost_vs_epsilon,
    surface = surface,
    show = TRUE,
    save = TRUE,
    title = FALSE,
    show_x_axis_title = TRUE,
    show_y_axis_title = TRUE,
    axis_title_cex = 1.4,
    tick_cex = 1,
    margins = c(3.8, 3.8, 1.5, 1),
    axis_mgp = c(2.6, 0.8, 0),
    add_grid = TRUE,
    label_delta_on_curve = TRUE,
    label_x = 0.9,
    add_legend = FALSE, cols = "black"
  )
  plot_misspecification_cost_canvas(
    misspec_grid,
    grid_obj_vs_delta = misspec_grid_cost_vs_delta,
    grid_obj_vs_epsilon = misspec_grid_cost_vs_epsilon,
    surface = surface,
    out_dir = base_out_dir,
    out_name = paste0(base_name, "_", surface, "_canvas"),
    show = TRUE,
    save = TRUE,
    width_in = canvas_width_in,
    height_in = canvas_height_in,
    panel_layout = canvas_panel_layout,
    column_widths = canvas_column_widths,
    row_heights = canvas_row_heights,
    contour_args = list(
      show_x_axis_title = TRUE,
      show_y_axis_title = TRUE,
      axis_title_cex = 1.1,
      tick_cex = 0.9,
      margins = c(3.2, 3.4, 0.8, 0.5),
      axis_mgp = c(2.2, 0.7, 0),
      outer_margins = c(0, 0, 0, 0),
      add_grid = TRUE,
      contour_nlevels = 12,
      contour_lwd = 1.3,
      contour_labcex = 0.9
    ),
    vs_delta_args = list(
      show_x_axis_title = TRUE,
      show_y_axis_title = TRUE,
      axis_title_cex = 1.1,
      tick_cex = 0.9,
      margins = c(3.2, 3.4, 0.8, 0.5),
      axis_mgp = c(2.2, 0.7, 0),
      outer_margins = c(0, 0, 0, 0),
      add_grid = TRUE,
      label_eps_on_curve = TRUE,
      label_x = vs_delta_label_x,
      add_legend = FALSE,
      cols = "black"
    ),
    vs_epsilon_args = list(
      show_x_axis_title = TRUE,
      show_y_axis_title = TRUE,
      axis_title_cex = 1.1,
      tick_cex = 0.9,
      margins = c(3.2, 3.4, 0.8, 0.5),
      axis_mgp = c(2.2, 0.7, 0),
      outer_margins = c(0, 0, 0, 0),
      add_grid = TRUE,
      label_delta_on_curve = TRUE,
      label_x = 0.9,
      add_legend = FALSE,
      cols = "black"
    )
  )
  
  tex_obj <- export_misspecification_cost_latex(
    misspec_grid_text,
    surface = table_surface,
    out_dir = base_out_dir,
    out_name = paste0("rmc_table_", tag),
    save = TRUE
  )
  
  invisible(list(
    main_grid = misspec_grid,
    grid_vs_delta = misspec_grid_cost_vs_delta,
    grid_vs_epsilon = misspec_grid_cost_vs_epsilon,
    text_grid = misspec_grid_text,
    tex = tex_obj
  ))
}

# Build the expensive grids once; rerun with force = TRUE when either the grid
# setup or the model parameters change.
misspec_example_common_params <- list(
  sigma = 1,
  mu = 1,
  u = 1,
  l = 1
)

misspec_example_specs <- list(
  list(
    b_value = 2,
    r = 1,
    vs_delta_label_x = 0.9,
    canvas_panel_layout = "3x1",
    canvas_row_heights = c(1, 1, 1)
  ),
  list(
    b_value = -2,
    r = 1,
    vs_delta_label_x = 0.2,
    canvas_panel_layout = "3x1",
    canvas_row_heights = c(1, 1, 1)
  )
)

run_misspecification_cost_example <- function(
    example_spec,
    common_params = misspec_example_common_params,
    force = FALSE
) {
  stopifnot(is.list(example_spec), !is.null(example_spec$b_value))
  
  model_param_names <- c("r", "sigma", "mu", "u", "l")
  plot_param_names <- c(
    "surface", "table_surface", "vs_delta_label_x",
    "canvas_panel_layout", "canvas_column_widths", "canvas_row_heights",
    "canvas_width_in", "canvas_height_in"
  )
  
  example_model_params <- utils::modifyList(
    common_params,
    example_spec[intersect(names(example_spec), model_param_names)]
  )
  example_model_params$r <- example_model_params$r %||% 1
  example_plot_params <- example_spec[intersect(names(example_spec), plot_param_names)]
  tag <- .resolve_misspecification_example_tag(
    tag = example_spec$tag %||% NULL,
    b_value = example_spec$b_value,
    r = example_model_params$r
  )
  
  paths <- do.call(
    generate_misspecification_cost_example_data,
    c(
      list(
        b_value = example_spec$b_value,
        tag = tag,
        force = force
      ),
      example_model_params
    )
  )
  
  # Reload cached .RData grids when only plot styling changes.
  result <- do.call(
    plot_misspecification_cost_example_set,
    c(
      list(tag = tag, b_value = example_spec$b_value),
      example_model_params,
      example_plot_params
    )
  )
  
  list(
    spec = utils::modifyList(
      example_model_params,
      utils::modifyList(example_spec, list(tag = tag))
    ),
    paths = paths,
    result = result
  )
}

misspec_example_runs <- lapply(
  misspec_example_specs,
  run_misspecification_cost_example
)
names(misspec_example_runs) <- vapply(
  misspec_example_runs,
  function(run) run$spec$tag,
  character(1)
)

cat(paste(
  vapply(misspec_example_runs, function(run) run$result$tex$latex, character(1)),
  collapse = "\n\n"
))



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
# b = 3
sweep_u <- comparative_sweeper(
  sweep_param  = "u",
  sweep_values = seq(0.01, 5, by = 0.01),
  b = 3, delta = 1, r = 1, eps = 0.5, sigma = 1, mu = 1, l = 1,
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
# b = -3
sweep_u <- comparative_sweeper(
  sweep_param  = "u",
  sweep_values = seq(0.01, 5, by = 0.01),
  b = -3, delta = 1, r = 1, eps = 0.5, sigma = 1, mu = 1, l = 1,
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
# b = 3
sweep_l <- comparative_sweeper(
  sweep_param  = "l",
  sweep_values = seq(1, 5, by = 0.01),
  b = 3, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1,
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
# b = -3
sweep_l <- comparative_sweeper(
  sweep_param  = "l",
  sweep_values = seq(1, 5, by = 0.01),
  b = -5, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1,
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

