# ============================================================
# STAT 438 — Best Linear Predictor + Partial Correlation (Given mu, S)
# Solves questions like:
# 1) β = Σ_YZ Σ_ZZ^{-1}
# 2) β0 = μ_Y − β μ_Z  and write y-hat equations
# 3) Σ_YY·Z = Σ_YY − Σ_YZ Σ_ZZ^{-1} Σ_ZY
# 4) ρ_{y1,y2·Z} = σ12·Z / sqrt(σ11·Z σ22·Z)
#
# Works for:
# - 2 responses (Y1,Y2) and 1+ predictors (Z1, Z2, ...)
# - or any number of responses (>=2) and predictors (>=1)
#
# You ONLY need: mean vector mu and covariance matrix S,
# plus which variables are Y and which are Z.
# ============================================================

stat438_muS_predictor <- function(mu, S,
                                  y_idx = 1:2,      # which positions in mu/S correspond to Y's
                                  z_idx = 3:4,      # which positions correspond to Z's
                                  digits = 6,
                                  show_equations = TRUE,
                                  names = NULL      # optional names for variables in order of mu/S
) {
  
  mu <- as.numeric(mu)
  S <- as.matrix(S)
  
  if (nrow(S) != ncol(S)) stop("S must be a square covariance matrix.")
  if (length(mu) != nrow(S)) stop("mu length must match dimensions of S.")
  if (length(intersect(y_idx, z_idx)) > 0) stop("y_idx and z_idx must not overlap.")
  if (length(y_idx) < 2) stop("Need at least 2 Y variables to compute partial correlation.")
  if (length(z_idx) < 1) stop("Need at least 1 Z variable.")
  
  p <- length(mu)
  
  if (is.null(names)) {
    var_names <- paste0("V", 1:p)
  } else {
    if (length(names) != p) stop("names must have the same length as mu.")
    var_names <- names
  }
  
  Y_names <- var_names[y_idx]
  Z_names <- var_names[z_idx]
  
  # ---- Partition mu ----
  mu_Y <- matrix(mu[y_idx], ncol = 1)
  mu_Z <- matrix(mu[z_idx], ncol = 1)
  
  # ---- Partition S into blocks ----
  S_YY <- S[y_idx, y_idx, drop = FALSE]
  S_YZ <- S[y_idx, z_idx, drop = FALSE]
  S_ZY <- t(S_YZ)
  S_ZZ <- S[z_idx, z_idx, drop = FALSE]
  
  # ---- Compute beta and beta0 ----
  beta <- S_YZ %*% solve(S_ZZ)            # (|Y| x |Z|)
  beta0 <- mu_Y - beta %*% mu_Z           # (|Y| x 1)
  
  # ---- Best linear predictor equations: yhat = beta0 + beta z ----
  # We'll format as:
  # y1_hat = beta0_1 + beta_11 z1 + beta_12 z2 + ...
  eq_lines <- character(0)
  for (i in seq_along(Y_names)) {
    rhs <- paste0(round(beta0[i, 1], digits))
    for (j in seq_along(Z_names)) {
      b <- beta[i, j]
      sign_txt <- ifelse(b >= 0, " + ", " - ")
      rhs <- paste0(rhs, sign_txt, round(abs(b), digits), " ", Z_names[j])
    }
    eq_lines <- c(eq_lines, paste0(Y_names[i], "_hat = ", rhs))
  }
  
  # ---- Conditional covariance Sigma_YY|Z ----
  Sigma_YY_given_Z <- S_YY - S_YZ %*% solve(S_ZZ) %*% S_ZY
  
  # ---- Partial correlation between first two Y variables given Z ----
  # rho = sigma12 / sqrt(sigma11*sigma22)
  sigma11 <- Sigma_YY_given_Z[1, 1]
  sigma22 <- Sigma_YY_given_Z[2, 2]
  sigma12 <- Sigma_YY_given_Z[1, 2]
  rho_partial <- as.numeric(sigma12 / sqrt(sigma11 * sigma22))
  
  # ---- Printing (EXAM STYLE) ----
  cat("\n================ STAT 438: Given mu & S =================\n")
  cat("Variables order:\n")
  cat(paste0("  ", paste(var_names, collapse = ", ")), "\n")
  cat("=========================================================\n\n")
  
  if (show_equations) {
    cat("Key formulas:\n")
    cat("  beta   = Sigma_YZ * Sigma_ZZ^{-1}\n")
    cat("  beta0  = mu_Y - beta * mu_Z\n")
    cat("  yhat   = beta0 + beta * z\n")
    cat("  Sigma_YY|Z = Sigma_YY - Sigma_YZ * Sigma_ZZ^{-1} * Sigma_ZY\n")
    cat("  rho_{y1,y2·Z} = sigma12|Z / sqrt(sigma11|Z * sigma22|Z)\n")
    cat("---------------------------------------------------------\n\n")
  }
  
  cat("1) Partitioned blocks:\n")
  cat("\nmu_Y =\n"); print(round(mu_Y, digits))
  cat("\nmu_Z =\n"); print(round(mu_Z, digits))
  
  cat("\nSigma_YY =\n"); print(round(S_YY, digits))
  cat("\nSigma_YZ =\n"); print(round(S_YZ, digits))
  cat("\nSigma_ZZ =\n"); print(round(S_ZZ, digits))
  
  cat("\n2) Regression coefficients:\n")
  cat("\nbeta = Sigma_YZ * Sigma_ZZ^{-1} =\n")
  rownames(beta) <- Y_names
  colnames(beta) <- Z_names
  print(round(beta, digits))
  
  cat("\nbeta0 = mu_Y - beta * mu_Z =\n")
  rownames(beta0) <- Y_names
  print(round(beta0, digits))
  
  cat("\nBest linear predictor (equations):\n")
  cat(paste0("  ", eq_lines, collapse = "\n"), "\n")
  
  cat("\n3) Conditional covariance matrix Sigma_YY|Z:\n")
  rownames(Sigma_YY_given_Z) <- Y_names
  colnames(Sigma_YY_given_Z) <- Y_names
  print(round(Sigma_YY_given_Z, digits))
  
  cat("\n4) Partial correlation between ", Y_names[1], " and ", Y_names[2],
      " eliminating Z:\n", sep = "")
  cat("rho = ", round(rho_partial, digits), "\n", sep = "")
  
  cat("=========================================================\n\n")
  
  invisible(list(
    var_names = var_names,
    y_names = Y_names,
    z_names = Z_names,
    mu_Y = mu_Y,
    mu_Z = mu_Z,
    S_YY = S_YY,
    S_YZ = S_YZ,
    S_ZZ = S_ZZ,
    beta = beta,
    beta0 = beta0,
    equations = eq_lines,
    Sigma_YY_given_Z = Sigma_YY_given_Z,
    rho_partial = rho_partial
  ))
}

# =======================
# Example (YOUR QUESTION)
# =======================
# mu <- c(2.287, 12.6, 0.347, 14.83)
# S <- matrix(c(
#   0.459, 0.254, -0.026, -0.244,
#   0.254, 27.465, -0.589, -0.267,
#  -0.026, -0.589, 0.03,  0.102,
#  -0.244, -0.267, 0.102, 6.854
# ), nrow = 4, byrow = TRUE)
#
# stat438_muS_predictor(mu, S,
#   y_idx = 1:2, z_idx = 3:4,
#   digits = 6,
#   names = c("Y1","Y2","Z1","Z2")
# )
