
# ============================================================
# STAT 438 — ALL EXAM FUNCTIONS (ONE FILE)
# Author: Generated for exam use
# ============================================================

# ---------------- Helper ----------------
`%||%` <- function(x, y) if (is.null(x)) y else x

# ================================
# 1) Best Linear Predictor (mu, S)
# ================================
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


# ================================
# 2) PCA Exam Function
# ================================
# ============================================================
# PCA Exam Function (UPDATED) — works for:
# ✅ Covariance matrix S (like CLEP question)
# ✅ Correlation matrix R (diag=1)
# ✅ Eigenvalues + eigenvectors + PC equations
# ✅ Proportion + cumulative % and #PCs for threshold
# ✅ Corr(Y_i, X_j) for BOTH covariance & correlation PCA
# ✅ Confidence intervals for eigenvalues:
#    - "chisq" (classic chi-square CI)
#    - "bonferroni_z" (the z-form used in your blue solution)
# ============================================================

pca_exam_from_R <- function(R,
                            threshold = 0.80,
                            n = NULL,
                            alpha = 0.05,
                            ci_for = 1:2,
                            ci_method = c("chisq", "bonferroni_z"),
                            corr_pairs = list(),
                            var_names = NULL,
                            digits = 7,
                            show_equations = TRUE) {
  
  options(scipen = 999)
  ci_method <- match.arg(ci_method)
  
  # ---- checks ----
  R <- as.matrix(R)
  if (nrow(R) != ncol(R)) stop("R must be a square matrix (p x p).")
  if (!is.numeric(R)) stop("R must be numeric.")
  if (any(!is.finite(R))) stop("R contains non-finite values.")
  if (is.null(n) && ci_method != "chisq") {
    # bonferroni_z needs n
    # chisq also needs n, but we will just skip CI if n is NULL
  }
  
  p <- nrow(R)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  if (length(var_names) != p) stop("var_names must have length p.")
  
  # Determine if it's correlation-like matrix
  is_corr_like <- all(abs(diag(R) - 1) < 1e-8)
  
  # ---- eigen decomposition ----
  eig <- eigen(R)
  lambda <- as.numeric(eig$values)
  E <- as.matrix(eig$vectors)
  
  # Ensure descending order (eigen() usually is, but keep safe)
  ord <- order(lambda, decreasing = TRUE)
  lambda <- lambda[ord]
  E <- E[, ord, drop = FALSE]
  
  colnames(E) <- paste0("Y", 1:p)
  rownames(E) <- var_names
  
  # ---- proportions ----
  prop <- lambda / sum(lambda)
  cumprop <- cumsum(prop)
  k80 <- which(cumprop >= threshold)[1]
  
  # ---- PC equations helper ----
  pc_equation <- function(i) {
    coefs <- E[, i]
    terms <- paste0(
      ifelse(coefs >= 0, "+ ", "- "),
      format(abs(round(coefs, digits)), nsmall = digits, scientific = FALSE),
      " ",
      var_names
    )
    # remove leading "+ "
    terms[1] <- sub("^\\+\\s+", "", terms[1])
    paste0("Y", i, " = ", paste(terms, collapse = " "))
  }
  
  # ---- Correlations Corr(Y_i, X_j) ----
  # For correlation-matrix PCA: Corr(Y_i, X_j) = sqrt(lambda_i) * e_{j,i}
  # For covariance-matrix PCA:  Corr(Y_i, X_j) = e_{j,i}*sqrt(lambda_i)/sqrt(sigma_jj)
  corr_results <- NULL
  if (length(corr_pairs) > 0) {
    corr_results <- sapply(corr_pairs, function(pair) {
      if (length(pair) != 2) stop("Each element of corr_pairs must be c(i, j).")
      i <- pair[1]  # PC index
      j <- pair[2]  # variable index
      if (i < 1 || i > p || j < 1 || j > p) stop("corr_pairs indices out of range.")
      if (is_corr_like) {
        sqrt(lambda[i]) * E[j, i]
      } else {
        (E[j, i] * sqrt(lambda[i])) / sqrt(R[j, j])
      }
    })
    names(corr_results) <- sapply(corr_pairs, function(pair) {
      paste0("r_{Y", pair[1], ",", var_names[pair[2]], "}")
    })
  }
  
  # ---- Confidence intervals for eigenvalues ----
  ci_table <- NULL
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1 || n <= 1) stop("n must be a single number > 1.")
    ci_for <- as.integer(ci_for)
    ci_for <- ci_for[ci_for >= 1 & ci_for <= p]
    if (length(ci_for) > 0) {
      ci_table <- data.frame(
        eigen = paste0("lambda", ci_for),
        LL = NA_real_,
        lambda_hat = lambda[ci_for],
        UL = NA_real_,
        stringsAsFactors = FALSE
      )
      
      if (ci_method == "chisq") {
        # Classic CI: ( (n-1)*lambda / chi^2_{1-a/2}, (n-1)*lambda / chi^2_{a/2} )
        df <- n - 1
        for (idx in seq_along(ci_for)) {
          i <- ci_for[idx]
          ci_table$LL[idx] <- (df * lambda[i]) / qchisq(1 - alpha/2, df)
          ci_table$UL[idx] <- (df * lambda[i]) / qchisq(alpha/2, df)
        }
      }
      
      if (ci_method == "bonferroni_z") {
        # Matches your blue sheet:
        # lambda_hat / (1 + z(alpha/(2m))*sqrt(2/n))  <= lambda <=  lambda_hat/(1 - z(..)*sqrt(2/n))
        m <- p
        zcrit <- qnorm(1 - alpha/(2*m))
        k <- zcrit * sqrt(2/n)
        for (idx in seq_along(ci_for)) {
          i <- ci_for[idx]
          ci_table$LL[idx] <- lambda[i] / (1 + k)
          ci_table$UL[idx] <- lambda[i] / (1 - k)
        }
      }
    }
  }
  
  # ---- printing ----
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  
  cat("\n==================== PCA (Exam Output) ====================\n")
  cat(sprintf("p = %d variables\n", p))
  cat(sprintf("Matrix type detected: %s\n", ifelse(is_corr_like, "Correlation (diag=1)", "Covariance")))
  if (!is.null(n)) cat(sprintf("n = %d, alpha = %.4f\n", n, alpha))
  cat("===========================================================\n\n")
  
  if (show_equations) {
    cat("Key facts:\n")
    cat("  Eigen decomposition: R e_i = lambda_i e_i\n")
    cat("  PC_i = e_i' X\n")
    cat("  Proportion = lambda_i / sum(lambda)\n")
    cat("  Corr(Y_i, X_j):\n")
    cat("    - if R is correlation: sqrt(lambda_i) * e_{j,i}\n")
    cat("    - if R is covariance:  e_{j,i}*sqrt(lambda_i)/sqrt(R_jj)\n")
    cat("  CI methods:\n")
    cat("    - chisq: (n-1)lambda / chi^2\n")
    cat("    - bonferroni_z: lambda /(1 ± z(alpha/(2p))*sqrt(2/n))\n")
    cat("-----------------------------------------------------------\n\n")
  }
  
  cat("1) Eigenvalues (lambda):\n")
  print(setNames(as.numeric(fmt(lambda)), paste0("lambda", 1:p)))
  
  cat("\n2) Eigenvectors (columns are PCs Y1..Yp):\n")
  print(apply(E, 2, fmt))
  
  cat("\n3) Proportion of variance explained:\n")
  prop_out <- rbind(
    Proportion = as.numeric(fmt(prop)),
    Cumulative = as.numeric(fmt(cumprop))
  )
  colnames(prop_out) <- paste0("Y", 1:p)
  print(prop_out)
  
  cat(sprintf("\n4) #PCs to reach %.0f%% cumulative variance: %d PCs (cum = %s)\n",
              100*threshold, k80, fmt(cumprop[k80])))
  
  cat("\n5) Principal component equations:\n")
  for (i in 1:p) cat("  ", pc_equation(i), "\n", sep = "")
  
  if (!is.null(corr_results)) {
    cat("\n6) Requested correlations Corr(Y_i, X_j):\n")
    print(setNames(as.numeric(fmt(corr_results)), names(corr_results)))
  } else if (length(corr_pairs) > 0) {
    cat("\n6) Requested correlations were not computed (check corr_pairs).\n")
  }
  
  if (!is.null(ci_table)) {
    cat("\n7) Confidence intervals for selected eigenvalues:\n")
    cat(sprintf("   Method: %s\n", ci_method))
    ci_print <- ci_table
    ci_print$LL <- as.numeric(fmt(ci_print$LL))
    ci_print$lambda_hat <- as.numeric(fmt(ci_print$lambda_hat))
    ci_print$UL <- as.numeric(fmt(ci_print$UL))
    print(ci_print, row.names = FALSE)
  } else if (!is.null(n) && length(ci_for) == 0) {
    cat("\n7) No eigenvalues selected for CI (ci_for empty after filtering).\n")
  } else {
    cat("\n7) CI not computed (provide n to compute CI).\n")
  }
  
  cat("===========================================================\n\n")
  
  invisible(list(
    p = p,
    n = n,
    matrix_type = ifelse(is_corr_like, "correlation", "covariance"),
    lambda = lambda,
    E = E,
    prop = prop,
    cumprop = cumprop,
    k_threshold = k80,
    corr = corr_results,
    ci = ci_table
  ))
}

# ============================================================
# Example: CLEP covariance-matrix question
# ============================================================
# n <- 87
# S <- matrix(c(
#   5691, 600, 217,
#   600,  126, 24,
#   217,   24, 23
# ), 3, 3, byrow=TRUE)
#
# pca_exam_from_R(
#   R = S,
#   n = 87,
#   alpha = 0.01,                 # 99% => alpha=0.01
#   ci_method = "bonferroni_z",    # match the blue solution style
#   ci_for = 1:3,
#   corr_pairs = list(c(1,1), c(1,2), c(1,3)),  # Corr(Y1, X1), Corr(Y1, X2), Corr(Y1, X3)
#   var_names = c("X1","X2","X3"),
#   digits = 7
# )


# ================================
# 3) Hotelling T2 (Contrasts)
# ================================
hotelling_T2_contrast <- function(xbar, S, n, C, alpha = 0.05, digits = 6, show_equations = TRUE) {
  
  # ----------------------------
  # Basic checks
  # ----------------------------
  xbar <- as.matrix(xbar)
  if (ncol(xbar) != 1) xbar <- matrix(xbar, ncol = 1)
  
  S <- as.matrix(S)
  C <- as.matrix(C)
  
  p <- nrow(S)
  q <- qr(C)$rank
  
  if (ncol(C) != p) stop("C must have p columns.")
  if (n <= q) stop("n must be larger than the number of contrasts q.")
  if (nrow(S) != ncol(S)) stop("S must be a square matrix.")
  
  # ----------------------------
  # Core calculations
  # ----------------------------
  delta_hat <- C %*% xbar
  S_delta   <- C %*% S %*% t(C)
  
  # Safe inverse (usually not needed in exams, but prevents crashes)
  if (qr(S_delta)$rank < nrow(S_delta)) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("S_delta is singular. Install MASS (install.packages('MASS')) or use full-rank contrasts.")
    }
    S_delta_i <- MASS::ginv(S_delta)
  } else {
    S_delta_i <- solve(S_delta)
  }
  
  T2 <- as.numeric(n * t(delta_hat) %*% S_delta_i %*% delta_hat)
  
  # ----------------------------
  # Sampling distribution + critical value
  # ----------------------------
  F_stat <- ((n - q) / (q * (n - 1))) * T2
  F_crit <- qf(1 - alpha, q, n - q)
  T2_crit <- (q * (n - 1) / (n - q)) * F_crit
  
  decision <- ifelse(T2 > T2_crit, "Reject H0", "Do NOT reject H0")
  
  # ----------------------------
  # (A) Hotelling simultaneous CIs for contrasts (components of C mu)
  # delta_i ± sqrt(T2_crit * S_delta[ii] / n)
  # ----------------------------
  se <- sqrt(diag(S_delta) / n)
  hw_hot <- sqrt(T2_crit) * se
  
  CI_hotelling <- cbind(
    Lower    = as.numeric(delta_hat) - hw_hot,
    Estimate = as.numeric(delta_hat),
    Upper    = as.numeric(delta_hat) + hw_hot
  )
  
  # ----------------------------
  # (B) Bonferroni simultaneous CIs for contrasts
  # delta_i ± t_{n-1, 1 - alpha/(2m)} * sqrt(S_delta[ii] / n)
  # where m = number of intervals = number of contrasts (rows of C)
  # ----------------------------
  m <- nrow(C)
  tcrit_bonf <- qt(1 - alpha / (2 * m), df = n - 1)
  hw_bonf <- tcrit_bonf * se
  
  CI_bonferroni <- cbind(
    Lower    = as.numeric(delta_hat) - hw_bonf,
    Estimate = as.numeric(delta_hat),
    Upper    = as.numeric(delta_hat) + hw_bonf
  )
  
  # ----------------------------
  # Print (EXAM FORMAT)
  # ----------------------------
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  
  cat("\n================ Hotelling's T^2 (Contrasts) ================\n")
  cat("H0 : C * mu = 0\n")
  cat(sprintf("n = %d, p = %d, q = %d, alpha = %.4f\n", n, p, q, alpha))
  cat("=============================================================\n\n")
  
  if (show_equations) {
    cat("Equations:\n")
    cat("  delta_hat = C * xbar\n")
    cat("  S_delta   = C * S * C'\n")
    cat("  T^2       = n * delta_hat' * (S_delta)^{-1} * delta_hat\n")
    cat("  F         = ((n-q)/(q(n-1))) * T^2  ~  F_{q, n-q}\n")
    cat("  T^2_crit  = (q(n-1)/(n-q)) * F_{q,n-q; 1-alpha}\n")
    cat("  Hotelling CI_i : delta_i ± sqrt(T^2_crit * S_delta[ii] / n)\n")
    cat("  Bonferroni CI_i: delta_i ± t_{n-1, 1-alpha/(2m)} * sqrt(S_delta[ii] / n),  m = #contrasts\n")
    cat("-------------------------------------------------------------\n\n")
  }
  
  cat("C =\n")
  print(C)
  
  cat("\n1) delta_hat = C * xbar:\n")
  print(matrix(fmt(delta_hat), ncol = 1))
  
  cat("\n2) S_delta = C * S * C':\n")
  print(matrix(fmt(S_delta), nrow = nrow(S_delta)))
  
  cat("\n3) T^2 statistic:\n")
  cat("T^2 =", fmt(T2), "\n")
  
  cat("\n4) Sampling distribution:\n")
  cat("((n - q)/(q(n - 1))) * T^2  ~  F_{q, n - q}\n")
  cat("q =", q, ", n =", n, "\n")
  
  cat("\n5) Critical value:\n")
  cat("T^2_crit =", fmt(T2_crit), "\n")
  
  cat("\n6) Decision:\n")
  cat(decision, "\n")
  
  cat("\n7A) (1-alpha) Hotelling simultaneous confidence intervals:\n")
  print(apply(CI_hotelling, 2, fmt))
  
  cat("\n7B) (1-alpha) Bonferroni simultaneous confidence intervals:\n")
  cat(sprintf("t_crit(Bonferroni) = %s  with df = %d and alpha/(2m) = %.6f\n",
              fmt(tcrit_bonf), n - 1, alpha / (2 * m)))
  print(apply(CI_bonferroni, 2, fmt))
  
  cat("=============================================================\n\n")
  
  # ----------------------------
  # Return everything
  # ----------------------------
  invisible(list(
    delta_hat = delta_hat,
    S_delta = S_delta,
    T2 = T2,
    F_stat = F_stat,
    T2_crit = T2_crit,
    decision = decision,
    CI_hotelling = CI_hotelling,
    CI_bonferroni = CI_bonferroni,
    tcrit_bonferroni = tcrit_bonf,
    q = q,
    m = m
  ))
}


# ================================
# 4) One-way MANOVA (Summary)
# ================================
# ============================================================
# STAT 438 — One-way MANOVA from summary stats (exam-style)
# Solves questions like your Q3:
# 1) Box's M test for equality of covariance matrices
# 2) Compute grand mean, B, W, T
# 3) Wilks' Lambda, F-approx + critical value + decision
#
# Inputs:
# - xbars: list of mean vectors (each length p)
# - Slist: list of covariance matrices (p x p) for each group
# - n:     vector of sample sizes for each group
# - alpha: significance level
#
# Notes:
# - Works with any p >= 2 and any g >= 2
# - Box's M uses a standard chi-square approximation
# - Wilks F-approx uses the common Rao's approximation (general)
#   and also prints the special-case sqrt(Lambda) form when p=2
# ============================================================

stat438_oneway_manova_summary <- function(xbars, Slist, n, alpha = 0.05, digits = 6,
                                          show_equations = TRUE) {
  
  # ---- checks ----
  if (!is.list(xbars) || !is.list(Slist)) stop("xbars and Slist must be lists.")
  g <- length(xbars)
  if (length(Slist) != g) stop("xbars and Slist must have the same length.")
  if (length(n) != g) stop("n must have length g (one sample size per group).")
  if (g < 2) stop("Need at least 2 groups.")
  
  # force matrices
  xbars <- lapply(xbars, function(x) matrix(as.numeric(x), ncol = 1))
  Slist <- lapply(Slist, function(S) as.matrix(S))
  
  p <- nrow(Slist[[1]])
  for (i in 1:g) {
    if (nrow(Slist[[i]]) != p || ncol(Slist[[i]]) != p) stop("All S matrices must be p x p.")
    if (nrow(xbars[[i]]) != p) stop("All mean vectors must have length p.")
    if (n[i] <= 1) stop("Each group size must be > 1.")
  }
  
  N <- sum(n)
  df_within <- sum(n - 1)
  
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  
  # ============================================================
  # Part A) Box's M test
  # ============================================================
  
  # pooled covariance
  Sp_num <- Reduce(`+`, lapply(1:g, function(i) (n[i] - 1) * Slist[[i]]))
  Sp <- Sp_num / df_within
  
  # determinants (use determinant() for numerical stability)
  logdet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)
  
  M <- df_within * logdet(Sp) - sum(sapply(1:g, function(i) (n[i]-1) * logdet(Slist[[i]])))
  
  # correction factor (standard Box M chi-square approx)
  # c = [ (2p^2 + 3p - 1) / (6(p+1)(g-1)) ] * [ sum 1/(n_i-1) - 1/df_within ]
  U_raw <- sum(1/(n - 1)) - 1/df_within
  c_fac <- ((2*p^2 + 3*p - 1) / (6*(p + 1)*(g - 1))) * U_raw
  
  D <- (1 - c_fac) * M
  df_box <- (g - 1) * p * (p + 1) / 2
  chi_crit <- qchisq(1 - alpha, df_box)
  box_decision <- ifelse(D > chi_crit, "Reject H0 (covariances NOT equal)", "Do NOT reject H0 (covariances can be treated equal)")
  
  # ============================================================
  # Part B) MANOVA SSCP matrices: grand mean, B, W, T
  # ============================================================
  
  # grand mean
  xbar <- Reduce(`+`, lapply(1:g, function(i) n[i] * xbars[[i]])) / N
  
  # W
  W <- Sp_num
  
  # B
  B <- Reduce(`+`, lapply(1:g, function(i) {
    d <- xbars[[i]] - xbar
    n[i] * (d %*% t(d))
  }))
  
  T <- B + W
  
  # ============================================================
  # Part C) Wilks' Lambda and F approximation
  # ============================================================
  
  Lambda <- det(W) / det(T)
  
  # General Rao's F approximation (works for general p,g)
  # Let s = sqrt( (p^2 (g-1)^2 - 4) / (p^2 + (g-1)^2 - 5) ) if denom>0 else 1
  # v = N - 1 - (p + g)/2
  # df1 = p(g-1)
  # df2 = s*v - (df1 - 2)/2
  # F = ((1 - Lambda^(1/s)) / (Lambda^(1/s))) * (df2/df1)
  df1 <- p * (g - 1)
  
  denom <- (p^2 + (g - 1)^2 - 5)
  if (denom > 0) {
    s <- sqrt((p^2 * (g - 1)^2 - 4) / denom)
  } else {
    s <- 1
  }
  
  v <- N - 1 - (p + g)/2
  # avoid negative/invalid df2 in tiny samples
  df2 <- max(1e-9, s * v - (df1 - 2)/2)
  
  Lambda_pow <- Lambda^(1/s)
  F_stat <- ((1 - Lambda_pow) / Lambda_pow) * (df2 / df1)
  F_crit <- qf(1 - alpha, df1, df2)
  wilks_decision <- ifelse(F_stat > F_crit, "Reject H0 (group means differ)", "Do NOT reject H0")
  
  # Special-case form often used in your course when p=2 (matches many model answers)
  special_F <- NULL
  special_df1 <- NULL
  special_df2 <- NULL
  special_Fcrit <- NULL
  if (p == 2) {
    # F = ((N-g-1)/(g-1)) * ((1 - sqrt(Lambda))/sqrt(Lambda))
    special_F <- ((N - g - 1) / (g - 1)) * ((1 - sqrt(Lambda)) / sqrt(Lambda))
    special_df1 <- 2 * (g - 1)
    special_df2 <- 2 * (N - g - 1)
    if (special_df2 > 0) special_Fcrit <- qf(1 - alpha, special_df1, special_df2)
  }
  
  # ============================================================
  # PRINT (EXAM STYLE)
  # ============================================================
  
  cat("\n================== One-way MANOVA (Summary) ==================\n")
  cat(sprintf("g = %d groups, p = %d responses, N = %d, alpha = %.3f\n", g, p, N, alpha))
  cat("==============================================================\n\n")
  
  if (show_equations) {
    cat("Formulas:\n")
    cat("  Sp = sum (n_i-1) S_i / sum(n_i-1)\n")
    cat("  M  = [sum(n_i-1)] ln|Sp| - sum (n_i-1) ln|S_i|\n")
    cat("  c  = [(2p^2+3p-1)/(6(p+1)(g-1))] * [sum 1/(n_i-1) - 1/sum(n_i-1)]\n")
    cat("  D  = (1-c) M  ~  Chi-square_{(g-1)p(p+1)/2}\n")
    cat("  xbar = (1/N) sum n_i xbar_i\n")
    cat("  W = sum (n_i-1) S_i\n")
    cat("  B = sum n_i (xbar_i - xbar)(xbar_i - xbar)'\n")
    cat("  T = B + W\n")
    cat("  Lambda = |W|/|T|\n")
    cat("--------------------------------------------------------------\n\n")
  }
  
  # 1) Box M
  cat("1) Box's M test for equal covariance matrices\n\n")
  cat("Pooled covariance Sp =\n")
  print(matrix(fmt(Sp), nrow = p))
  
  cat("\nU_raw = sum 1/(n_i-1) - 1/sum(n_i-1) = ", fmt(U_raw), "\n", sep="")
  cat("c (correction factor) = ", fmt(c_fac), "\n", sep="")
  cat("M = ", fmt(M), "\n", sep="")
  cat("D = (1-c)M = ", fmt(D), "\n", sep="")
  cat("Critical value (Chi-square) = ", fmt(chi_crit), " with df = ", df_box, "\n", sep="")
  cat("Decision: ", box_decision, "\n", sep="")
  
  # 2) SSCP
  cat("\n2) One-way MANOVA SSCP matrices\n\n")
  cat("Grand mean xbar =\n")
  print(matrix(fmt(xbar), ncol = 1))
  
  cat("\nBetween SSCP (B) =\n")
  print(matrix(fmt(B), nrow = p))
  
  cat("\nWithin SSCP (W) =\n")
  print(matrix(fmt(W), nrow = p))
  
  cat("\nTotal SSCP (T=B+W) =\n")
  print(matrix(fmt(T), nrow = p))
  
  # 3) Wilks
  cat("\n3) Wilks' Lambda test\n\n")
  cat("Lambda = |W|/|T| = ", fmt(Lambda), "\n", sep="")
  cat("Rao F-approx: F = ", fmt(F_stat), "  with df = (", fmt(df1), ", ", fmt(df2), ")\n", sep="")
  cat("Critical value = ", fmt(F_crit), "\n", sep="")
  cat("Decision: ", wilks_decision, "\n", sep="")
  
  if (!is.null(special_F) && !is.null(special_Fcrit)) {
    cat("\n(Extra) Course special-case when p=2:\n")
    cat("F = ((N-g-1)/(g-1)) * ((1 - sqrt(Lambda))/sqrt(Lambda))\n")
    cat("F_special = ", fmt(special_F), " with df = (", special_df1, ", ", special_df2, ")\n", sep="")
    cat("F_crit_special = ", fmt(special_Fcrit), "\n", sep="")
  }
  
  cat("==============================================================\n\n")
  
  invisible(list(
    g = g, p = p, n = n, N = N,
    Sp = Sp, M = M, U_raw = U_raw, c = c_fac, D = D, df_box = df_box,
    chi_crit = chi_crit, box_decision = box_decision,
    xbar = xbar, B = B, W = W, T = T,
    Lambda = Lambda,
    F_stat = F_stat, df1 = df1, df2 = df2, F_crit = F_crit, wilks_decision = wilks_decision,
    F_special = special_F, F_special_df1 = special_df1, F_special_df2 = special_df2, F_special_crit = special_Fcrit
  ))
}

# =======================
# Example (YOUR QUESTION)
# =======================
# x1 <- c(96, 92.3)
# x2 <- c(90.25, 88.25)
# x3 <- c(77.33, 73.33)
#
# S1 <- matrix(c(1, 0.5, 0.5, 0.33), 2, 2, byrow=TRUE)
# S2 <- matrix(c(30.92, -8.42, -8.42, 4.91), 2, 2, byrow=TRUE)
# S3 <- matrix(c(4.33, 1.83, 1.83, 2.33), 2, 2, byrow=TRUE)
#
# stat438_oneway_manova_summary(
#   xbars = list(x1, x2, x3),
#   Slist = list(S1, S2, S3),
#   n = c(3,3,3),
#   alpha = 0.05,
#   digits = 6
# )


# ================================
# 5) STAT 438 Q4 Regression
# ================================
# ============================================================
# STAT 438 — Q4 Exam Function (FULL WORKING SCRIPT)
# Multivariate Linear Regression (2 responses: y1,y2)
# Predictors: z1 required, z2 optional
#
# Prints exam-ready steps:
# 1) OLS B_hat
# 2) Fitted values + residuals
# 3) SSCP decomposition check
# 4) MLE regression using mean vector + covariance matrix
# 5) Partial correlation r_{y1,y2 · Z}
# 6) MLE of Sigma:  Sigma_hat = (1/n) E'E
# 7) 95% CI for mean response E(Y_i|z0)
# 8) 95% PI for response Y_i at z0
#
# IMPORTANT:
# Your course PDF may use df2 = n - r - m - 1 (common) or df2 = n - r - m.
# Use df2_mode to match the PDF:
#   df2_mode = "n-r-m-1"  (often matches the provided solution)
#   df2_mode = "n-r-m"
# ============================================================

stat438_q4 <- function(z1, y1, y2,
                       z2 = NULL,
                       digits = 4,
                       alpha = 0.05,
                       z0 = NULL,                 # e.g. c(1,6,10) if intercept+z1+z2
                       response_index = 1,        # 1 for Y1, 2 for Y2
                       do_intervals = TRUE,
                       df2_mode = c("n-r-m", "n-r-m")) {
  
  options(scipen = 999)
  df2_mode <- match.arg(df2_mode)
  
  # ----------------------------
  # Input checks
  # ----------------------------
  if (length(z1) != length(y1) || length(y1) != length(y2)) {
    stop("z1, y1, and y2 must have the same length.")
  }
  if (!is.null(z2) && length(z2) != length(z1)) {
    stop("z2 must have the same length as z1.")
  }
  if (!(response_index %in% c(1, 2))) stop("response_index must be 1 (Y1) or 2 (Y2).")
  
  n <- length(z1)
  Y <- cbind(y1, y2)
  m <- ncol(Y)  # number of responses (m = 2)
  
  # Predictor matrix (no intercept yet)
  Zvars <- if (is.null(z2)) cbind(z1) else cbind(z1, z2)
  colnames(Zvars) <- if (is.null(z2)) "z1" else c("z1", "z2")
  
  # Design matrix with intercept
  Z <- cbind(1, Zvars)
  colnames(Z)[1] <- "Intercept"
  
  r <- ncol(Z) - 1  # number of predictors (excluding intercept)
  
  # ----------------------------
  # 1) OLS Estimates: B_hat = (Z'Z)^(-1) Z'Y
  # ----------------------------
  ZtZ_inv <- solve(t(Z) %*% Z)
  B_hat <- ZtZ_inv %*% t(Z) %*% Y
  rownames(B_hat) <- colnames(Z)
  colnames(B_hat) <- c("y1", "y2")
  
  cat("1) Least squares estimates of the parameters, B-hat\n\n")
  cat("B_hat = (Z'Z)^(-1) Z'Y =\n")
  print(round(B_hat, digits))
  
  # ----------------------------
  # 2) Fitted values & residuals
  # ----------------------------
  Y_hat <- Z %*% B_hat
  E_hat <- Y - Y_hat
  
  cat("\n2) Fitted values and residuals\n\n")
  cat("Y_hat = Z B_hat =\n")
  print(round(Y_hat, digits))
  
  cat("\nE_hat = Y - Y_hat =\n")
  print(round(E_hat, digits))
  
  # ----------------------------
  # 3) SSCP decomposition check
  # ----------------------------
  YY   <- t(Y) %*% Y
  YHYH <- t(Y_hat) %*% Y_hat
  EE   <- t(E_hat) %*% E_hat
  
  cat("\n3) Verify the SSCP decomposition\n\n")
  cat("Y'Y =\n")
  print(round(YY, digits))
  
  cat("\nY_hat'Y_hat + E_hat'E_hat =\n")
  print(round(YHYH + EE, digits))
  
  # ----------------------------
  # 4) MLE regression using mean vector and covariance matrix
  # ----------------------------
  joint <- cbind(y1, y2, Zvars)
  mu <- colMeans(joint)
  S  <- stats::cov(joint)
  
  mu_y <- matrix(mu[1:2], ncol = 1)
  mu_z <- matrix(mu[3:(2 + r)], ncol = 1)
  
  S_YY <- S[1:2, 1:2, drop = FALSE]
  S_YZ <- S[1:2, 3:(2 + r), drop = FALSE]
  S_ZZ <- S[3:(2 + r), 3:(2 + r), drop = FALSE]
  
  B_MLE  <- S_YZ %*% solve(S_ZZ)
  beta0_MLE <- mu_y - B_MLE %*% mu_z
  
  rownames(B_MLE) <- c("y1", "y2")
  colnames(B_MLE) <- colnames(Zvars)
  
  cat("\n4) MLE regression using mean vector and covariance matrix\n\n")
  cat("mu =\n")
  print(round(mu, digits))
  
  cat("\nS =\n")
  print(round(S, digits))
  
  cat("\nB_MLE = S_YZ S_ZZ^(-1) =\n")
  print(round(B_MLE, digits))
  
  cat("\nbeta0_MLE = mu_Y - B_MLE mu_Z =\n")
  print(round(beta0_MLE, digits))
  
  cat("\nEstimated regression equations:\n")
  if (is.null(z2)) {
    cat(paste0("y1_hat = ", round(beta0_MLE[1], digits),
               " + ", round(B_MLE[1, "z1"], digits), " z1\n"))
    cat(paste0("y2_hat = ", round(beta0_MLE[2], digits),
               " + ", round(B_MLE[2, "z1"], digits), " z1\n"))
  } else {
    cat(paste0("y1_hat = ", round(beta0_MLE[1], digits),
               " + ", round(B_MLE[1, "z1"], digits), " z1",
               " + ", round(B_MLE[1, "z2"], digits), " z2\n"))
    cat(paste0("y2_hat = ", round(beta0_MLE[2], digits),
               " + ", round(B_MLE[2, "z1"], digits), " z1",
               " + ", round(B_MLE[2, "z2"], digits), " z2\n"))
  }
  
  # ----------------------------
  # 5) Partial correlation r_{y1,y2 · Z}
  # ----------------------------
  Sigma_YY_given_Z <- S_YY - S_YZ %*% solve(S_ZZ) %*% t(S_YZ)
  r_partial <- Sigma_YY_given_Z[1, 2] / sqrt(Sigma_YY_given_Z[1, 1] * Sigma_YY_given_Z[2, 2])
  
  cat("\n5) Partial correlation coefficient\n\n")
  cat("Sigma_YY|Z = S_YY - S_YZ S_ZZ^(-1) S_ZY =\n")
  print(round(Sigma_YY_given_Z, digits))
  cat("\nr_{y1,y2 · ", if (is.null(z2)) "z1" else "(z1,z2)", "} =\n", sep = "")
  cat(round(r_partial, 2), "\n")
  
  # ----------------------------
  # 6) MLE of Sigma: Sigma_hat = (1/n) E'E
  # ----------------------------
  Sigma_hat <- (1 / n) * (t(E_hat) %*% E_hat)
  
  cat("\n6) Find the maximum likelihood estimate of Sigma\n\n")
  cat("Sigma_hat = (1/n) * E' E =\n")
  print(round(Sigma_hat, digits))
  
  # ----------------------------
  # 7) 95% CI for mean response and 8) 95% PI for response at z0
  # ----------------------------
  CI_mean <- NULL
  PI_pred <- NULL
  
  if (do_intervals) {
    
    if (is.null(z0)) {
      stop("Provide z0 including intercept. Example: z0 = c(1,6,10).")
    }
    
    z0 <- matrix(as.numeric(z0), ncol = 1)
    if (length(z0) != ncol(Z)) {
      stop(paste0("z0 must have length ", ncol(Z), " (Intercept + predictors). Example: c(1,6,10)."))
    }
    
    # Degrees of freedom for F
    df1 <- m
    df2 <- if (df2_mode == "n-r-m") (n - r - m) else (n - r - m - 1)
    if (df2 <= 0) stop("df2 must be positive. Check n, r, m, or df2_mode.")
    
    Fcrit <- qf(1 - alpha, df1, df2)
    
    denom <- if (df2_mode == "n-r-m") (n - r - m) else (n - r - m - 1)
    cfac <- (m * (n - r - 1) / denom) * Fcrit
    
    # Scale Sigma_hat as in the slide: (n/(n-r-1)) * sigma_hat_ii
    Sigma_tilde <- (n / (n - r - 1)) * Sigma_hat
    
    i <- response_index
    sigma_ii <- Sigma_tilde[i, i]
    
    # mean at z0: z0' b_i
    b_i <- matrix(B_hat[, i], ncol = 1)
    mean_hat <- as.numeric(t(z0) %*% b_i)
    
    h <- as.numeric(t(z0) %*% ZtZ_inv %*% z0)
    
    hw_CI <- sqrt(cfac * h * sigma_ii)
    hw_PI <- sqrt(cfac * (1 + h) * sigma_ii)
    
    CI_mean <- c(mean_hat - hw_CI, mean_hat + hw_CI)
    PI_pred <- c(mean_hat - hw_PI, mean_hat + hw_PI)
    
    cat("\n7) Find a 95% confidence interval for the mean response E(Y", i, ") at z0\n\n", sep = "")
    cat("z0 =\n"); print(round(z0, digits))
    cat("df2_mode = ", df2_mode, "  =>  F ~ F_{", df1, ",", df2, "}\n", sep = "")
    cat("E(Y", i, "|z0) = z0' b_hat(i) = ", round(mean_hat, digits), "\n", sep = "")
    cat("95% CI = (", round(CI_mean[1], digits), ", ", round(CI_mean[2], digits), ")\n", sep = "")
    
    cat("\n8) Find a 95% prediction interval for the response Y", i, " at z0\n\n", sep = "")
    cat("95% PI = (", round(PI_pred[1], digits), ", ", round(PI_pred[2], digits), ")\n", sep = "")
  }
  
  invisible(list(
    B_hat = B_hat,
    Y_hat = Y_hat,
    E_hat = E_hat,
    Sigma_hat = Sigma_hat,
    SSCP_YtY = YY,
    SSCP_YhatYhat = YHYH,
    SSCP_EtE = EE,
    mu = mu,
    S = S,
    B_MLE = B_MLE,
    beta0_MLE = beta0_MLE,
    Sigma_YY_given_Z = Sigma_YY_given_Z,
    partial_corr = r_partial,
    z0 = z0,
    CI_mean = CI_mean,
    PI_pred = PI_pred,
    df2_mode = df2_mode
  ))
}

# =======================
# Example (YOUR TABLE)
# =======================
# z1 <- c(2,5,7,8,9,12,11,10,15)
# z2 <- c(-4,-1,0,3,7,6,12,17,15)
# y1 <- c(0.3,0.5,1.2,3,4.5,6.1,7,9,7.5)
# y2 <- c(10,11,15,14,17,13,19,16,18)
#
# stat438_q4(z1=z1, z2=z2, y1=y1, y2=y2, z0=c(1,6,10), response_index=1,
#            digits=5, df2_mode="n-r-m-1")


# ================================
# 6) Paired Hotelling T2
# ================================
# ============================================================
# Paired Hotelling's T^2 (Before vs After, multivariate)
# Exam-style output: dbar, Sd, Sd^-1, T2, F, crit, decision,
# and BOTH simultaneous CIs:
#   (A) Hotelling T^2 simultaneous CIs
#   (B) Bonferroni simultaneous CIs
# ============================================================

# Small helper: "x %||% y" means if x is NULL use y
`%||%` <- function(x, y) if (is.null(x)) y else x

paired_hotelling_T2 <- function(before, after,
                                alpha = 0.05,
                                diff = c("before-after", "after-before"),
                                digits = 6,
                                show_equations = TRUE) {
  
  diff <- match.arg(diff)
  
  # ---- checks ----
  if (!is.matrix(before)) before <- as.matrix(before)
  if (!is.matrix(after))  after  <- as.matrix(after)
  
  if (nrow(before) != nrow(after) || ncol(before) != ncol(after)) {
    stop("before and after must have the same dimensions (n x p).")
  }
  
  n <- nrow(before)
  p <- ncol(before)
  
  if (n <= p) stop("Need n > p for the Hotelling T^2 to be valid (so F conversion works).")
  
  # ---- differences ----
  D <- if (diff == "before-after") before - after else after - before
  
  # ---- mean difference vector and covariance of differences ----
  dbar <- colMeans(D)
  Sd   <- stats::cov(D)          # sample covariance uses (n-1) by default
  Sd_inv <- solve(Sd)
  
  # ---- Hotelling T^2 test ----
  T2 <- as.numeric(n * t(dbar) %*% Sd_inv %*% dbar)
  
  # ---- Convert to F ----
  F_stat <- ((n - p) / (p * (n - 1))) * T2
  df1 <- p
  df2 <- n - p
  
  F_crit <- stats::qf(1 - alpha, df1, df2)
  T2_crit <- (p * (n - 1) / (n - p)) * F_crit
  
  decision <- if (T2 > T2_crit) "Reject H0" else "Do NOT reject H0"
  
  # ---- (A) Hotelling simultaneous CIs ----
  # delta_j: dbar_j ± sqrt(T2_crit * s_jj / n)
  se <- sqrt(diag(Sd) / n)
  hw_hotelling <- sqrt(T2_crit) * se
  
  CI_hotelling <- cbind(
    Lower = dbar - hw_hotelling,
    Estimate = dbar,
    Upper = dbar + hw_hotelling
  )
  
  # ---- (B) Bonferroni simultaneous CIs ----
  # delta_j: dbar_j ± t_{n-1, 1 - alpha/(2p)} * sqrt(s_jj / n)
  tcrit_bonf <- stats::qt(1 - alpha / (2 * p), df = n - 1)
  hw_bonf <- tcrit_bonf * se
  
  CI_bonferroni <- cbind(
    Lower = dbar - hw_bonf,
    Estimate = dbar,
    Upper = dbar + hw_bonf
  )
  
  # ---- printing helpers ----
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  vnames <- colnames(before) %||% paste0("V", 1:p)
  
  cat("\n============================================================\n")
  cat("Paired Hotelling's T^2 Test (Before vs After)\n")
  cat("============================================================\n")
  cat(sprintf("n = %d, p = %d, alpha = %.4f\n", n, p, alpha))
  cat(sprintf("Difference defined as: %s\n\n", diff))
  
  if (show_equations) {
    cat("Equations:\n")
    cat("  dbar = (1/n) * sum d_i\n")
    cat("  Sd   = (1/(n-1)) * sum (d_i - dbar)(d_i - dbar)'\n")
    cat("  T^2  = n * dbar' * Sd^{-1} * dbar\n")
    cat("  F    = ((n-p)/(p(n-1))) * T^2  ~  F_{p, n-p}\n")
    cat("  T^2_crit = (p(n-1)/(n-p)) * F_{p,n-p; 1-alpha}\n")
    cat("  Hotelling CI_j: dbar_j ± sqrt(T^2_crit * s_jj / n)\n")
    cat("  Bonferroni CI_j: dbar_j ± t_{n-1, 1-alpha/(2p)} * sqrt(s_jj / n)\n")
    cat("------------------------------------------------------------\n\n")
  }
  
  cat("1) Mean difference vector (dbar):\n")
  print(setNames(as.numeric(fmt(dbar)), vnames))
  
  cat("\n2) Covariance matrix of differences (Sd):\n")
  print(apply(Sd, 2, fmt) |>
          matrix(nrow = p, dimnames = list(vnames, vnames)))
  
  cat("\n3) Inverse covariance (Sd^{-1}):\n")
  print(apply(Sd_inv, 2, fmt) |>
          matrix(nrow = p, dimnames = list(vnames, vnames)))
  
  cat("\n4) Test statistics:\n")
  cat(sprintf("  T^2      = %s\n", fmt(T2)))
  cat(sprintf("  F        = %s  with df = (%d, %d)\n", fmt(F_stat), df1, df2))
  cat(sprintf("  F_crit   = %s\n", fmt(F_crit)))
  cat(sprintf("  T^2_crit = %s\n", fmt(T2_crit)))
  
  cat("\n5) Decision:\n")
  cat(sprintf("  %s  (since T^2 %s T^2_crit)\n",
              decision, ifelse(T2 > T2_crit, ">", "<=")))
  
  cat("\n6A) (1-alpha) Hotelling simultaneous CIs for delta:\n")
  CIh_out <- apply(CI_hotelling, 2, fmt)
  rownames(CIh_out) <- vnames
  print(CIh_out)
  
  cat("\n6B) (1-alpha) Bonferroni simultaneous CIs for delta:\n")
  cat(sprintf("    t_crit(Bonferroni) = %s  with df = %d and alpha/(2p) = %.6f\n",
              fmt(tcrit_bonf), n - 1, alpha / (2 * p)))
  CIb_out <- apply(CI_bonferroni, 2, fmt)
  rownames(CIb_out) <- vnames
  print(CIb_out)
  
  cat("============================================================\n\n")
  
  invisible(list(
    n = n, p = p, alpha = alpha, diff = diff,
    D = D, dbar = dbar, Sd = Sd, Sd_inv = Sd_inv,
    T2 = T2, F_stat = F_stat, df1 = df1, df2 = df2,
    F_crit = F_crit, T2_crit = T2_crit, decision = decision,
    CI_hotelling = CI_hotelling,
    CI_bonferroni = CI_bonferroni,
    tcrit_bonferroni = tcrit_bonf
  ))
}


# ================================
# END OF FILE
# ================================


# ============================================================
# UPDATED FUNCTIONS (2025-12-31)
# - Updated PCA function (prints Corr(Y1,Xk) exactly like exam formula)
# - New general multivariate regression function that supports:
#     * any number of responses (Y columns)
#     * any number of predictors (X columns)
#     * either raw data OR given (mu,S)
#     * optional intercept and optional z0; prints your requested note if z0 is missing
# - Backward-compatible wrapper stat438_q4() (old name)
# ============================================================

# ============================================================
# STAT 438 — PCA Exam Function (FULL SCRIPT, COPY/PASTE)
# Works for:
# ✅ Covariance matrix S (like CLEP question)
# ✅ Correlation matrix R (diag=1)
# ✅ Eigenvalues + eigenvectors + PC equations
# ✅ Proportion + cumulative % and #PCs for threshold
# ✅ Correlations Corr(Y_i, X_k) EXACTLY like your blue formula:
#     rho_{Y_i, X_k} = e_{k,i} * sqrt(lambda_i) / sqrt(sigma_kk)
#   (and for correlation-matrix PCA where sigma_kk = 1, becomes sqrt(lambda_i)*e_{k,i})
# ✅ Confidence intervals for eigenvalues:
#    - "chisq" (classic)
#    - "bonferroni_z" (your sheet’s z-form)
#
# EXTRA SAFETY:
# ✅ Warns if matrix looks "scaled" (e.g., S/(n-1)) — if eigenvalues look too small.
# ✅ Default corr_pairs computes Corr(Y1, Xk) for ALL k (like your screenshot)
# ============================================================

pca_exam_from_R <- function(R,
                            threshold = 0.80,
                            n = NULL,
                            alpha = 0.05,
                            ci_for = 1:2,
                            ci_method = c("chisq", "bonferroni_z"),
                            # correlations to compute, e.g. list(c(1,1), c(1,2), c(1,3)) for Y1 with all Xk
                            corr_pairs = NULL,
                            var_names = NULL,
                            digits = 7,
                            show_equations = TRUE,
                            warn_scaled = TRUE) {

  options(scipen = 999)
  ci_method <- match.arg(ci_method)

  # ----------------------------
  # Checks
  # ----------------------------
  R <- as.matrix(R)
  if (nrow(R) != ncol(R)) stop("R must be a square matrix (p x p).")
  if (!is.numeric(R)) stop("R must be numeric.")
  if (any(!is.finite(R))) stop("R contains non-finite values.")

  p <- nrow(R)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  if (length(var_names) != p) stop("var_names must have length p.")

  # Determine if correlation-like matrix
  is_corr_like <- all(abs(diag(R) - 1) < 1e-10)

  # Optional warning: user accidentally passed scaled matrix (common mistake: S/(n-1))
  if (warn_scaled && !is.null(n) && !is_corr_like) {
    # If you divide S by (n-1), the trace becomes much smaller but structure stays same.
    # Heuristic: if eigenvalues are about 1/(n-1) of what you expect, it's hard to detect
    # without the original. We only warn when the trace is "unusually small" vs max entry.
    tr <- sum(diag(R))
    maxabs <- max(abs(R))
    if (tr < 1.5 * maxabs) {
      warning(
        "Heuristic warning: this matrix's trace is not much larger than its largest entry.\n",
        "If your exam gives a large covariance S and your eigenvalues look too small,\n",
        "make sure you typed S EXACTLY as given (not divided by n-1)."
      )
    }
  }

  # ----------------------------
  # Eigen decomposition
  # ----------------------------
  eig <- eigen(R)
  lambda <- as.numeric(eig$values)
  E <- as.matrix(eig$vectors)

  # Ensure descending order
  ord <- order(lambda, decreasing = TRUE)
  lambda <- lambda[ord]
  E <- E[, ord, drop = FALSE]

  colnames(E) <- paste0("Y", 1:p)
  rownames(E) <- var_names

  # ----------------------------
  # Proportions
  # ----------------------------
  prop <- lambda / sum(lambda)
  cumprop <- cumsum(prop)
  k_needed <- which(cumprop >= threshold)[1]

  # ----------------------------
  # PC equations
  # ----------------------------
  pc_equation <- function(i) {
    coefs <- E[, i]
    terms <- paste0(
      ifelse(coefs >= 0, "+ ", "- "),
      format(abs(round(coefs, digits)), nsmall = digits, scientific = FALSE),
      " ",
      var_names
    )
    terms[1] <- sub("^\\+\\s+", "", terms[1])
    paste0("Y", i, " = ", paste(terms, collapse = " "))
  }

  # ----------------------------
  # Correlations Corr(Y_i, X_k)
  # EXAM formula:
  #   rho_{Y_i, X_k} = e_{k,i} * sqrt(lambda_i) / sqrt(sigma_kk)
  # For correlation PCA, sigma_kk = 1, so rho = sqrt(lambda_i) * e_{k,i}
  # ----------------------------
  if (is.null(corr_pairs)) {
    # Default: compute Corr(Y1, Xk) for all variables (matches your screenshot)
    corr_pairs <- lapply(1:p, function(k) c(1, k))
  }

  corr_results <- sapply(corr_pairs, function(pair) {
    if (length(pair) != 2) stop("Each element of corr_pairs must be c(i, k).")
    i <- pair[1] # PC index (Y_i)
    k <- pair[2] # variable index (X_k)
    if (i < 1 || i > p || k < 1 || k > p) stop("corr_pairs indices out of range.")
    sigma_kk <- R[k, k]
    (E[k, i] * sqrt(lambda[i])) / sqrt(sigma_kk)
  })

  names(corr_results) <- sapply(corr_pairs, function(pair) {
    paste0("rho_{Y", pair[1], ",", var_names[pair[2]], "}")
  })

  # Breakdown for Corr(Y1, Xk)
  corr_breakdown <- NULL
  if (any(sapply(corr_pairs, function(x) x[1] == 1))) {
    idxY1 <- which(sapply(corr_pairs, function(x) x[1] == 1))
    corr_breakdown <- data.frame(
      k = sapply(corr_pairs[idxY1], function(x) x[2]),
      variable = var_names[sapply(corr_pairs[idxY1], function(x) x[2])],
      e_k1 = E[sapply(corr_pairs[idxY1], function(x) x[2]), 1],
      lambda1 = lambda[1],
      sigma_kk = diag(R)[sapply(corr_pairs[idxY1], function(x) x[2])],
      rho = as.numeric(corr_results[idxY1]),
      stringsAsFactors = FALSE
    )
  }

  # ----------------------------
  # Confidence intervals for eigenvalues
  # ----------------------------
  ci_table <- NULL
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1 || n <= 1) stop("n must be a single number > 1.")
    ci_for <- as.integer(ci_for)
    ci_for <- ci_for[ci_for >= 1 & ci_for <= p]

    if (length(ci_for) > 0) {
      ci_table <- data.frame(
        eigen = paste0("lambda", ci_for),
        LL = NA_real_,
        lambda_hat = lambda[ci_for],
        UL = NA_real_,
        stringsAsFactors = FALSE
      )

      if (ci_method == "chisq") {
        df <- n - 1
        for (idx in seq_along(ci_for)) {
          i <- ci_for[idx]
          ci_table$LL[idx] <- (df * lambda[i]) / qchisq(1 - alpha/2, df)
          ci_table$UL[idx] <- (df * lambda[i]) / qchisq(alpha/2, df)
        }
      } else if (ci_method == "bonferroni_z") {
        zcrit <- qnorm(1 - alpha/(2 * p))
        kfac <- zcrit * sqrt(2 / n)
        for (idx in seq_along(ci_for)) {
          i <- ci_for[idx]
          ci_table$LL[idx] <- lambda[i] / (1 + kfac)
          ci_table$UL[idx] <- lambda[i] / (1 - kfac)
        }
      }
    }
  }

  # ----------------------------
  # Printing
  # ----------------------------
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)

  cat("\n==================== PCA (Exam Output) ====================\n")
  cat(sprintf("p = %d variables\n", p))
  cat(sprintf("Matrix type detected: %s\n",
              ifelse(is_corr_like, "Correlation (diag=1)", "Covariance")))
  if (!is.null(n)) cat(sprintf("n = %d, alpha = %.4f\n", n, alpha))
  cat("===========================================================\n\n")

  if (show_equations) {
    cat("Key facts:\n")
    cat("  Eigen decomposition: R e_i = lambda_i e_i\n")
    cat("  PC_i = e_i' X\n")
    cat("  Proportion = lambda_i / sum(lambda)\n")
    cat("  Correlation (EXAM FORMULA): rho_{Y_i,X_k} = e_{k,i}*sqrt(lambda_i)/sqrt(sigma_kk)\n")
    cat("    - if R is correlation: sigma_kk = 1 => rho = sqrt(lambda_i)*e_{k,i}\n")
    cat("  CI methods:\n")
    cat("    - chisq: (n-1)lambda / chi^2\n")
    cat("    - bonferroni_z: lambda /(1 ± z(alpha/(2p))*sqrt(2/n))\n")
    cat("-----------------------------------------------------------\n\n")
  }

  cat("1) Eigenvalues (lambda):\n")
  print(setNames(as.numeric(fmt(lambda)), paste0("lambda", 1:p)))

  cat("\n2) Eigenvectors (columns are PCs Y1..Yp):\n")
  print(apply(E, 2, fmt))

  cat("\n3) Proportion of variance explained:\n")
  prop_out <- rbind(
    Proportion = as.numeric(fmt(prop)),
    Cumulative = as.numeric(fmt(cumprop))
  )
  colnames(prop_out) <- paste0("Y", 1:p)
  print(prop_out)

  cat(sprintf("\n4) #PCs to reach %.0f%% cumulative variance: %d PCs (cum = %s)\n",
              100 * threshold, k_needed, fmt(cumprop[k_needed])))

  cat("\n5) Principal component equations:\n")
  for (i in 1:p) cat("  ", pc_equation(i), "\n", sep = "")

  cat("\n6) Correlations Corr(Y_i, X_k) (using exam formula):\n")
  print(setNames(as.numeric(fmt(corr_results)), names(corr_results)))

  if (!is.null(corr_breakdown)) {
    cat("\n6A) Breakdown for Corr(Y1, Xk):  rho_{Y1,Xk} = e_{k1}*sqrt(lambda1)/sqrt(sigma_kk)\n")
    cb <- corr_breakdown
    cb$e_k1 <- as.numeric(fmt(cb$e_k1))
    cb$lambda1 <- as.numeric(fmt(cb$lambda1))
    cb$sigma_kk <- as.numeric(fmt(cb$sigma_kk))
    cb$rho <- as.numeric(fmt(cb$rho))
    print(cb, row.names = FALSE)
  }

  if (!is.null(ci_table)) {
    cat("\n7) Confidence intervals for selected eigenvalues:\n")
    cat(sprintf("   Method: %s\n", ci_method))
    ci_print <- ci_table
    ci_print$LL <- as.numeric(fmt(ci_print$LL))
    ci_print$lambda_hat <- as.numeric(fmt(ci_print$lambda_hat))
    ci_print$UL <- as.numeric(fmt(ci_print$UL))
    print(ci_print, row.names = FALSE)
  } else {
    cat("\n7) CI not computed (provide n).\n")
  }

  cat("===========================================================\n\n")

  invisible(list(
    p = p,
    n = n,
    matrix_type = ifelse(is_corr_like, "correlation", "covariance"),
    lambda = lambda,
    E = E,
    prop = prop,
    cumprop = cumprop,
    k_threshold = k_needed,
    corr = corr_results,
    corr_breakdown_Y1 = corr_breakdown,
    ci = ci_table
  ))
}


# ============================================================
# STAT 438 — Multivariate Regression (EXAM SCRIPT) — ONE FUNCTION
# Works in BOTH exam styles:
#
# (A) RAW DATA GIVEN:
#     Provide Y (n x m) and X (n x r)
#
# (B) ONLY mu and S GIVEN:
#     Provide mu, S, n, and ALSO tell the function which variables are Y vs X
#     via y_idx and x_idx (indices inside mu/S).
#
# OUTPUT (exam style):
# 1) OLS B_hat (only if raw data)
# 2) Fitted values + residuals (only if raw data)
# 3) SSCP check (only if raw data)
# 4) MLE regression via mu & S  (always: from sample if raw data; or directly if mu/S given)
# 5) Partial correlation matrix Corr(Y | X) + one requested pair
# 6) Sigma_hat (MLE):
#       - from residuals if raw data
#       - from conditional covariance Sigma_YY|X if only mu/S
# 7-8) CI/PI (optional):
#     - only available when raw data is given AND intercept=TRUE AND z0 is provided
#     - if z0 is NULL => DOES NOT STOP, prints a note at the end exactly as requested
# ============================================================

stat438_reg_exam <- function(
  # --- Mode A: raw data ---
  Y = NULL,                      # (n x m)
  X = NULL,                      # (n x r)

  # --- Mode B: mu/S given ---
  mu = NULL,                     # vector of means for ALL variables (Y and X)
  S  = NULL,                     # covariance matrix for ALL variables (same order as mu)
  n  = NULL,                     # sample size (required in mu/S mode; printed in raw mode)
  y_idx = NULL,                  # indices of response variables inside mu/S
  x_idx = NULL,                  # indices of predictor variables inside mu/S

  # --- General options ---
  intercept = TRUE,
  digits = 4,
  alpha = 0.05,
  z0 = NULL,                     # if intercept=TRUE: c(1, x1_0, x2_0, ...)
  response_index = 1,            # which response column for CI/PI (raw-data mode only)
  do_intervals = TRUE,
  df2_mode = c("n-r-m-1", "n-r-m"),

  # Partial correlation reporting
  partial_pair = c(1, 2),        # which two responses (by column index in Y-order)
  show_partial_matrix = TRUE,
  show_equations = TRUE
) {
  options(scipen = 999)
  df2_mode <- match.arg(df2_mode)

  # ----------------------------
  # Helpers
  # ----------------------------
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)

  safe_solve <- function(A) {
    A <- as.matrix(A)
    if (qr(A)$rank < nrow(A)) stop("Matrix is singular (cannot invert). Check your data/inputs.")
    solve(A)
  }

  make_corr <- function(Sigma) {
    Sigma <- as.matrix(Sigma)
    D <- diag(1 / sqrt(diag(Sigma)))
    D %*% Sigma %*% D
  }

  cat_header <- function(title) {
    cat("\n============================================================\n")
    cat(title, "\n")
    cat("============================================================\n")
  }

  # ----------------------------
  # Decide which mode to use
  # ----------------------------
  using_muS <- (!is.null(mu) && !is.null(S))

  if (using_muS) {
    # ============ MODE B ============
    if (is.null(n)) stop("When mu and S are provided, you must also provide n (sample size).")
    if (is.null(y_idx) || is.null(x_idx)) {
      stop("When mu and S are provided, you must also provide y_idx and x_idx (indices inside mu/S).")
    }

    mu <- as.numeric(mu)
    S  <- as.matrix(S)
    if (length(mu) != nrow(S) || nrow(S) != ncol(S)) stop("mu length must match dimensions of S.")

    y_idx <- as.integer(y_idx)
    x_idx <- as.integer(x_idx)
    if (anyDuplicated(c(y_idx, x_idx)) > 0) stop("y_idx and x_idx must NOT overlap.")
    if (min(c(y_idx, x_idx)) < 1 || max(c(y_idx, x_idx)) > length(mu)) stop("y_idx/x_idx out of range.")

    m <- length(y_idx)
    r <- length(x_idx)
    if (m < 2) stop("Need at least 2 responses (length(y_idx) >= 2).")
    if (r < 1) stop("Need at least 1 predictor (length(x_idx) >= 1).")

    # Names (optional)
    if (is.null(names(mu))) names(mu) <- paste0("V", seq_along(mu))
    y_names <- names(mu)[y_idx]
    x_names <- names(mu)[x_idx]

    mu_Y <- matrix(mu[y_idx], ncol = 1)
    mu_X <- matrix(mu[x_idx], ncol = 1)

    S_YY <- S[y_idx, y_idx, drop = FALSE]
    S_YX <- S[y_idx, x_idx, drop = FALSE]
    S_XX <- S[x_idx, x_idx, drop = FALSE]

    # MLE regression (key exam formula)
    B_MLE <- S_YX %*% safe_solve(S_XX)                 # (m x r)
    beta0_MLE <- mu_Y - B_MLE %*% mu_X                 # (m x 1)

    rownames(B_MLE) <- y_names
    colnames(B_MLE) <- x_names
    rownames(beta0_MLE) <- y_names
    colnames(beta0_MLE) <- "Intercept"

    # Conditional covariance and partial correlations
    Sigma_YY_given_X <- S_YY - S_YX %*% safe_solve(S_XX) %*% t(S_YX)
    R_YY_given_X <- make_corr(Sigma_YY_given_X)
    rownames(R_YY_given_X) <- y_names
    colnames(R_YY_given_X) <- y_names

    # Best available "error covariance" in mu/S mode
    Sigma_hat <- Sigma_YY_given_X
    rownames(Sigma_hat) <- y_names
    colnames(Sigma_hat) <- y_names

    # Print
    cat_header("STAT 438 — Multivariate Regression (Using PROVIDED mu and S)")
    cat(sprintf("n = %d, m = %d responses, r = %d predictors, alpha = %.4f\n\n", n, m, r, alpha))

    if (show_equations) {
      cat("Equations (mu/S mode):\n")
      cat("  B_MLE      = S_YX S_XX^{-1}\n")
      cat("  beta0_MLE  = mu_Y - B_MLE mu_X\n")
      cat("  Sigma_YY|X = S_YY - S_YX S_XX^{-1} S_XY\n")
      cat("------------------------------------------------------------\n\n")
    }

    cat("Given mu =\n"); print(round(mu, digits))
    cat("\nGiven S =\n"); print(round(S, digits))

    cat("\n4) MLE regression using mean vector and covariance matrix\n\n")
    cat("B_MLE = S_YX S_XX^{-1} =\n")
    print(round(B_MLE, digits))

    cat("\nbeta0_MLE = mu_Y - B_MLE mu_X =\n")
    print(round(beta0_MLE, digits))

    cat("\nEstimated regression equations (MLE form):\n")
    for (i in 1:m) {
      rhs <- paste0(round(beta0_MLE[i, 1], digits))
      for (j in 1:r) {
        b <- B_MLE[i, j]
        rhs <- paste0(rhs,
                      ifelse(b >= 0, " + ", " - "),
                      round(abs(b), digits), " ", x_names[j])
      }
      cat(paste0(y_names[i], "_hat = ", rhs, "\n"))
    }

    cat("\n5) Partial correlation\n\n")
    cat("Sigma_YY|X =\n")
    print(round(Sigma_YY_given_X, digits))

    if (show_partial_matrix) {
      cat("\nCorr(Y | X) =\n")
      print(round(R_YY_given_X, digits))
    }

    a <- partial_pair[1]; b <- partial_pair[2]
    if (a < 1 || a > m || b < 1 || b > m) stop("partial_pair out of range for responses.")
    cat(sprintf("\nRequested r_{%s,%s · X} = %s\n", y_names[a], y_names[b], fmt(R_YY_given_X[a, b])))

    cat("\n6) Error covariance (best available in mu/S mode)\n\n")
    cat("Sigma_hat (reported as Sigma_YY|X) =\n")
    print(round(Sigma_hat, digits))

    if (do_intervals) {
      cat("\n7–8) Intervals skipped:\n")
      cat("Raw data were not provided, so Z'Z and residual-based interval formulas cannot be computed here.\n")
    }

    cat("\n============================== END ==============================\n\n")

    return(invisible(list(
      mode = "muS",
      n = n, m = m, r = r,
      y_names = y_names, x_names = x_names,
      mu = mu, S = S,
      B_MLE = B_MLE, beta0_MLE = beta0_MLE,
      Sigma_YY_given_X = Sigma_YY_given_X,
      R_YY_given_X = R_YY_given_X,
      Sigma_hat = Sigma_hat
    )))
  }

  # ============ MODE A ============
  if (is.null(Y) || is.null(X)) stop("Raw-data mode requires BOTH Y and X (matrices/data.frames).")

  Y <- as.matrix(Y)
  X <- as.matrix(X)
  if (!is.numeric(Y) || !is.numeric(X)) stop("Y and X must be numeric.")
  n <- nrow(Y)
  m <- ncol(Y)
  r <- ncol(X)
  if (m < 2) stop("Need at least 2 responses: ncol(Y) >= 2.")
  if (r < 1) stop("Need at least 1 predictor: ncol(X) >= 1.")
  if (nrow(X) != n) stop("X must have the same number of rows as Y.")

  if (is.null(colnames(Y))) colnames(Y) <- paste0("y", 1:m)
  if (is.null(colnames(X))) colnames(X) <- paste0("z", 1:r)

  if (response_index < 1 || response_index > m) stop("response_index out of range.")

  # Design matrix Z
  if (intercept) {
    Z <- cbind(1, X)
    colnames(Z)[1] <- "Intercept"
    intercept_msg <- "Intercept was NOT provided explicitly → added automatically (intercept = TRUE)."
  } else {
    Z <- X
    intercept_msg <- "Model fitted WITHOUT intercept (intercept = FALSE)."
  }
  pZ <- ncol(Z)

  cat_header("STAT 438 — Multivariate Regression (Using RAW DATA)")
  cat(sprintf("n = %d, m = %d responses, r = %d predictors, alpha = %.4f\n", n, m, r, alpha))
  cat("Design matrix:\n")
  cat("  ", intercept_msg, "\n", sep = "")
  cat("  Predictors: ", paste(colnames(X), collapse = ", "), "\n\n", sep = "")

  if (show_equations) {
    cat("Equations:\n")
    cat("  B_hat = (Z'Z)^(-1) Z'Y\n")
    cat("  Y_hat = Z B_hat\n")
    cat("  E_hat = Y - Y_hat\n")
    cat("  SSCP: Y'Y = Yhat'Yhat + E'E\n")
    cat("  (Joint) B_MLE = S_YX S_XX^{-1},  beta0_MLE = mu_Y - B_MLE mu_X\n")
    cat("  Sigma_hat(MLE) = (1/n) E'E\n")
    cat("------------------------------------------------------------\n\n")
  }

  # 1) OLS
  ZtZ_inv <- safe_solve(t(Z) %*% Z)
  B_hat <- ZtZ_inv %*% t(Z) %*% Y
  rownames(B_hat) <- colnames(Z)
  colnames(B_hat) <- colnames(Y)

  cat("1) Least squares estimates, B-hat\n\n")
  cat("B_hat = (Z'Z)^(-1) Z'Y =\n")
  print(round(B_hat, digits))

  # 2) fitted and residuals
  Y_hat <- Z %*% B_hat
  E_hat <- Y - Y_hat

  cat("\n2) Fitted values and residuals\n\n")
  cat("Y_hat = Z B_hat =\n")
  print(round(Y_hat, digits))

  cat("\nE_hat = Y - Y_hat =\n")
  print(round(E_hat, digits))

  # 3) SSCP check
  YY   <- t(Y) %*% Y
  YHYH <- t(Y_hat) %*% Y_hat
  EE   <- t(E_hat) %*% E_hat

  cat("\n3) Verify the SSCP decomposition\n\n")
  cat("Y'Y =\n")
  print(round(YY, digits))

  cat("\nY_hat'Y_hat + E_hat'E_hat =\n")
  print(round(YHYH + EE, digits))

  # 4) MLE regression via joint (Y,X)
  joint <- cbind(Y, X)
  mu_joint <- colMeans(joint)
  S_joint  <- stats::cov(joint)

  mu_Y <- matrix(mu_joint[1:m], ncol = 1)
  mu_X <- matrix(mu_joint[(m + 1):(m + r)], ncol = 1)

  S_YY <- S_joint[1:m, 1:m, drop = FALSE]
  S_YX <- S_joint[1:m, (m + 1):(m + r), drop = FALSE]
  S_XX <- S_joint[(m + 1):(m + r), (m + 1):(m + r), drop = FALSE]

  B_MLE <- S_YX %*% safe_solve(S_XX)
  beta0_MLE <- mu_Y - B_MLE %*% mu_X

  rownames(B_MLE) <- colnames(Y)
  colnames(B_MLE) <- colnames(X)
  rownames(beta0_MLE) <- colnames(Y)
  colnames(beta0_MLE) <- "Intercept"

  cat("\n4) MLE regression using mean vector and covariance matrix (from raw data)\n\n")
  cat("mu (joint) =\n"); print(round(mu_joint, digits))
  cat("\nS (joint) =\n"); print(round(S_joint, digits))

  cat("\nB_MLE = S_YX S_XX^{-1} =\n")
  print(round(B_MLE, digits))

  cat("\nbeta0_MLE = mu_Y - B_MLE mu_X =\n")
  print(round(beta0_MLE, digits))

  cat("\nEstimated regression equations (MLE form):\n")
  for (i in 1:m) {
    rhs <- paste0(round(beta0_MLE[i, 1], digits))
    for (j in 1:r) {
      b <- B_MLE[i, j]
      rhs <- paste0(rhs,
                    ifelse(b >= 0, " + ", " - "),
                    round(abs(b), digits), " ", colnames(X)[j])
    }
    cat(paste0(colnames(Y)[i], "_hat = ", rhs, "\n"))
  }

  # 5) partial correlation
  Sigma_YY_given_X <- S_YY - S_YX %*% safe_solve(S_XX) %*% t(S_YX)
  R_YY_given_X <- make_corr(Sigma_YY_given_X)
  rownames(R_YY_given_X) <- colnames(Y)
  colnames(R_YY_given_X) <- colnames(Y)

  a <- partial_pair[1]; b <- partial_pair[2]
  if (a < 1 || a > m || b < 1 || b > m) stop("partial_pair out of range for responses.")

  cat("\n5) Partial correlation\n\n")
  cat("Sigma_YY|X =\n")
  print(round(Sigma_YY_given_X, digits))

  if (show_partial_matrix) {
    cat("\nCorr(Y | X) =\n")
    print(round(R_YY_given_X, digits))
  }

  cat(sprintf("\nRequested r_{%s,%s · X} = %s\n", colnames(Y)[a], colnames(Y)[b], fmt(R_YY_given_X[a, b])))

  # 6) Sigma_hat (MLE) from residuals
  Sigma_hat <- (1 / n) * (t(E_hat) %*% E_hat)
  rownames(Sigma_hat) <- colnames(Y)
  colnames(Sigma_hat) <- colnames(Y)

  cat("\n6) Maximum likelihood estimate of Sigma\n\n")
  cat("Sigma_hat = (1/n) * E' E =\n")
  print(round(Sigma_hat, digits))

  # 7-8) CI / PI
  CI_mean <- NULL
  PI_pred <- NULL
  intervals_done <- FALSE

  if (do_intervals) {
    if (!intercept) {
      cat("\n7–8) Intervals skipped:\n")
      cat("Model has NO intercept (intercept=FALSE). Your course CI/PI formulas typically assume an intercept.\n")
    } else if (is.null(z0)) {
      cat("\n7–8) Intervals skipped:\n")
      cat("No z0 was provided for this question.\n")
    } else {
      z0 <- matrix(as.numeric(z0), ncol = 1)

      if (length(z0) != pZ) {
        cat("\n7–8) Intervals skipped:\n")
        cat("z0 length does NOT match the design matrix columns:\n")
        cat("  ", paste(colnames(Z), collapse = ", "), "\n", sep = "")
        cat(sprintf("So z0 must have length %d.\n", pZ))
      } else {
        df1 <- m
        df2 <- if (df2_mode == "n-r-m") (n - r - m) else (n - r - m - 1)
        if (df2 <= 0) {
          cat("\n7–8) Intervals skipped:\n")
          cat("df2 <= 0. Check n, r, m, and df2_mode.\n")
        } else {
          Fcrit <- qf(1 - alpha, df1, df2)
          denom <- if (df2_mode == "n-r-m") (n - r - m) else (n - r - m - 1)
          cfac  <- (m * (n - r - 1) / denom) * Fcrit

          # Course scaling frequently used in solutions:
          Sigma_tilde <- (n / (n - r - 1)) * Sigma_hat

          i <- response_index
          sigma_ii <- Sigma_tilde[i, i]

          b_i <- matrix(B_hat[, i], ncol = 1)
          mean_hat <- as.numeric(t(z0) %*% b_i)
          h <- as.numeric(t(z0) %*% ZtZ_inv %*% z0)

          hw_CI <- sqrt(cfac * h * sigma_ii)
          hw_PI <- sqrt(cfac * (1 + h) * sigma_ii)

          CI_mean <- c(mean_hat - hw_CI, mean_hat + hw_CI)
          PI_pred <- c(mean_hat - hw_PI, mean_hat + hw_PI)
          intervals_done <- TRUE

          cat(sprintf("\n7) %.0f%% CI for mean response E(%s | z0)\n\n", 100 * (1 - alpha), colnames(Y)[i]))
          cat("z0 =\n"); print(round(z0, digits))
          cat(sprintf("df2_mode = %s  =>  F ~ F_{%d, %g}\n", df2_mode, df1, df2))
          cat(sprintf("E(%s | z0) = z0' b_hat(i) = %s\n", colnames(Y)[i], fmt(mean_hat)))
          cat(sprintf("CI = (%s, %s)\n", fmt(CI_mean[1]), fmt(CI_mean[2])))

          cat(sprintf("\n8) %.0f%% PI for response %s at z0\n\n", 100 * (1 - alpha), colnames(Y)[i]))
          cat(sprintf("PI = (%s, %s)\n", fmt(PI_pred[1]), fmt(PI_pred[2])))
        }
      }
    }
  }

  # Your requested note when z0 is missing
  if (do_intervals && !intervals_done) {
    cat("\nNOTE (as requested):\n")
    cat("Intervals (CI/PI) were NOT computed.\n")
    if (is.null(z0) && intercept) {
      cat("Reason: No z0 (with intercept) was provided for this question.\n")
      cat("If you want to include the intercept, provide z0 like this:\n")
      cat("  Example: z0 = c(1, 1, 1)\n")
      cat("  (Here: 1 = intercept, then the predictor values follow.)\n")
    } else if (!intercept) {
      cat("Reason: Model has no intercept.\n")
    }
  }

  cat("\n============================== END ==============================\n\n")

  invisible(list(
    mode = "raw",
    n = n, m = m, r = r, intercept = intercept,
    Y = Y, X = X, Z = Z,
    ZtZ_inv = ZtZ_inv,
    B_hat = B_hat, Y_hat = Y_hat, E_hat = E_hat,
    SSCP_YtY = YY, SSCP_YhatYhat = YHYH, SSCP_EtE = EE,
    mu_joint = mu_joint, S_joint = S_joint,
    B_MLE = B_MLE, beta0_MLE = beta0_MLE,
    Sigma_YY_given_X = Sigma_YY_given_X,
    R_YY_given_X = R_YY_given_X,
    Sigma_hat = Sigma_hat,
    z0 = z0,
    CI_mean = CI_mean,
    PI_pred = PI_pred,
    df2_mode = df2_mode
  ))
}


# ------------------------------------------------------------
# Backward-compatible wrapper (old name):
# If you used stat438_q4(z1, y1, y2, z2=NULL, ...) before,
# you can still call it, and it will route to stat438_reg_exam().
# ------------------------------------------------------------
stat438_q4 <- function(z1, y1, y2, z2 = NULL,
                       intercept = TRUE,
                       digits = 4,
                       alpha = 0.05,
                       z0 = NULL,
                       response_index = 1,
                       do_intervals = TRUE,
                       df2_mode = c("n-r-m-1", "n-r-m")) {
  df2_mode <- match.arg(df2_mode)
  X <- if (is.null(z2)) cbind(z1) else cbind(z1, z2)
  colnames(X) <- if (is.null(z2)) "z1" else c("z1","z2")
  Y <- cbind(y1, y2)
  colnames(Y) <- c("y1","y2")
  stat438_reg_exam(
    Y = Y, X = X,
    intercept = intercept,
    digits = digits,
    alpha = alpha,
    z0 = z0,
    response_index = response_index,
    do_intervals = do_intervals,
    df2_mode = df2_mode,
    partial_pair = c(1,2)
  )
}



# =====================================================================
# UPDATED / ADDED FUNCTIONS (appended so they override earlier versions)
# =====================================================================

# ============================================================
# STAT 438 — Multivariate Regression (EXAM SCRIPT) — ONE FUNCTION
# Works in BOTH exam styles:
#
# (A) RAW DATA GIVEN:
#     Provide Y (n x m) and X (n x r)
#
# (B) ONLY mu and S GIVEN:
#     Provide mu, S, n, and ALSO tell the function which variables are Y vs X
#     via y_idx and x_idx (indices inside mu/S).
#
# OUTPUT (exam style):
# 1) OLS B_hat (only if raw data)
# 2) Fitted values + residuals (only if raw data)
# 3) SSCP check (only if raw data)
# 4) MLE regression via mu & S  (always: from sample if raw data; or directly if mu/S given)
# 5) Partial correlation matrix Corr(Y | X) + one requested pair
# 6) Sigma_hat (MLE):
#       - from residuals if raw data
#       - from conditional covariance Sigma_YY|X if only mu/S
# 7-8) CI/PI (optional):
#     - only available when raw data is given AND intercept=TRUE AND z0 is provided
#     - if z0 is NULL => DOES NOT STOP, prints a note at the end exactly as you requested
#
# ============================================================

stat438_reg_exam <- function(
    # --- Mode A: raw data ---
  Y = NULL,                      # (n x m)
  X = NULL,                      # (n x r)
  
  # --- Mode B: mu/S given ---
  mu = NULL,                     # vector of means for ALL variables (Y and X)
  S  = NULL,                     # covariance matrix for ALL variables (same order as mu)
  n  = NULL,                     # sample size (required for CI/PI critical values if you ever extend; also printed)
  y_idx = NULL,                  # indices of response variables inside mu/S
  x_idx = NULL,                  # indices of predictor variables inside mu/S
  
  # --- General options ---
  intercept = TRUE,
  digits = 4,
  alpha = 0.05,
  z0 = NULL,                     # if intercept=TRUE: c(1, x1_0, x2_0, ...)
  response_index = 1,            # which response column for CI/PI (raw-data mode only)
  do_intervals = TRUE,
  df2_mode = c("n-r-m-1", "n-r-m"),
  
  # Partial correlation reporting
  partial_pair = c(1, 2),        # which two responses (by column index in Y-order)
  show_partial_matrix = TRUE,
  show_equations = TRUE
) {
  options(scipen = 999)
  df2_mode <- match.arg(df2_mode)
  
  # ----------------------------
  # Helpers
  # ----------------------------
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  
  safe_solve <- function(A) {
    A <- as.matrix(A)
    if (qr(A)$rank < nrow(A)) stop("Matrix is singular (cannot invert). Check your data/contrasts.")
    solve(A)
  }
  
  make_corr <- function(Sigma) {
    Sigma <- as.matrix(Sigma)
    D <- diag(1 / sqrt(diag(Sigma)))
    D %*% Sigma %*% D
  }
  
  cat_header <- function(title) {
    cat("\n============================================================\n")
    cat(title, "\n")
    cat("============================================================\n")
  }
  
  # ----------------------------
  # Decide which mode to use
  # ----------------------------
  using_muS <- (!is.null(mu) && !is.null(S))
  
  if (using_muS) {
    # ============ MODE B ============
    # mu, S, n, y_idx, x_idx must be provided
    if (is.null(n)) stop("When mu and S are provided, you must also provide n (sample size).")
    if (is.null(y_idx) || is.null(x_idx)) {
      stop("When mu and S are provided, you must also provide y_idx and x_idx (indices inside mu/S).")
    }
    
    mu <- as.numeric(mu)
    S  <- as.matrix(S)
    
    if (length(mu) != nrow(S) || nrow(S) != ncol(S)) stop("mu length must match dimensions of S (square).")
    
    y_idx <- as.integer(y_idx)
    x_idx <- as.integer(x_idx)
    
    if (anyDuplicated(c(y_idx, x_idx)) > 0) stop("y_idx and x_idx must NOT overlap.")
    if (min(c(y_idx, x_idx)) < 1 || max(c(y_idx, x_idx)) > length(mu)) stop("y_idx/x_idx out of range for mu/S.")
    
    m <- length(y_idx)
    r <- length(x_idx)
    
    if (m < 2) stop("Need at least 2 responses (length(y_idx) >= 2).")
    if (r < 1) stop("Need at least 1 predictor (length(x_idx) >= 1).")
    
    # Names (optional)
    if (is.null(names(mu))) {
      names(mu) <- paste0("V", seq_along(mu))
    }
    y_names <- names(mu)[y_idx]
    x_names <- names(mu)[x_idx]
    
    mu_Y <- matrix(mu[y_idx], ncol = 1)
    mu_X <- matrix(mu[x_idx], ncol = 1)
    
    S_YY <- S[y_idx, y_idx, drop = FALSE]
    S_YX <- S[y_idx, x_idx, drop = FALSE]
    S_XX <- S[x_idx, x_idx, drop = FALSE]
    
    # MLE regression (this is the key exam formula when only mu/S are given)
    B_MLE <- S_YX %*% safe_solve(S_XX)                 # (m x r)
    beta0_MLE <- mu_Y - B_MLE %*% mu_X                 # (m x 1)
    
    rownames(B_MLE) <- y_names
    colnames(B_MLE) <- x_names
    rownames(beta0_MLE) <- y_names
    colnames(beta0_MLE) <- "Intercept"
    
    # Conditional covariance and partial correlations
    Sigma_YY_given_X <- S_YY - S_YX %*% safe_solve(S_XX) %*% t(S_YX)
    R_YY_given_X <- make_corr(Sigma_YY_given_X)
    rownames(R_YY_given_X) <- y_names
    colnames(R_YY_given_X) <- y_names
    
    # Sigma_hat (best available in mu/S mode)
    # With only mu and S, we do NOT have residuals, so we report Sigma_YY|X as the error covariance.
    Sigma_hat <- Sigma_YY_given_X
    rownames(Sigma_hat) <- y_names
    colnames(Sigma_hat) <- y_names
    
    # Print
    cat_header("STAT 438 — Multivariate Regression (Using PROVIDED mu and S)")
    cat(sprintf("n = %d, m = %d responses, r = %d predictors, alpha = %.4f\n\n", n, m, r, alpha))
    
    if (show_equations) {
      cat("Equations (mu/S mode):\n")
      cat("  B_MLE     = S_YX S_XX^{-1}\n")
      cat("  beta0_MLE = mu_Y - B_MLE mu_X\n")
      cat("  Sigma_YY|X = S_YY - S_YX S_XX^{-1} S_XY\n")
      cat("------------------------------------------------------------\n\n")
    }
    
    cat("Given mu =\n"); print(round(mu, digits))
    cat("\nGiven S =\n"); print(round(S, digits))
    
    cat("\n4) MLE regression using mean vector and covariance matrix\n\n")
    cat("B_MLE = S_YX S_XX^{-1} =\n")
    print(round(B_MLE, digits))
    
    cat("\nbeta0_MLE = mu_Y - B_MLE mu_X =\n")
    print(round(beta0_MLE, digits))
    
    cat("\nEstimated regression equations (MLE form):\n")
    for (i in 1:m) {
      rhs <- paste0(round(beta0_MLE[i, 1], digits))
      for (j in 1:r) {
        b <- B_MLE[i, j]
        rhs <- paste0(rhs,
                      ifelse(b >= 0, " + ", " - "),
                      round(abs(b), digits), " ", x_names[j])
      }
      cat(paste0(y_names[i], "_hat = ", rhs, "\n"))
    }
    
    cat("\n5) Partial correlation\n\n")
    cat("Sigma_YY|X =\n")
    print(round(Sigma_YY_given_X, digits))
    
    if (show_partial_matrix) {
      cat("\nCorr(Y | X) =\n")
      print(round(R_YY_given_X, digits))
    }
    
    a <- partial_pair[1]; b <- partial_pair[2]
    if (a < 1 || a > m || b < 1 || b > m) stop("partial_pair out of range for responses.")
    cat(sprintf("\nRequested r_{%s,%s · X} = %s\n",
                y_names[a], y_names[b], fmt(R_YY_given_X[a, b])))
    
    cat("\n6) Error covariance (best available in mu/S mode)\n\n")
    cat("Sigma_hat (reported as Sigma_YY|X) =\n")
    print(round(Sigma_hat, digits))
    
    # Intervals note (not computed in mu/S mode here)
    if (do_intervals) {
      cat("\n7–8) Intervals skipped:\n")
      cat("Raw data were not provided, so Z'Z and residual-based interval formulas cannot be computed here.\n")
    }
    
    cat("\n============================== END ==============================\n\n")
    
    invisible(list(
      mode = "muS",
      n = n, m = m, r = r,
      y_names = y_names, x_names = x_names,
      mu = mu, S = S,
      B_MLE = B_MLE, beta0_MLE = beta0_MLE,
      Sigma_YY_given_X = Sigma_YY_given_X,
      R_YY_given_X = R_YY_given_X,
      Sigma_hat = Sigma_hat
    ))
  }
  
  # ============ MODE A ============
  # Raw data required: Y and X
  if (is.null(Y) || is.null(X)) stop("Raw-data mode requires BOTH Y and X (matrices/data.frames).")
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  if (!is.numeric(Y) || !is.numeric(X)) stop("Y and X must be numeric.")
  n <- nrow(Y)
  m <- ncol(Y)
  r <- ncol(X)
  if (m < 2) stop("Need at least 2 responses: ncol(Y) >= 2.")
  if (r < 1) stop("Need at least 1 predictor: ncol(X) >= 1.")
  if (nrow(X) != n) stop("X must have the same number of rows as Y.")
  
  if (is.null(colnames(Y))) colnames(Y) <- paste0("y", 1:m)
  if (is.null(colnames(X))) colnames(X) <- paste0("z", 1:r)
  
  if (response_index < 1 || response_index > m) stop("response_index out of range.")
  
  # Design matrix Z
  if (intercept) {
    Z <- cbind(1, X)
    colnames(Z)[1] <- "Intercept"
    intercept_msg <- "Intercept was NOT provided explicitly → added automatically (intercept = TRUE)."
  } else {
    Z <- X
    intercept_msg <- "Model fitted WITHOUT intercept (intercept = FALSE)."
  }
  pZ <- ncol(Z)
  
  cat_header("STAT 438 — Multivariate Regression (Using RAW DATA)")
  cat(sprintf("n = %d, m = %d responses, r = %d predictors, alpha = %.4f\n", n, m, r, alpha))
  cat("Design matrix:\n")
  cat("  ", intercept_msg, "\n", sep = "")
  cat("  Predictors: ", paste(colnames(X), collapse = ", "), "\n\n", sep = "")
  
  if (show_equations) {
    cat("Equations:\n")
    cat("  B_hat = (Z'Z)^(-1) Z'Y\n")
    cat("  Y_hat = Z B_hat\n")
    cat("  E_hat = Y - Y_hat\n")
    cat("  SSCP: Y'Y = Yhat'Yhat + E'E\n")
    cat("  (Joint) B_MLE = S_YX S_XX^{-1},  beta0_MLE = mu_Y - B_MLE mu_X\n")
    cat("  Sigma_hat(MLE) = (1/n) E'E\n")
    cat("------------------------------------------------------------\n\n")
  }
  
  # 1) OLS
  ZtZ_inv <- safe_solve(t(Z) %*% Z)
  B_hat <- ZtZ_inv %*% t(Z) %*% Y
  rownames(B_hat) <- colnames(Z)
  colnames(B_hat) <- colnames(Y)
  
  cat("1) Least squares estimates, B-hat\n\n")
  cat("B_hat = (Z'Z)^(-1) Z'Y =\n")
  print(round(B_hat, digits))
  
  # 2) fitted and residuals
  Y_hat <- Z %*% B_hat
  E_hat <- Y - Y_hat
  
  cat("\n2) Fitted values and residuals\n\n")
  cat("Y_hat = Z B_hat =\n")
  print(round(Y_hat, digits))
  
  cat("\nE_hat = Y - Y_hat =\n")
  print(round(E_hat, digits))
  
  # 3) SSCP check
  YY   <- t(Y) %*% Y
  YHYH <- t(Y_hat) %*% Y_hat
  EE   <- t(E_hat) %*% E_hat
  
  cat("\n3) Verify the SSCP decomposition\n\n")
  cat("Y'Y =\n")
  print(round(YY, digits))
  
  cat("\nY_hat'Y_hat + E_hat'E_hat =\n")
  print(round(YHYH + EE, digits))
  
  # 4) MLE regression via joint (Y,X)
  joint <- cbind(Y, X)
  mu_joint <- colMeans(joint)
  S_joint  <- stats::cov(joint)
  
  mu_Y <- matrix(mu_joint[1:m], ncol = 1)
  mu_X <- matrix(mu_joint[(m + 1):(m + r)], ncol = 1)
  
  S_YY <- S_joint[1:m, 1:m, drop = FALSE]
  S_YX <- S_joint[1:m, (m + 1):(m + r), drop = FALSE]
  S_XX <- S_joint[(m + 1):(m + r), (m + 1):(m + r), drop = FALSE]
  
  B_MLE <- S_YX %*% safe_solve(S_XX)
  beta0_MLE <- mu_Y - B_MLE %*% mu_X
  
  rownames(B_MLE) <- colnames(Y)
  colnames(B_MLE) <- colnames(X)
  rownames(beta0_MLE) <- colnames(Y)
  colnames(beta0_MLE) <- "Intercept"
  
  cat("\n4) MLE regression using mean vector and covariance matrix (from raw data)\n\n")
  cat("mu (joint) =\n"); print(round(mu_joint, digits))
  cat("\nS (joint) =\n"); print(round(S_joint, digits))
  
  cat("\nB_MLE = S_YX S_XX^{-1} =\n")
  print(round(B_MLE, digits))
  
  cat("\nbeta0_MLE = mu_Y - B_MLE mu_X =\n")
  print(round(beta0_MLE, digits))
  
  cat("\nEstimated regression equations (MLE form):\n")
  for (i in 1:m) {
    rhs <- paste0(round(beta0_MLE[i, 1], digits))
    for (j in 1:r) {
      b <- B_MLE[i, j]
      rhs <- paste0(rhs,
                    ifelse(b >= 0, " + ", " - "),
                    round(abs(b), digits), " ", colnames(X)[j])
    }
    cat(paste0(colnames(Y)[i], "_hat = ", rhs, "\n"))
  }
  
  # 5) partial correlation
  Sigma_YY_given_X <- S_YY - S_YX %*% safe_solve(S_XX) %*% t(S_YX)
  R_YY_given_X <- make_corr(Sigma_YY_given_X)
  rownames(R_YY_given_X) <- colnames(Y)
  colnames(R_YY_given_X) <- colnames(Y)
  
  a <- partial_pair[1]; b <- partial_pair[2]
  if (a < 1 || a > m || b < 1 || b > m) stop("partial_pair out of range for responses.")
  
  cat("\n5) Partial correlation\n\n")
  cat("Sigma_YY|X =\n")
  print(round(Sigma_YY_given_X, digits))
  
  if (show_partial_matrix) {
    cat("\nCorr(Y | X) =\n")
    print(round(R_YY_given_X, digits))
  }
  
  cat(sprintf("\nRequested r_{%s,%s · X} = %s\n",
              colnames(Y)[a], colnames(Y)[b], fmt(R_YY_given_X[a, b])))
  
  # 6) Sigma_hat (MLE) from residuals
  Sigma_hat <- (1 / n) * (t(E_hat) %*% E_hat)
  rownames(Sigma_hat) <- colnames(Y)
  colnames(Sigma_hat) <- colnames(Y)
  
  cat("\n6) Maximum likelihood estimate of Sigma\n\n")
  cat("Sigma_hat = (1/n) * E' E =\n")
  print(round(Sigma_hat, digits))
  
  # 7-8) CI / PI
  CI_mean <- NULL
  PI_pred <- NULL
  intervals_done <- FALSE
  
  if (do_intervals) {
    if (!intercept) {
      cat("\n7–8) Intervals skipped:\n")
      cat("Model has NO intercept (intercept=FALSE). Your course CI/PI formulas typically assume an intercept.\n")
    } else if (is.null(z0)) {
      # YOUR REQUEST: don't stop — just print a note at the end
      cat("\n7–8) Intervals skipped:\n")
      cat("No z0 was provided for this question.\n")
    } else {
      z0 <- matrix(as.numeric(z0), ncol = 1)
      
      if (length(z0) != pZ) {
        cat("\n7–8) Intervals skipped:\n")
        cat("z0 length does NOT match the design matrix columns:\n")
        cat("  ", paste(colnames(Z), collapse = ", "), "\n", sep = "")
        cat(sprintf("So z0 must have length %d.\n", pZ))
      } else {
        df1 <- m
        df2 <- if (df2_mode == "n-r-m") (n - r - m) else (n - r - m - 1)
        if (df2 <= 0) {
          cat("\n7–8) Intervals skipped:\n")
          cat("df2 <= 0. Check n, r, m, and df2_mode.\n")
        } else {
          Fcrit <- qf(1 - alpha, df1, df2)
          denom <- if (df2_mode == "n-r-m") (n - r - m) else (n - r - m - 1)
          cfac  <- (m * (n - r - 1) / denom) * Fcrit
          
          # Course scaling frequently used in solutions:
          Sigma_tilde <- (n / (n - r - 1)) * Sigma_hat
          
          i <- response_index
          sigma_ii <- Sigma_tilde[i, i]
          
          b_i <- matrix(B_hat[, i], ncol = 1)
          mean_hat <- as.numeric(t(z0) %*% b_i)
          h <- as.numeric(t(z0) %*% ZtZ_inv %*% z0)
          
          hw_CI <- sqrt(cfac * h * sigma_ii)
          hw_PI <- sqrt(cfac * (1 + h) * sigma_ii)
          
          CI_mean <- c(mean_hat - hw_CI, mean_hat + hw_CI)
          PI_pred <- c(mean_hat - hw_PI, mean_hat + hw_PI)
          intervals_done <- TRUE
          
          cat(sprintf("\n7) %.0f%% CI for mean response E(%s | z0)\n\n", 100 * (1 - alpha), colnames(Y)[i]))
          cat("z0 =\n"); print(round(z0, digits))
          cat(sprintf("df2_mode = %s  =>  F ~ F_{%d, %g}\n", df2_mode, df1, df2))
          cat(sprintf("E(%s | z0) = z0' b_hat(i) = %s\n", colnames(Y)[i], fmt(mean_hat)))
          cat(sprintf("CI = (%s, %s)\n", fmt(CI_mean[1]), fmt(CI_mean[2])))
          
          cat(sprintf("\n8) %.0f%% PI for response %s at z0\n\n", 100 * (1 - alpha), colnames(Y)[i]))
          cat(sprintf("PI = (%s, %s)\n", fmt(PI_pred[1]), fmt(PI_pred[2])))
        }
      }
    }
  }
  
  # Final note exactly as requested
  if (do_intervals && !intervals_done) {
    cat("\nNOTE (as requested):\n")
    cat("Intervals (CI/PI) were NOT computed.\n")
    if (is.null(z0) && intercept) {
      cat("Reason: No z0 (with intercept) was provided for this question.\n")
      cat("If you want to include the intercept, provide z0 like this:\n")
      cat("  Example: z0 = c(1, 1, 1)\n")
      cat("  (Here: 1 = intercept, then the predictor values follow.)\n")
    } else if (!intercept) {
      cat("Reason: Model has no intercept.\n")
    }
  }
  
  cat("\n============================== END ==============================\n\n")
  
  invisible(list(
    mode = "raw",
    n = n, m = m, r = r, intercept = intercept,
    Y = Y, X = X, Z = Z,
    ZtZ_inv = ZtZ_inv,
    B_hat = B_hat, Y_hat = Y_hat, E_hat = E_hat,
    SSCP_YtY = YY, SSCP_YhatYhat = YHYH, SSCP_EtE = EE,
    mu_joint = mu_joint, S_joint = S_joint,
    B_MLE = B_MLE, beta0_MLE = beta0_MLE,
    Sigma_YY_given_X = Sigma_YY_given_X,
    R_YY_given_X = R_YY_given_X,
    Sigma_hat = Sigma_hat,
    z0 = z0,
    CI_mean = CI_mean,
    PI_pred = PI_pred,
    df2_mode = df2_mode
  ))
}

# ============================================================
# QUICK USAGE EXAMPLES
# ============================================================

# ---------- Example 1: RAW DATA (your simple vectors) ----------
# z1 <- c(2,3,3,6,7,9)
# y1 <- c(10,5,7,19,11,18)
# y2 <- c(15,9,3,25,7,13)
# X <- cbind(z1)
# Y <- cbind(y1, y2)
# stat438_reg_exam(Y=Y, X=X, z0=NULL)             # runs fine, prints note about z0
# stat438_reg_exam(Y=Y, X=X, z0=c(1,4), response_index=1)  # intervals at z1=4

# ---------- Example 2: RAW DATA with MANY predictors & responses ----------
# X <- cbind(z1, z2, z3)
# Y <- cbind(y1, y2, y3)
# stat438_reg_exam(Y=Y, X=X, z0=c(1, 2, -1, 0), response_index=2, partial_pair=c(1,3))

# ---------- Example 3: mu/S MODE ----------
# Suppose mu/S are for variables in order: (Y1,Y2,Z1,Z2)
# mu <- c(Y1=2.0, Y2=3.0, Z1=5.0, Z2=7.0)
# S  <- matrix(c(
#   4, 1, 2, 0,
#   1, 3, 1, 1,
#   2, 1, 5, 2,
#   0, 1, 2, 6
# ), 4,4, byrow=TRUE)
#
# stat438_reg_exam(mu=mu, S=S, n=40, y_idx=1:2, x_idx=3:4, partial_pair=c(1,2))


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


# ============================================================
# STAT 438 — PCA Exam Function (FULL SCRIPT, COPY/PASTE)
# Works for:
# ✅ Covariance matrix S (like CLEP question)
# ✅ Correlation matrix R (diag=1)
# ✅ Eigenvalues + eigenvectors + PC equations
# ✅ Proportion + cumulative % and #PCs for threshold
# ✅ Correlations Corr(Y_i, X_k) EXACTLY like your blue formula:
#     rho_{Y_i, X_k} = e_{k,i} * sqrt(lambda_i) / sqrt(sigma_kk)
#   (and for correlation-matrix PCA where sigma_kk = 1, becomes sqrt(lambda_i)*e_{k,i})
# ✅ Confidence intervals for eigenvalues:
#    - "chisq" (classic)
#    - "bonferroni_z" (your sheet’s z-form)
#
# EXTRA SAFETY:
# ✅ Warns if matrix looks "scaled" by (n-1)
# ✅ Always prints the Corr(Y1, Xk) requested in your screenshot style
# ============================================================

pca_exam_from_R <- function(R,
                            threshold = 0.80,
                            n = NULL,
                            alpha = 0.05,
                            ci_for = 1:2,
                            ci_method = c("chisq", "bonferroni_z"),
                            # correlations to compute, e.g. list(c(1,1), c(1,2), c(1,3)) for Y1 with all Xk
                            corr_pairs = NULL,
                            var_names = NULL,
                            digits = 7,
                            show_equations = TRUE,
                            warn_scaled = TRUE) {
  
  options(scipen = 999)
  ci_method <- match.arg(ci_method)
  
  # ----------------------------
  # Checks
  # ----------------------------
  R <- as.matrix(R)
  if (nrow(R) != ncol(R)) stop("R must be a square matrix (p x p).")
  if (!is.numeric(R)) stop("R must be numeric.")
  if (any(!is.finite(R))) stop("R contains non-finite values.")
  
  p <- nrow(R)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  if (length(var_names) != p) stop("var_names must have length p.")
  
  # Determine if correlation-like matrix
  is_corr_like <- all(abs(diag(R) - 1) < 1e-10)
  
  # Optional warning: user accidentally passed S/(n-1)
  if (warn_scaled && !is.null(n) && !is_corr_like) {
    tr <- sum(diag(R))
    # For covariance S, trace could be large (problem dependent).
    # If trace is suspiciously small relative to typical magnitude of entries,
    # warn that it might be divided by (n-1).
    maxabs <- max(abs(R))
    if (tr < 0.5 * maxabs && maxabs > 1) {
      warning("Your covariance matrix has a small trace relative to its entries.\n",
              "This sometimes happens if you accidentally used S/(n-1). ",
              "If your eigenvalues look too small, re-enter S exactly as given in the exam.")
    }
  }
  
  # ----------------------------
  # Eigen decomposition
  # ----------------------------
  eig <- eigen(R)
  lambda <- as.numeric(eig$values)
  E <- as.matrix(eig$vectors)
  
  # Ensure descending order
  ord <- order(lambda, decreasing = TRUE)
  lambda <- lambda[ord]
  E <- E[, ord, drop = FALSE]
  
  colnames(E) <- paste0("Y", 1:p)
  rownames(E) <- var_names
  
  # ----------------------------
  # Proportions
  # ----------------------------
  prop <- lambda / sum(lambda)
  cumprop <- cumsum(prop)
  k_needed <- which(cumprop >= threshold)[1]
  
  # ----------------------------
  # PC equations
  # ----------------------------
  pc_equation <- function(i) {
    coefs <- E[, i]
    terms <- paste0(
      ifelse(coefs >= 0, "+ ", "- "),
      format(abs(round(coefs, digits)), nsmall = digits, scientific = FALSE),
      " ",
      var_names
    )
    terms[1] <- sub("^\\+\\s+", "", terms[1])
    paste0("Y", i, " = ", paste(terms, collapse = " "))
  }
  
  # ----------------------------
  # Correlations Corr(Y_i, X_k)
  # EXACT exam formula (blue):
  #   rho_{Y_i, X_k} = e_{k,i} * sqrt(lambda_i) / sqrt(sigma_kk)
  # For correlation PCA, sigma_kk = 1, so rho = sqrt(lambda_i) * e_{k,i}
  # ----------------------------
  if (is.null(corr_pairs)) {
    # Default: compute Corr(Y1, Xk) for all variables (matches screenshot)
    corr_pairs <- lapply(1:p, function(k) c(1, k))
  }
  
  corr_results <- sapply(corr_pairs, function(pair) {
    if (length(pair) != 2) stop("Each element of corr_pairs must be c(i, k).")
    i <- pair[1] # PC index (Y_i)
    k <- pair[2] # variable index (X_k)
    if (i < 1 || i > p || k < 1 || k > p) stop("corr_pairs indices out of range.")
    sigma_kk <- R[k, k]
    (E[k, i] * sqrt(lambda[i])) / sqrt(sigma_kk)
  })
  
  names(corr_results) <- sapply(corr_pairs, function(pair) {
    paste0("rho_{Y", pair[1], ",", var_names[pair[2]], "}")
  })
  
  # Also build an "exam-style breakdown" for Corr(Y1, Xk)
  corr_breakdown <- NULL
  if (any(sapply(corr_pairs, function(x) x[1] == 1))) {
    idxY1 <- which(sapply(corr_pairs, function(x) x[1] == 1))
    corr_breakdown <- data.frame(
      k = sapply(corr_pairs[idxY1], function(x) x[2]),
      variable = var_names[sapply(corr_pairs[idxY1], function(x) x[2])],
      e_k1 = E[sapply(corr_pairs[idxY1], function(x) x[2]), 1],
      lambda1 = lambda[1],
      sigma_kk = diag(R)[sapply(corr_pairs[idxY1], function(x) x[2])],
      rho = as.numeric(corr_results[idxY1]),
      stringsAsFactors = FALSE
    )
  }
  
  # ----------------------------
  # Confidence intervals for eigenvalues
  # ----------------------------
  ci_table <- NULL
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1 || n <= 1) stop("n must be a single number > 1.")
    ci_for <- as.integer(ci_for)
    ci_for <- ci_for[ci_for >= 1 & ci_for <= p]
    
    if (length(ci_for) > 0) {
      ci_table <- data.frame(
        eigen = paste0("lambda", ci_for),
        LL = NA_real_,
        lambda_hat = lambda[ci_for],
        UL = NA_real_,
        stringsAsFactors = FALSE
      )
      
      if (ci_method == "chisq") {
        df <- n - 1
        for (idx in seq_along(ci_for)) {
          i <- ci_for[idx]
          ci_table$LL[idx] <- (df * lambda[i]) / qchisq(1 - alpha/2, df)
          ci_table$UL[idx] <- (df * lambda[i]) / qchisq(alpha/2, df)
        }
      } else if (ci_method == "bonferroni_z") {
        # your sheet:
        # lambda /(1 ± z(alpha/(2p))*sqrt(2/n))
        zcrit <- qnorm(1 - alpha/(2 * p))
        kfac <- zcrit * sqrt(2 / n)
        for (idx in seq_along(ci_for)) {
          i <- ci_for[idx]
          ci_table$LL[idx] <- lambda[i] / (1 + kfac)
          ci_table$UL[idx] <- lambda[i] / (1 - kfac)
        }
      }
    }
  }
  
  # ----------------------------
  # Printing
  # ----------------------------
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  
  cat("\n==================== PCA (Exam Output) ====================\n")
  cat(sprintf("p = %d variables\n", p))
  cat(sprintf("Matrix type detected: %s\n",
              ifelse(is_corr_like, "Correlation (diag=1)", "Covariance")))
  if (!is.null(n)) cat(sprintf("n = %d, alpha = %.4f\n", n, alpha))
  cat("===========================================================\n\n")
  
  if (show_equations) {
    cat("Key facts:\n")
    cat("  Eigen decomposition: R e_i = lambda_i e_i\n")
    cat("  PC_i = e_i' X\n")
    cat("  Proportion = lambda_i / sum(lambda)\n")
    cat("  Correlation (EXAM FORMULA): rho_{Y_i,X_k} = e_{k,i}*sqrt(lambda_i)/sqrt(sigma_kk)\n")
    cat("    - if R is correlation: sigma_kk = 1 => rho = sqrt(lambda_i)*e_{k,i}\n")
    cat("  CI methods:\n")
    cat("    - chisq: (n-1)lambda / chi^2\n")
    cat("    - bonferroni_z: lambda /(1 ± z(alpha/(2p))*sqrt(2/n))\n")
    cat("-----------------------------------------------------------\n\n")
  }
  
  cat("1) Eigenvalues (lambda):\n")
  print(setNames(as.numeric(fmt(lambda)), paste0("lambda", 1:p)))
  
  cat("\n2) Eigenvectors (columns are PCs Y1..Yp):\n")
  print(apply(E, 2, fmt))
  
  cat("\n3) Proportion of variance explained:\n")
  prop_out <- rbind(
    Proportion = as.numeric(fmt(prop)),
    Cumulative = as.numeric(fmt(cumprop))
  )
  colnames(prop_out) <- paste0("Y", 1:p)
  print(prop_out)
  
  cat(sprintf("\n4) #PCs to reach %.0f%% cumulative variance: %d PCs (cum = %s)\n",
              100 * threshold, k_needed, fmt(cumprop[k_needed])))
  
  cat("\n5) Principal component equations:\n")
  for (i in 1:p) cat("  ", pc_equation(i), "\n", sep = "")
  
  cat("\n6) Correlations Corr(Y_i, X_k) (using exam formula):\n")
  print(setNames(as.numeric(fmt(corr_results)), names(corr_results)))
  
  # Extra: explicitly show Corr(Y1,Xk) in the same structure as the screenshot
  if (!is.null(corr_breakdown)) {
    cat("\n6A) Breakdown for Corr(Y1, Xk):  rho_{Y1,Xk} = e_{k1}*sqrt(lambda1)/sqrt(sigma_kk)\n")
    corr_breakdown_print <- corr_breakdown
    corr_breakdown_print$e_k1 <- as.numeric(fmt(corr_breakdown_print$e_k1))
    corr_breakdown_print$lambda1 <- as.numeric(fmt(corr_breakdown_print$lambda1))
    corr_breakdown_print$sigma_kk <- as.numeric(fmt(corr_breakdown_print$sigma_kk))
    corr_breakdown_print$rho <- as.numeric(fmt(corr_breakdown_print$rho))
    print(corr_breakdown_print, row.names = FALSE)
  }
  
  if (!is.null(ci_table)) {
    cat("\n7) Confidence intervals for selected eigenvalues:\n")
    cat(sprintf("   Method: %s\n", ci_method))
    ci_print <- ci_table
    ci_print$LL <- as.numeric(fmt(ci_print$LL))
    ci_print$lambda_hat <- as.numeric(fmt(ci_print$lambda_hat))
    ci_print$UL <- as.numeric(fmt(ci_print$UL))
    print(ci_print, row.names = FALSE)
  } else {
    cat("\n7) CI not computed (provide n).\n")
  }
  
  cat("===========================================================\n\n")
  
  invisible(list(
    p = p,
    n = n,
    matrix_type = ifelse(is_corr_like, "correlation", "covariance"),
    lambda = lambda,
    E = E,
    prop = prop,
    cumprop = cumprop,
    k_threshold = k_needed,
    corr = corr_results,
    corr_breakdown_Y1 = corr_breakdown,
    ci = ci_table
  ))
}

# ============================================================
# HOW TO USE IT ON YOUR CLEP QUESTION (covariance S)
# ============================================================
# n <- 87
# S_exam <- matrix(c(
#   5691, 600, 217,
#   600,  126, 24,
#   217,   24, 23
# ), 3, 3, byrow = TRUE)
#
# pca_exam_from_R(
#   R = S_exam,
#   n = 87,
#   alpha = 0.01,
#   ci_method = "bonferroni_z",
#   ci_for = 1:3,
#   var_names = c("X1","X2","X3"),
#   digits = 7
# )
#
# This will ALSO print:
#   rho_{Y1,X1}, rho_{Y1,X2}, rho_{Y1,X3}
# using: e_{k1}*sqrt(lambda1)/sqrt(sigma_kk)


hotelling_T2_contrast <- function(xbar, S, n, C, alpha = 0.05, digits = 6, show_equations = TRUE) {
  
  # ----------------------------
  # Basic checks
  # ----------------------------
  xbar <- as.matrix(xbar)
  if (ncol(xbar) != 1) xbar <- matrix(xbar, ncol = 1)
  
  S <- as.matrix(S)
  C <- as.matrix(C)
  
  p <- nrow(S)
  q <- qr(C)$rank
  
  if (ncol(C) != p) stop("C must have p columns.")
  if (n <= q) stop("n must be larger than the number of contrasts q.")
  if (nrow(S) != ncol(S)) stop("S must be a square matrix.")
  
  # ----------------------------
  # Core calculations
  # ----------------------------
  delta_hat <- C %*% xbar
  S_delta   <- C %*% S %*% t(C)
  
  # Safe inverse (usually not needed in exams, but prevents crashes)
  if (qr(S_delta)$rank < nrow(S_delta)) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("S_delta is singular. Install MASS (install.packages('MASS')) or use full-rank contrasts.")
    }
    S_delta_i <- MASS::ginv(S_delta)
  } else {
    S_delta_i <- solve(S_delta)
  }
  
  T2 <- as.numeric(n * t(delta_hat) %*% S_delta_i %*% delta_hat)
  
  # ----------------------------
  # Sampling distribution + critical value
  # ----------------------------
  F_stat <- ((n - q) / (q * (n - 1))) * T2
  F_crit <- qf(1 - alpha, q, n - q)
  T2_crit <- (q * (n - 1) / (n - q)) * F_crit
  
  decision <- ifelse(T2 > T2_crit, "Reject H0", "Do NOT reject H0")
  
  # ----------------------------
  # (A) Hotelling simultaneous CIs for contrasts (components of C mu)
  # delta_i ± sqrt(T2_crit * S_delta[ii] / n)
  # ----------------------------
  se <- sqrt(diag(S_delta) / n)
  hw_hot <- sqrt(T2_crit) * se
  
  CI_hotelling <- cbind(
    Lower    = as.numeric(delta_hat) - hw_hot,
    Estimate = as.numeric(delta_hat),
    Upper    = as.numeric(delta_hat) + hw_hot
  )
  
  # ----------------------------
  # (B) Bonferroni simultaneous CIs for contrasts
  # delta_i ± t_{n-1, 1 - alpha/(2m)} * sqrt(S_delta[ii] / n)
  # where m = number of intervals = number of contrasts (rows of C)
  # ----------------------------
  m <- nrow(C)
  tcrit_bonf <- qt(1 - alpha / (2 * m), df = n - 1)
  hw_bonf <- tcrit_bonf * se
  
  CI_bonferroni <- cbind(
    Lower    = as.numeric(delta_hat) - hw_bonf,
    Estimate = as.numeric(delta_hat),
    Upper    = as.numeric(delta_hat) + hw_bonf
  )
  
  # ----------------------------
  # Print (EXAM FORMAT)
  # ----------------------------
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  
  cat("\n================ Hotelling's T^2 (Contrasts) ================\n")
  cat("H0 : C * mu = 0\n")
  cat(sprintf("n = %d, p = %d, q = %d, alpha = %.4f\n", n, p, q, alpha))
  cat("=============================================================\n\n")
  
  if (show_equations) {
    cat("Equations:\n")
    cat("  delta_hat = C * xbar\n")
    cat("  S_delta   = C * S * C'\n")
    cat("  T^2       = n * delta_hat' * (S_delta)^{-1} * delta_hat\n")
    cat("  F         = ((n-q)/(q(n-1))) * T^2  ~  F_{q, n-q}\n")
    cat("  T^2_crit  = (q(n-1)/(n-q)) * F_{q,n-q; 1-alpha}\n")
    cat("  Hotelling CI_i : delta_i ± sqrt(T^2_crit * S_delta[ii] / n)\n")
    cat("  Bonferroni CI_i: delta_i ± t_{n-1, 1-alpha/(2m)} * sqrt(S_delta[ii] / n),  m = #contrasts\n")
    cat("-------------------------------------------------------------\n\n")
  }
  
  cat("C =\n")
  print(C)
  
  cat("\n1) delta_hat = C * xbar:\n")
  print(matrix(fmt(delta_hat), ncol = 1))
  
  cat("\n2) S_delta = C * S * C':\n")
  print(matrix(fmt(S_delta), nrow = nrow(S_delta)))
  
  cat("\n3) T^2 statistic:\n")
  cat("T^2 =", fmt(T2), "\n")
  
  cat("\n4) Sampling distribution:\n")
  cat("((n - q)/(q(n - 1))) * T^2  ~  F_{q, n - q}\n")
  cat("q =", q, ", n =", n, "\n")
  
  cat("\n5) Critical value:\n")
  cat("T^2_crit =", fmt(T2_crit), "\n")
  
  cat("\n6) Decision:\n")
  cat(decision, "\n")
  
  cat("\n7A) (1-alpha) Hotelling simultaneous confidence intervals:\n")
  print(apply(CI_hotelling, 2, fmt))
  
  cat("\n7B) (1-alpha) Bonferroni simultaneous confidence intervals:\n")
  cat(sprintf("t_crit(Bonferroni) = %s  with df = %d and alpha/(2m) = %.6f\n",
              fmt(tcrit_bonf), n - 1, alpha / (2 * m)))
  print(apply(CI_bonferroni, 2, fmt))
  
  cat("=============================================================\n\n")
  
  # ----------------------------
  # Return everything
  # ----------------------------
  invisible(list(
    delta_hat = delta_hat,
    S_delta = S_delta,
    T2 = T2,
    F_stat = F_stat,
    T2_crit = T2_crit,
    decision = decision,
    CI_hotelling = CI_hotelling,
    CI_bonferroni = CI_bonferroni,
    tcrit_bonferroni = tcrit_bonf,
    q = q,
    m = m
  ))
}


# ============================================================
# STAT 438 — One-way MANOVA from summary stats (exam-style)
# Solves questions like your Q3:
# 1) Box's M test for equality of covariance matrices
# 2) Compute grand mean, B, W, T
# 3) Wilks' Lambda, F-approx + critical value + decision
#
# Inputs:
# - xbars: list of mean vectors (each length p)
# - Slist: list of covariance matrices (p x p) for each group
# - n:     vector of sample sizes for each group
# - alpha: significance level
#
# Notes:
# - Works with any p >= 2 and any g >= 2
# - Box's M uses a standard chi-square approximation
# - Wilks F-approx uses the common Rao's approximation (general)
#   and also prints the special-case sqrt(Lambda) form when p=2
# ============================================================

stat438_oneway_manova_summary <- function(xbars, Slist, n, alpha = 0.05, digits = 6,
                                          show_equations = TRUE) {
  
  # ---- checks ----
  if (!is.list(xbars) || !is.list(Slist)) stop("xbars and Slist must be lists.")
  g <- length(xbars)
  if (length(Slist) != g) stop("xbars and Slist must have the same length.")
  if (length(n) != g) stop("n must have length g (one sample size per group).")
  if (g < 2) stop("Need at least 2 groups.")
  
  # force matrices
  xbars <- lapply(xbars, function(x) matrix(as.numeric(x), ncol = 1))
  Slist <- lapply(Slist, function(S) as.matrix(S))
  
  p <- nrow(Slist[[1]])
  for (i in 1:g) {
    if (nrow(Slist[[i]]) != p || ncol(Slist[[i]]) != p) stop("All S matrices must be p x p.")
    if (nrow(xbars[[i]]) != p) stop("All mean vectors must have length p.")
    if (n[i] <= 1) stop("Each group size must be > 1.")
  }
  
  N <- sum(n)
  df_within <- sum(n - 1)
  
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  
  # ============================================================
  # Part A) Box's M test
  # ============================================================
  
  # pooled covariance
  Sp_num <- Reduce(`+`, lapply(1:g, function(i) (n[i] - 1) * Slist[[i]]))
  Sp <- Sp_num / df_within
  
  # determinants (use determinant() for numerical stability)
  logdet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)
  
  M <- df_within * logdet(Sp) - sum(sapply(1:g, function(i) (n[i]-1) * logdet(Slist[[i]])))
  
  # correction factor (standard Box M chi-square approx)
  # c = [ (2p^2 + 3p - 1) / (6(p+1)(g-1)) ] * [ sum 1/(n_i-1) - 1/df_within ]
  U_raw <- sum(1/(n - 1)) - 1/df_within
  c_fac <- ((2*p^2 + 3*p - 1) / (6*(p + 1)*(g - 1))) * U_raw
  
  D <- (1 - c_fac) * M
  df_box <- (g - 1) * p * (p + 1) / 2
  chi_crit <- qchisq(1 - alpha, df_box)
  box_decision <- ifelse(D > chi_crit, "Reject H0 (covariances NOT equal)", "Do NOT reject H0 (covariances can be treated equal)")
  
  # ============================================================
  # Part B) MANOVA SSCP matrices: grand mean, B, W, T
  # ============================================================
  
  # grand mean
  xbar <- Reduce(`+`, lapply(1:g, function(i) n[i] * xbars[[i]])) / N
  
  # W
  W <- Sp_num
  
  # B
  B <- Reduce(`+`, lapply(1:g, function(i) {
    d <- xbars[[i]] - xbar
    n[i] * (d %*% t(d))
  }))
  
  T <- B + W
  
  # ============================================================
  # Part C) Wilks' Lambda and F approximation
  # ============================================================
  
  Lambda <- det(W) / det(T)
  
  # General Rao's F approximation (works for general p,g)
  # Let s = sqrt( (p^2 (g-1)^2 - 4) / (p^2 + (g-1)^2 - 5) ) if denom>0 else 1
  # v = N - 1 - (p + g)/2
  # df1 = p(g-1)
  # df2 = s*v - (df1 - 2)/2
  # F = ((1 - Lambda^(1/s)) / (Lambda^(1/s))) * (df2/df1)
  df1 <- p * (g - 1)
  
  denom <- (p^2 + (g - 1)^2 - 5)
  if (denom > 0) {
    s <- sqrt((p^2 * (g - 1)^2 - 4) / denom)
  } else {
    s <- 1
  }
  
  v <- N - 1 - (p + g)/2
  # avoid negative/invalid df2 in tiny samples
  df2 <- max(1e-9, s * v - (df1 - 2)/2)
  
  Lambda_pow <- Lambda^(1/s)
  F_stat <- ((1 - Lambda_pow) / Lambda_pow) * (df2 / df1)
  F_crit <- qf(1 - alpha, df1, df2)
  wilks_decision <- ifelse(F_stat > F_crit, "Reject H0 (group means differ)", "Do NOT reject H0")
  
  # Special-case form often used in your course when p=2 (matches many model answers)
  special_F <- NULL
  special_df1 <- NULL
  special_df2 <- NULL
  special_Fcrit <- NULL
  if (p == 2) {
    # F = ((N-g-1)/(g-1)) * ((1 - sqrt(Lambda))/sqrt(Lambda))
    special_F <- ((N - g - 1) / (g - 1)) * ((1 - sqrt(Lambda)) / sqrt(Lambda))
    special_df1 <- 2 * (g - 1)
    special_df2 <- 2 * (N - g - 1)
    if (special_df2 > 0) special_Fcrit <- qf(1 - alpha, special_df1, special_df2)
  }
  
  # ============================================================
  # PRINT (EXAM STYLE)
  # ============================================================
  
  cat("\n================== One-way MANOVA (Summary) ==================\n")
  cat(sprintf("g = %d groups, p = %d responses, N = %d, alpha = %.3f\n", g, p, N, alpha))
  cat("==============================================================\n\n")
  
  if (show_equations) {
    cat("Formulas:\n")
    cat("  Sp = sum (n_i-1) S_i / sum(n_i-1)\n")
    cat("  M  = [sum(n_i-1)] ln|Sp| - sum (n_i-1) ln|S_i|\n")
    cat("  c  = [(2p^2+3p-1)/(6(p+1)(g-1))] * [sum 1/(n_i-1) - 1/sum(n_i-1)]\n")
    cat("  D  = (1-c) M  ~  Chi-square_{(g-1)p(p+1)/2}\n")
    cat("  xbar = (1/N) sum n_i xbar_i\n")
    cat("  W = sum (n_i-1) S_i\n")
    cat("  B = sum n_i (xbar_i - xbar)(xbar_i - xbar)'\n")
    cat("  T = B + W\n")
    cat("  Lambda = |W|/|T|\n")
    cat("--------------------------------------------------------------\n\n")
  }
  
  # 1) Box M
  cat("1) Box's M test for equal covariance matrices\n\n")
  cat("Pooled covariance Sp =\n")
  print(matrix(fmt(Sp), nrow = p))
  
  cat("\nU_raw = sum 1/(n_i-1) - 1/sum(n_i-1) = ", fmt(U_raw), "\n", sep="")
  cat("c (correction factor) = ", fmt(c_fac), "\n", sep="")
  cat("M = ", fmt(M), "\n", sep="")
  cat("D = (1-c)M = ", fmt(D), "\n", sep="")
  cat("Critical value (Chi-square) = ", fmt(chi_crit), " with df = ", df_box, "\n", sep="")
  cat("Decision: ", box_decision, "\n", sep="")
  
  # 2) SSCP
  cat("\n2) One-way MANOVA SSCP matrices\n\n")
  cat("Grand mean xbar =\n")
  print(matrix(fmt(xbar), ncol = 1))
  
  cat("\nBetween SSCP (B) =\n")
  print(matrix(fmt(B), nrow = p))
  
  cat("\nWithin SSCP (W) =\n")
  print(matrix(fmt(W), nrow = p))
  
  cat("\nTotal SSCP (T=B+W) =\n")
  print(matrix(fmt(T), nrow = p))
  
  # 3) Wilks
  cat("\n3) Wilks' Lambda test\n\n")
  cat("Lambda = |W|/|T| = ", fmt(Lambda), "\n", sep="")
  cat("Rao F-approx: F = ", fmt(F_stat), "  with df = (", fmt(df1), ", ", fmt(df2), ")\n", sep="")
  cat("Critical value = ", fmt(F_crit), "\n", sep="")
  cat("Decision: ", wilks_decision, "\n", sep="")
  
  if (!is.null(special_F) && !is.null(special_Fcrit)) {
    cat("\n(Extra) Course special-case when p=2:\n")
    cat("F = ((N-g-1)/(g-1)) * ((1 - sqrt(Lambda))/sqrt(Lambda))\n")
    cat("F_special = ", fmt(special_F), " with df = (", special_df1, ", ", special_df2, ")\n", sep="")
    cat("F_crit_special = ", fmt(special_Fcrit), "\n", sep="")
  }
  
  cat("==============================================================\n\n")
  
  invisible(list(
    g = g, p = p, n = n, N = N,
    Sp = Sp, M = M, U_raw = U_raw, c = c_fac, D = D, df_box = df_box,
    chi_crit = chi_crit, box_decision = box_decision,
    xbar = xbar, B = B, W = W, T = T,
    Lambda = Lambda,
    F_stat = F_stat, df1 = df1, df2 = df2, F_crit = F_crit, wilks_decision = wilks_decision,
    F_special = special_F, F_special_df1 = special_df1, F_special_df2 = special_df2, F_special_crit = special_Fcrit
  ))
}

# =======================
# Example (YOUR QUESTION)
# =======================
# x1 <- c(96, 92.3)
# x2 <- c(90.25, 88.25)
# x3 <- c(77.33, 73.33)
#
# S1 <- matrix(c(1, 0.5, 0.5, 0.33), 2, 2, byrow=TRUE)
# S2 <- matrix(c(30.92, -8.42, -8.42, 4.91), 2, 2, byrow=TRUE)
# S3 <- matrix(c(4.33, 1.83, 1.83, 2.33), 2, 2, byrow=TRUE)
#
# stat438_oneway_manova_summary(
#   xbars = list(x1, x2, x3),
#   Slist = list(S1, S2, S3),
#   n = c(3,3,3),
#   alpha = 0.05,
#   digits = 6
# )


# ============================================================
# STAT 438 — Chi-square Q-Q Plot + One-sample Hotelling T^2
# (FULL FUNCTION — copy/paste)
#
# Solves questions like:
# 1) Construct chi-square Q-Q plot table (D2 and Q)
# 2) Compute Hotelling T^2 for H0: mu = mu0
# 3) Sampling distribution of T^2
# 4) Critical value at alpha
# 5) Decision
# ============================================================

stat438_qq_hotelling_1sample <- function(X,
                                         mu0,
                                         alpha = 0.05,
                                         qq_prob = c("i-0.5_over_n", "i_over_n1"),
                                         digits = 4,
                                         make_plot = TRUE,
                                         show_equations = TRUE) {
  options(scipen = 999)
  qq_prob <- match.arg(qq_prob)
  
  # ----------------------------
  # Input checks
  # ----------------------------
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (any(!is.finite(X))) stop("X contains non-finite values.")
  n <- nrow(X)
  p <- ncol(X)
  if (n <= p) stop("Need n > p for Hotelling T^2 (and invertible S).")
  
  mu0 <- as.numeric(mu0)
  if (length(mu0) != p) stop(paste0("mu0 must have length p = ", p, "."))
  
  # ----------------------------
  # 1) Sample mean and covariance
  # ----------------------------
  xbar <- colMeans(X)
  S <- cov(X)
  
  # Invert S safely
  if (qr(S)$rank < p) stop("Sample covariance S is singular (not invertible).")
  S_inv <- solve(S)
  
  # ----------------------------
  # 2) Chi-square Q-Q plot values
  # D_i^2 = (x_i - xbar)' S^{-1} (x_i - xbar)
  # Q_i = chi^2_p( (i-0.5)/n )  (default exam choice)
  # ----------------------------
  D2 <- mahalanobis(X, center = xbar, cov = S)
  ord <- order(D2)
  D2_sorted <- D2[ord]
  
  probs <- if (qq_prob == "i-0.5_over_n") {
    ((1:n) - 0.5) / n
  } else {
    (1:n) / (n + 1)
  }
  Q <- qchisq(probs, df = p)
  
  qq_table <- data.frame(
    i = 1:n,
    D2 = D2_sorted,
    Q = Q
  )
  
  # ----------------------------
  # 3) Hotelling one-sample T^2 test
  # T^2 = n (xbar - mu0)' S^{-1} (xbar - mu0)
  # (n-p)/(p(n-1)) T^2 ~ F_{p, n-p}
  # T^2_crit = (p(n-1)/(n-p)) F_{p,n-p;1-alpha}
  # ----------------------------
  d <- matrix(xbar - mu0, ncol = 1)
  T2 <- as.numeric(n * t(d) %*% S_inv %*% d)
  
  F_stat <- ((n - p) / (p * (n - 1))) * T2
  F_crit <- qf(1 - alpha, df1 = p, df2 = n - p)
  T2_crit <- (p * (n - 1) / (n - p)) * F_crit
  
  decision <- if (T2 > T2_crit) "Reject H0" else "Do NOT reject H0"
  
  # ----------------------------
  # Printing (exam-ready)
  # ----------------------------
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  
  cat("\n================ STAT 438: Q-Q Plot + Hotelling T^2 ================\n")
  cat(sprintf("n = %d, p = %d, alpha = %.4f\n", n, p, alpha))
  cat("H0: mu = mu0\n")
  cat("====================================================================\n\n")
  
  if (show_equations) {
    cat("Equations:\n")
    cat("  xbar = (1/n) Σ x_i\n")
    cat("  S    = sample covariance\n")
    cat("  D_i^2 = (x_i - xbar)' S^{-1} (x_i - xbar)\n")
    cat("  Q_i   = χ²_p(prob_i)\n")
    cat("  T^2   = n (xbar - mu0)' S^{-1} (xbar - mu0)\n")
    cat("  ((n-p)/(p(n-1))) T^2  ~  F_{p, n-p}\n")
    cat("  T^2_crit = (p(n-1)/(n-p)) F_{p,n-p; 1-alpha}\n")
    cat("--------------------------------------------------------------------\n\n")
  }
  
  cat("A) Sample mean xbar:\n")
  print(matrix(fmt(xbar), ncol = 1))
  
  cat("\nB) Sample covariance S:\n")
  print(matrix(fmt(S), nrow = p))
  
  cat("\nC) Chi-square Q-Q table (sorted D2 and theoretical Q):\n")
  qq_print <- qq_table
  qq_print$D2 <- as.numeric(fmt(qq_print$D2))
  qq_print$Q  <- as.numeric(fmt(qq_print$Q))
  print(qq_print, row.names = FALSE)
  
  cat("\nD) Hotelling T^2 test for mu0:\n")
  cat("mu0 =\n"); print(matrix(fmt(mu0), ncol = 1))
  cat("\nT^2 =", fmt(T2), "\n")
  
  cat("\nE) Sampling distribution:\n")
  cat("((n-p)/(p(n-1))) * T^2  ~  F_{p, n-p}\n")
  cat(sprintf("=> F_{%d, %d}\n", p, n - p))
  
  cat("\nF) Critical value:\n")
  cat("T^2_crit =", fmt(T2_crit), "\n")
  
  cat("\nG) Decision:\n")
  cat(decision, "\n")
  
  cat("====================================================================\n\n")
  
  # ----------------------------
  # Optional plot
  # ----------------------------
  if (make_plot) {
    plot(Q, D2_sorted,
         pch = 19,
         xlab = expression(chi[p]^2~"quantiles (Q)"),
         ylab = expression("Ordered " * D^2),
         main = "Chi-square Q-Q Plot")
    abline(lm(D2_sorted ~ Q), lty = 2)
  }
  
  invisible(list(
    n = n,
    p = p,
    alpha = alpha,
    xbar = xbar,
    S = S,
    D2 = D2,
    D2_sorted = D2_sorted,
    Q = Q,
    qq_table = qq_table,
    mu0 = mu0,
    T2 = T2,
    F_stat = F_stat,
    T2_crit = T2_crit,
    decision = decision,
    order = ord
  ))
}

# =======================
# Example usage (your question)
# =======================
# X <- rbind(
#   c(9,12),
#   c(3,4),
#   c(6,6),
#   c(5,4),
#   c(7,9),
#   c(4,3),
#   c(6,7),
#   c(8,5)
# )
# stat438_qq_hotelling_1sample(X, mu0 = c(5,7), alpha = 0.05)


# ============================================================
# Paired Hotelling's T^2 (Before vs After, multivariate)
# Exam-style output: dbar, Sd, Sd^-1, T2, F, crit, decision,
# and BOTH simultaneous CIs:
#   (A) Hotelling T^2 simultaneous CIs
#   (B) Bonferroni simultaneous CIs
# ============================================================

# Small helper: "x %||% y" means if x is NULL use y
`%||%` <- function(x, y) if (is.null(x)) y else x

paired_hotelling_T2 <- function(before, after,
                                alpha = 0.05,
                                diff = c("before-after", "after-before"),
                                digits = 6,
                                show_equations = TRUE) {
  
  diff <- match.arg(diff)
  
  # ---- checks ----
  if (!is.matrix(before)) before <- as.matrix(before)
  if (!is.matrix(after))  after  <- as.matrix(after)
  
  if (nrow(before) != nrow(after) || ncol(before) != ncol(after)) {
    stop("before and after must have the same dimensions (n x p).")
  }
  
  n <- nrow(before)
  p <- ncol(before)
  
  if (n <= p) stop("Need n > p for the Hotelling T^2 to be valid (so F conversion works).")
  
  # ---- differences ----
  D <- if (diff == "before-after") before - after else after - before
  
  # ---- mean difference vector and covariance of differences ----
  dbar <- colMeans(D)
  Sd   <- stats::cov(D)          # sample covariance uses (n-1) by default
  Sd_inv <- solve(Sd)
  
  # ---- Hotelling T^2 test ----
  T2 <- as.numeric(n * t(dbar) %*% Sd_inv %*% dbar)
  
  # ---- Convert to F ----
  F_stat <- ((n - p) / (p * (n - 1))) * T2
  df1 <- p
  df2 <- n - p
  
  F_crit <- stats::qf(1 - alpha, df1, df2)
  T2_crit <- (p * (n - 1) / (n - p)) * F_crit
  
  decision <- if (T2 > T2_crit) "Reject H0" else "Do NOT reject H0"
  
  # ---- (A) Hotelling simultaneous CIs ----
  # delta_j: dbar_j ± sqrt(T2_crit * s_jj / n)
  se <- sqrt(diag(Sd) / n)
  hw_hotelling <- sqrt(T2_crit) * se
  
  CI_hotelling <- cbind(
    Lower = dbar - hw_hotelling,
    Estimate = dbar,
    Upper = dbar + hw_hotelling
  )
  
  # ---- (B) Bonferroni simultaneous CIs ----
  # delta_j: dbar_j ± t_{n-1, 1 - alpha/(2p)} * sqrt(s_jj / n)
  tcrit_bonf <- stats::qt(1 - alpha / (2 * p), df = n - 1)
  hw_bonf <- tcrit_bonf * se
  
  CI_bonferroni <- cbind(
    Lower = dbar - hw_bonf,
    Estimate = dbar,
    Upper = dbar + hw_bonf
  )
  
  # ---- printing helpers ----
  fmt <- function(x) format(round(x, digits), nsmall = digits, scientific = FALSE)
  vnames <- colnames(before) %||% paste0("V", 1:p)
  
  cat("\n============================================================\n")
  cat("Paired Hotelling's T^2 Test (Before vs After)\n")
  cat("============================================================\n")
  cat(sprintf("n = %d, p = %d, alpha = %.4f\n", n, p, alpha))
  cat(sprintf("Difference defined as: %s\n\n", diff))
  
  if (show_equations) {
    cat("Equations:\n")
    cat("  dbar = (1/n) * sum d_i\n")
    cat("  Sd   = (1/(n-1)) * sum (d_i - dbar)(d_i - dbar)'\n")
    cat("  T^2  = n * dbar' * Sd^{-1} * dbar\n")
    cat("  F    = ((n-p)/(p(n-1))) * T^2  ~  F_{p, n-p}\n")
    cat("  T^2_crit = (p(n-1)/(n-p)) * F_{p,n-p; 1-alpha}\n")
    cat("  Hotelling CI_j: dbar_j ± sqrt(T^2_crit * s_jj / n)\n")
    cat("  Bonferroni CI_j: dbar_j ± t_{n-1, 1-alpha/(2p)} * sqrt(s_jj / n)\n")
    cat("------------------------------------------------------------\n\n")
  }
  
  cat("1) Mean difference vector (dbar):\n")
  print(setNames(as.numeric(fmt(dbar)), vnames))
  
  cat("\n2) Covariance matrix of differences (Sd):\n")
  print(apply(Sd, 2, fmt) |>
          matrix(nrow = p, dimnames = list(vnames, vnames)))
  
  cat("\n3) Inverse covariance (Sd^{-1}):\n")
  print(apply(Sd_inv, 2, fmt) |>
          matrix(nrow = p, dimnames = list(vnames, vnames)))
  
  cat("\n4) Test statistics:\n")
  cat(sprintf("  T^2      = %s\n", fmt(T2)))
  cat(sprintf("  F        = %s  with df = (%d, %d)\n", fmt(F_stat), df1, df2))
  cat(sprintf("  F_crit   = %s\n", fmt(F_crit)))
  cat(sprintf("  T^2_crit = %s\n", fmt(T2_crit)))
  
  cat("\n5) Decision:\n")
  cat(sprintf("  %s  (since T^2 %s T^2_crit)\n",
              decision, ifelse(T2 > T2_crit, ">", "<=")))
  
  cat("\n6A) (1-alpha) Hotelling simultaneous CIs for delta:\n")
  CIh_out <- apply(CI_hotelling, 2, fmt)
  rownames(CIh_out) <- vnames
  print(CIh_out)
  
  cat("\n6B) (1-alpha) Bonferroni simultaneous CIs for delta:\n")
  cat(sprintf("    t_crit(Bonferroni) = %s  with df = %d and alpha/(2p) = %.6f\n",
              fmt(tcrit_bonf), n - 1, alpha / (2 * p)))
  CIb_out <- apply(CI_bonferroni, 2, fmt)
  rownames(CIb_out) <- vnames
  print(CIb_out)
  
  cat("============================================================\n\n")
  
  invisible(list(
    n = n, p = p, alpha = alpha, diff = diff,
    D = D, dbar = dbar, Sd = Sd, Sd_inv = Sd_inv,
    T2 = T2, F_stat = F_stat, df1 = df1, df2 = df2,
    F_crit = F_crit, T2_crit = T2_crit, decision = decision,
    CI_hotelling = CI_hotelling,
    CI_bonferroni = CI_bonferroni,
    tcrit_bonferroni = tcrit_bonf
  ))
}
