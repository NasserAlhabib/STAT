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
