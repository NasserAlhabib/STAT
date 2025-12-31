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
