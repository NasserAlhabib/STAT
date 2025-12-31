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
