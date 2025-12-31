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
