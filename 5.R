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
