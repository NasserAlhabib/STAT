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
