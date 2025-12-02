# algorithms.R - Low-Rank Approximation Suite (SVD, NMF, CUR)

# Helper function: Calculate RMSE
calc_rmse <- function(X, X_recon) {
  sqrt(mean((X - X_recon)^2))
}

# Helper function: Calculate relative Frobenius error
calc_rel_error <- function(X, X_recon) {
  norm(X - X_recon, type = "F") / norm(X, type = "F")
}

# ==========================================
# 1. SVD Algorithms
# ==========================================

# Normal: Standard Truncated SVD (based on R's built-in LAPACK)
algo_svd_normal <- function(X, k) {
  # Enforce type check
  if(!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  s <- svd(X, nu = k, nv = k)
  
  # Reconstruct matrix
  D <- diag(s$d[1:k], nrow = k, ncol = k)
  recon <- s$u %*% D %*% t(s$v)
  
  list(recon = recon, components = s$u, name = "Deterministic SVD")
}

# Randomized: Randomized SVD (Halko et al.)
algo_svd_randomized <- function(X, k, n_oversamples = 10, n_iter = 2) {
  if(!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X); m <- ncol(X)
  p <- min(k + n_oversamples, m)
  
  # 1. Random projection to find range Q
  Omega <- matrix(rnorm(m * p), nrow = m)
  Y <- X %*% Omega
  
  # Power Iteration to enhance accuracy
  for (j in 1:n_iter) {
    Y <- X %*% crossprod(X, Y)
  }
  
  qr_res <- qr(Y)
  Q <- qr.Q(qr_res)
  
  # 2. Project to lower dimension B
  B <- crossprod(Q, X) # Q'X
  
  # 3. Compute SVD of B
  s <- svd(B, nu = k, nv = k)
  
  # 4. Recover U = Q * Ub
  U <- Q %*% s$u
  
  recon <- U %*% diag(s$d[1:k], nrow = k) %*% t(s$v)
  
  list(recon = recon, components = U, name = "Randomized SVD")
}

# ==========================================
# 2. NMF Algorithms
# ==========================================

# Core update logic (Multiplicative Updates)
run_nmf_core <- function(X, W, H, max_iter = 200, tol = 1e-4) {
  errors <- c()
  for (i in 1:max_iter) {
    # Update H
    WH <- W %*% H
    H <- H * (crossprod(W, X)) / (crossprod(W, WH) + 1e-9)
    
    # Update W
    WH <- W %*% H
    W <- W * (X %*% t(H)) / (WH %*% t(H) + 1e-9)
    
    # Error check
    if (i %% 10 == 0) {
      e <- norm(X - W %*% H, type = "F")
      errors <- c(errors, e)
      if (length(errors) > 1 && abs(errors[length(errors)] - errors[length(errors)-1]) < tol) break
    }
  }
  list(W = W, H = H, errors = errors)
}

# Normal: NMF with Random Uniform Initialization
algo_nmf_normal <- function(X, k) {
  X <- abs(X) # NMF requires non-negative input
  n <- nrow(X); m <- ncol(X)
  
  # Naive random initialization
  W_init <- matrix(runif(n*k), nrow=n)
  H_init <- matrix(runif(k*m), nrow=k)
  
  res <- run_nmf_core(X, W_init, H_init)
  list(recon = res$W %*% res$H, components = res$W, name = "Standard NMF (Random Init)")
}

# Randomized: NMF initialized via SVD (SVD-NMF)
# Use SVD results as a superior starting point for NMF (Warm Start)
algo_nmf_randomized <- function(X, k) {
  X <- abs(X)
  
  # Use Randomized SVD to get initialization points (abs to ensure non-negativity)
  svd_init <- algo_svd_randomized(X, k, n_iter = 1) 
  W_init <- abs(svd_init$components)
  
  # H could be estimated via least squares or taking V directly
  # Here we simplify by taking random to avoid singularity, focusing on W feature capture
  H_init <- matrix(runif(k * ncol(X)), nrow = k) 
  
  res <- run_nmf_core(X, W_init, H_init)
  list(recon = res$W %*% res$H, components = res$W, name = "Optimized NMF (rSVD Init)")
}

# ==========================================
# 3. CUR Matrix Decomposition
# ==========================================
# A approx C * U * R

compute_leverage_scores <- function(X, k) {
  # Compute leverage scores using SVD
  s <- svd(X, nu = k, nv = k)
  
  # Column leverage scores (based on V)
  lev_col <- rowSums(s$v[, 1:k]^2) / k
  # Row leverage scores (based on U)
  lev_row <- rowSums(s$u[, 1:k]^2) / k
  
  list(col = lev_col, row = lev_row)
}

# Normal: Deterministic CUR (Select Top Leverage Scores)
algo_cur_normal <- function(X, k) {
  # Selection count is usually a multiple of k, set to k+5 here to ensure rank
  c_num <- k + 5
  r_num <- k + 5
  c_num <- min(c_num, ncol(X)); r_num <- min(r_num, nrow(X))
  
  probs <- compute_leverage_scores(X, k)
  
  # Deterministically select the highest scores
  col_idx <- order(probs$col, decreasing = TRUE)[1:c_num]
  row_idx <- order(probs$row, decreasing = TRUE)[1:r_num]
  
  C <- X[, col_idx, drop = FALSE]
  R <- X[row_idx, , drop = FALSE]
  
  # Compute middle matrix U = pinv(C) * X * pinv(R)
  # Use manual pseudo-inverse to reduce dependencies
  pinv <- function(M) {
    s <- svd(M); tol <- 1e-5
    d_inv <- 1 / s$d; d_inv[s$d < tol] <- 0
    s$v %*% diag(d_inv, nrow=length(d_inv)) %*% t(s$u)
  }
  
  U_core <- pinv(C) %*% X %*% pinv(R)
  
  recon <- C %*% U_core %*% R
  list(recon = recon, components = C, name = "Deterministic CUR (Top Leverage)")
}

# Randomized: Randomized CUR (Probability sampling based on leverage scores)
algo_cur_randomized <- function(X, k) {
  c_num <- k + 5
  r_num <- k + 5
  c_num <- min(c_num, ncol(X)); r_num <- min(r_num, nrow(X))
  
  probs <- compute_leverage_scores(X, k)
  
  # Probabilistic sampling
  col_idx <- sample(1:ncol(X), c_num, prob = probs$col, replace = TRUE)
  row_idx <- sample(1:nrow(X), r_num, prob = probs$row, replace = TRUE)
  
  C <- X[, col_idx, drop = FALSE]
  R <- X[row_idx, , drop = FALSE]
  
  pinv <- function(M) {
    s <- svd(M); tol <- 1e-5
    d_inv <- 1 / s$d; d_inv[s$d < tol] <- 0
    s$v %*% diag(d_inv, nrow=length(d_inv)) %*% t(s$u)
  }
  
  U_core <- pinv(C) %*% X %*% pinv(R)
  recon <- C %*% U_core %*% R
  
  list(recon = recon, components = C, name = "Randomized CUR (Weighted Sampling)")
}

