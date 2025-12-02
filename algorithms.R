# algorithms.R - Low-Rank Approximation Suite (SVD, NMF, CUR)

set.seed(1) # For reproductibility

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

pinv <- function(A, tol = 1e-10) {
  s = svd(A)
  d = s$d
  
  # Invert non-zero singular values
  d_inv <- ifelse(d > tol, 1 / d, 0)
  
  # Pseudoinverse formula
  A_pinv <- s$v %*% diag(d_inv) %*% t(s$u)
  return(A_pinv)
}

# CUR
cur_decompose <- function(A, k) {
  
  # Frobenius norm squared
  frob2 <- sum(A^2)
  
  # ---- Step 1: Column sampling probabilities ----
  col_probs <- colSums(A^2) / frob2 # probability of each column
  col_idx <- sample(1:ncol(A), k, replace = TRUE, prob = col_probs) # select k columns with higher proba of selecting higher proba columns
  C <- A[, col_idx, drop = FALSE]
  
  # Rescale selected columns
  for (i in 1:k) {
    C[, i] <- C[, i] / sqrt(k * col_probs[col_idx[i]])  #boost rare columns
  }
  
  # ---- Step 2: Row sampling probabilities ----
  row_probs <- rowSums(A^2) / frob2
  row_idx <- sample(1:nrow(A), k, replace = TRUE, prob = row_probs)
  R <- A[row_idx, , drop = FALSE]
  
  # Rescale selected rows
  for (i in 1:k) {
    R[i, ] <- R[i, ] / sqrt(k * row_probs[row_idx[i]])
  }
  
  # ---- Step 3: Compute U from intersection ----
  W <- A[row_idx, col_idx, drop = FALSE]
  
  # Compute pseudo-inverse of W
  pinvW <- pinv(W)
  
  # Middle matrix U
  U <- pinvW
  
  # ---- Step 4: CUR reconstruction ----
  A_cur <- C %*% U %*% R
  
  list(C = C, U = U, R = R, CUR = A_cur,
       col_idx = col_idx, row_idx = row_idx)
}

