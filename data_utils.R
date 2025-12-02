# data_utils.R - Enforced Numeric Matrix Version
load_data <- function(dataset_name) {
  set.seed(42)
  
  if (dataset_name == "Movie Ratings (Simulated)") {
    n_users <- 100
    n_movies <- 50
    # Use 0.0 to ensure float matrix
    X <- matrix(0.0, nrow = n_users, ncol = n_movies)
    
    U_true <- matrix(abs(rnorm(n_users * 5, mean = 3, sd = 1)), nrow = n_users)
    V_true <- matrix(abs(rnorm(5 * n_movies, mean = 3, sd = 1)), nrow = 5)
    X_complete <- U_true %*% V_true
    
    # Add sparsity
    for (i in 1:n_users) {
      for (j in 1:n_movies) {
        if (runif(1) < 0.3) {
          # Force conversion to float
          X[i, j] <- as.numeric(pmax(1.0, pmin(5.0, X_complete[i, j] + rnorm(1, 0, 0.5))))
        }
      }
    }
    
    # ==== Double check: Force conversion to numeric matrix ====
    X <- as.matrix(X)
    storage.mode(X) <- "double"  # Explicitly double precision float
    X[is.na(X)] <- 0.0
    
    return(list(
      X = X,
      title = "Simulated Movie Ratings",
      description = paste("A", n_users, "x", n_movies, "numeric sparse matrix.")
    ))
    
  } else if (dataset_name == "Face Images (Simulated)") {
    n_pixels <- 64 * 64
    n_faces <- 50
    k_true <- 15
    
    A <- matrix(abs(rnorm(n_pixels * k_true)), nrow = n_pixels)
    B <- matrix(abs(rnorm(k_true * n_faces)), nrow = k_true)
    X <- A %*% B + matrix(rnorm(n_pixels * n_faces, 0, 0.1), nrow = n_pixels)
    
    # ==== Double check: Force conversion to numeric matrix ====
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    X[is.na(X)] <- 0.0
    
    return(list(
      X = X,
      title = "Simulated Face Images",
      description = paste("64x64 pixels, 50 samples. Low-rank + noise.")
    ))
  }
}

