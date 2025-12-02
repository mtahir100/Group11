# viz_utils.R
library(ggplot2)
library(reshape2)

plot_matrix_heatmap <- function(mat, title = "") {
  # Downsample to prevent rendering lag if matrix is too large
  if (nrow(mat) > 150) mat <- mat[sample(nrow(mat), 150), ]
  if (ncol(mat) > 150) mat <- mat[, sample(ncol(mat), 150)]
  
  df <- reshape2::melt(as.matrix(mat))
  colnames(df) <- c("Row", "Col", "Value")
  
  # Automatic color scaling
  lim <- quantile(abs(df$Value), 0.98, na.rm = TRUE)
  
  ggplot(df, aes(Col, Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", 
                         midpoint = 0, limits = c(-lim, lim)) +
    scale_y_reverse() +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          panel.grid = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
}

