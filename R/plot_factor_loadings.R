plot_factor_loadings <- function(mave, d){
  # Extract the 'd' MAVE object
  # QR retract and varimax rotate it
  # Plot the varimax rotated plot. 
  
  beta <- mave$dir[[d]]
  # Orthonormalize (preserve rownames) then varimax
  Q <- qr.Q(qr(beta))
  rownames(Q) <- rownames(beta)
  
  if(d > 1){
    beta_final <- Q |>
      as.matrix() |>
      varimax()
    beta.f <- unclass(beta_final$loadings) |>
      data.frame()
  } else{
    beta.f <- Q
  }
  # Clean row names:
  # 1. remove surrounding backticks if present
  # 2. remove "_m" suffix
  clean_names <- rownames(beta)
  clean_names <- gsub("^`|`$", "", clean_names)   # strip leading/trailing backticks
  clean_names <- sub("_m$", "", clean_names)      # drop trailing _m
  rownames(beta.f) <- clean_names
  colnames(beta.f) <- paste0("beta", 1:ncol(beta.f))
  # Convert to long format
  df_long <- data.frame(
    Cytokine = rep(rownames(beta.f), times = ncol(beta.f)),
    Factor   = rep(colnames(beta.f),  each = nrow(beta.f)),
    Loading  = as.vector(as.matrix(beta.f))
  )
  
  # Label text = actual loadings
  df_long$label <- sprintf("%.2f", df_long$Loading)
  
  # Dynamic text color: white on dark cells, black otherwise
  df_long$text_color <- ifelse(abs(df_long$Loading) > 0.4, "white", "black")
  
  # Make Factor an ordered factor, reverse levels so X1 is on top
  df_long$Factor <- factor(df_long$Factor, levels = rev(colnames(beta.f)))
  
  # Now plot
  ggplot(df_long, aes(x = Cytokine, y = Factor, fill = Loading)) +
    geom_tile(color = "grey80") +
    geom_text(aes(label = label, color = text_color), size = 5) +
    scale_color_identity() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    coord_fixed() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      legend.position = "none"
    ) + 
    ggtitle(paste0("d=",d)) + 
    theme(plot.title = element_text(size = 10))
  
}