library(ggplot2)
library(dplyr)

# 1. Create the data frame with your exact values
plot_data <- data.frame(
  Exposure = rep(c("PFOS", "PFOA", "PFNA", "PFHxS"), 2),
  Method = rep(c("CSDR (d=1)", "SDR (d=1)"), each = 4),
  Value = c(
    diag(Pi(beta)), # Causal SDR values
    diag(Pi(beta_assoc))  # Standard SDR values
  )
)

# 2. Reorder factor levels so "Causal SDR" appears on the top row, 
# and the exposures maintain their original order left-to-right.
plot_data$Method <- factor(plot_data$Method, levels = c("SDR (d=1)", "CSDR (d=1)"))
plot_data$Exposure <- factor(plot_data$Exposure, levels = c("PFOS", "PFOA", "PFNA", "PFHxS"))

# 3. Create the heatmap
heatmap <- ggplot(plot_data, aes(x = Exposure, y = Method, fill = Value)) +
  # Draw the squares with a small white border
  geom_tile(color = "white", linewidth = 1) +
  
  # Add the numeric values inside the squares (rounded to 3 decimal places)
  # Dynamically switch text color to white if the background is too dark
  geom_text(aes(label = sprintf("%.3f", Value), 
                color = Value > 0.6), 
            size = 6, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("black", "white")) +
  
  # Force the color scale to be exactly between 0 and 1
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1)) +
  
  # Clean, slide-ready theme
  theme_minimal(base_size = 18) +
  labs(
    #title = "Subspace importance score",
    title = NULL,
    x = NULL,
    y = NULL,
    fill = "Importance\n(0 to 1)"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black")
  )

ggsave(filename = "results/jasa-initial-submission/ATL-AA/diag_plot_notitle.pdf",
       plot = heatmap, width = 8, height = 4)


library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Define the methods and dimensions
methods_list <- c("PCA", "pCCA", "MAVE", "Oracle-MAVE", "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE")
dims_list <- paste0("Dim ", 1:10)

# Hardcode the exact mean values from the screenshot row by row
means_matrix <- matrix(c(
  0.002, 0.003, 0.013, 0.018, 0.969, 0.979, 0.007, 0.004, 0.003, 0.002, # PCA
  0.546, 0.716, 0.134, 0.134, 0.002, 0.002, 0.132, 0.122, 0.132, 0.081, # pCCA
  0.530, 0.663, 0.370, 0.022, 0.002, 0.001, 0.117, 0.110, 0.113, 0.072, # MAVE
  1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, # Oracle-MAVE
  0.880, 0.947, 0.050, 0.006, 0.000, 0.000, 0.033, 0.030, 0.031, 0.022, # RA-MAVE
  0.860, 0.928, 0.042, 0.005, 0.000, 0.001, 0.048, 0.045, 0.042, 0.029, # DR-MAVE
  0.770, 0.874, 0.044, 0.007, 0.001, 0.001, 0.088, 0.082, 0.079, 0.054, # PO-MAVE
  0.500, 0.580, 0.155, 0.163, 0.002, 0.002, 0.164, 0.171, 0.157, 0.105  # RP-MAVE
), byrow = TRUE, nrow = 8, dimnames = list(methods_list, dims_list))

# Convert the matrix into a long-format data frame for ggplot
plot_data <- as.data.frame(as.table(means_matrix))
colnames(plot_data) <- c("Method", "Dimension", "Value")

# 2. Reorder factor levels 
# rev() ensures "PCA" appears on the top row instead of the bottom
plot_data$Method <- factor(plot_data$Method, levels = rev(methods_list))
plot_data$Dimension <- factor(plot_data$Dimension, levels = dims_list)

# 3. Create the heatmap
heatmap <- ggplot(plot_data, aes(x = Dimension, y = Method, fill = Value)) +
  # Draw the squares with a small white border
  geom_tile(color = "white", linewidth = 1) +
  
  # Add the numeric values inside the squares (rounded to 3 decimal places)
  # Dynamically switch text color to white if the background is too dark
  geom_text(aes(label = sprintf("%.3f", Value), 
                color = Value > 0.6), 
            size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("black", "white")) +
  
  # Force the color scale to be exactly between 0 and 1
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1)) +
  
  # Clean, slide-ready theme
  theme_minimal(base_size = 18) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    fill = "Importance\n(0 to 1)"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    # Angled x-axis text to prevent overlap with 10 columns
    axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", color = "black")
  )

# Display the plot
print(heatmap)

# 4. Save the plot
ggsave(filename = "results/jasa-initial-submission/ATL-AA/diag_plot_n=500.pdf",
        plot = heatmap, width = 12, height = 5)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Define the methods and dimensions
methods_list <- c("PCA", "pCCA", "MAVE", "Oracle-MAVE", "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE")
dims_list <- paste0("X", 1:10)

# Hardcode the exact mean values from the n=5000 screenshot row by row
means_matrix <- matrix(c(
  0.002, 0.003, 0.011, 0.015, 0.972, 0.982, 0.007, 0.004, 0.003, 0.002, # PCA
  0.819, 0.884, 0.053, 0.053, 0.001, 0.001, 0.053, 0.053, 0.053, 0.031, # pCCA
  0.363, 0.782, 0.764, 0.010, 0.002, 0.000, 0.023, 0.021, 0.022, 0.013, # MAVE
  1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, # Oracle-MAVE
  0.999, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, # RA-MAVE
  0.996, 0.998, 0.002, 0.000, 0.000, 0.000, 0.001, 0.001, 0.001, 0.001, # DR-MAVE
  0.945, 0.978, 0.001, 0.001, 0.000, 0.000, 0.018, 0.021, 0.022, 0.013, # PO-MAVE
  0.827, 0.878, 0.050, 0.053, 0.001, 0.001, 0.052, 0.053, 0.052, 0.033  # RP-MAVE
), byrow = TRUE, nrow = 8, dimnames = list(methods_list, dims_list))

# Convert the matrix into a long-format data frame for ggplot
plot_data <- as.data.frame(as.table(means_matrix))
colnames(plot_data) <- c("Method", "Dimension", "Value")

# 2. Reorder factor levels 
# rev() ensures "PCA" appears on the top row instead of the bottom
plot_data$Method <- factor(plot_data$Method, levels = rev(methods_list))
plot_data$Dimension <- factor(plot_data$Dimension, levels = dims_list)

# 3. Create the heatmap
heatmap <- ggplot(plot_data, aes(x = Dimension, y = Method, fill = Value)) +
  # Draw the squares with a small white border
  geom_tile(color = "white", linewidth = 1) +
  
  # Add the numeric values inside the squares (rounded to 3 decimal places)
  # Dynamically switch text color to white if the background is too dark
  geom_text(aes(label = sprintf("%.3f", Value), 
                color = Value > 0.6), 
            size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("black", "white")) +
  
  # Force the color scale to be exactly between 0 and 1
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1)) +
  
  # Clean, slide-ready theme
  theme_minimal(base_size = 18) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    fill = "Importance\n(0 to 1)"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    # Angled x-axis text to prevent overlap with 10 columns
    axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", color = "black")
  )

# Display the plot
print(heatmap)

# 4. Save the plot
ggsave(filename = "results/jasa-initial-submission/ATL-AA/diag_plot_n=5000.pdf",
        plot = heatmap, width = 12, height = 5)