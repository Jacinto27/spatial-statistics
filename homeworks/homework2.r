library(ggplot2)
library(dplyr)
library(viridis)

# setwd("C:/Users/Cun/Downloads/Spatial Statistics/Data")
load("class files/data/RainData.RData")

# Convert coordinates from meters to kilometers
dati <- dati %>% mutate(XUTM = XUTM / 1000, YUTM = YUTM / 1000)

# Define helper functions
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))

# Leave-One-Out Cross Validation for IDW
loocv_idw <- function(df_sf, p_vals = seq(1, 10, by = 0.2)) {
  z_vec <- st_drop_geometry(df_sf)$totanno
  scores <- numeric(length(p_vals))
  
  for (j in seq_along(p_vals)) {
    p <- p_vals[j]
    errs <- numeric(nrow(df_sf))
    
    for (i in seq_len(nrow(df_sf))) {
      train <- df_sf[-i, ]
      test  <- df_sf[i, ]
      pred  <- idw(totanno ~ 1, train, test, idp = p)$var1.pred
      errs[i] <- pred - z_vec[i]
    }
    scores[j] <- rmse(errs)
  }
  
  data.frame(p = p_vals, RMSE = scores) %>% arrange(RMSE)
}

# Create dense regular grid function
create_dense_grid <- function(df, n_cells = 100) {
  x_range <- range(df$XUTM)
  y_range <- range(df$YUTM)
  
  # Add 5% buffer around data
  x_buffer <- diff(x_range) * 0.05
  y_buffer <- diff(y_range) * 0.05
  
  x_seq <- seq(x_range[1] - x_buffer, x_range[2] + x_buffer, length.out = n_cells)
  y_seq <- seq(y_range[1] - y_buffer, y_range[2] + y_buffer, length.out = n_cells)
  
  dense_grid <- expand.grid(XUTM = x_seq, YUTM = y_seq)
  return(dense_grid)
}

# WORKING rainfall map function
create_rainfall_map <- function(df, best_p, year, n_cells = 100) {
  # Create dense grid
  dense_grid <- create_dense_grid(df, n_cells)
  
  # Convert to sf objects
  df_sf <- st_as_sf(df, coords = c("XUTM", "YUTM"), crs = NA)
  grid_sf <- st_as_sf(dense_grid, coords = c("XUTM", "YUTM"), crs = NA)
  
  # IDW interpolation
  idw_res <- idw(totanno ~ 1, df_sf, grid_sf, idp = best_p)
  
  interp_df <- data.frame(
    x = st_coordinates(grid_sf)[, "X"],
    y = st_coordinates(grid_sf)[, "Y"],
    pred = idw_res$var1.pred
  )
  
  cat("Year:", year, "| Rainfall range:", round(range(interp_df$pred), 1), "\n")
  
  # Create the map
  p <- ggplot() +
    geom_raster(data = interp_df, aes(x = x, y = y, fill = pred)) +
    geom_point(data = df, aes(x = XUTM, y = YUTM), color = "red", size = 1.5) +
    scale_fill_viridis_c(option = "C", name = "Rainfall\n(mm)") +
    labs(
      title = paste("Rainfall Map -", year, "(IDW, p =", best_p, ")"),
      x = "X (km)", 
      y = "Y (km)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

# Process all years
years <- c(1975, 1979, 1982)

for (yr in years) {
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("Processing year:", yr, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  df <- dati %>% filter(anno == yr)
  
  # 1. Exploratory Data Analysis
  cat("→ Running exploratory analysis...\n")
  print(summary(df$totanno))
  
  # Spatial distribution
  p_dist <- ggplot(df, aes(x = XUTM, y = YUTM, color = totanno)) +
    geom_point(size = 3) +
    scale_color_viridis_c(option = "C") +
    labs(title = paste("Rainfall Distribution -", yr),
         x = "X (km)", y = "Y (km)") +
    theme_minimal()
  print(p_dist)
  
  # 2. Find optimal p
  cat("→ Finding optimal p...\n")
  df_sf <- st_as_sf(df, coords = c("XUTM", "YUTM"), crs = NA)
  cv_results <- loocv_idw(df_sf)
  print(cv_results)
  
  best_p <- cv_results$p[which.min(cv_results$RMSE)]
  cat("✅ Optimal p =", best_p, "(RMSE =", round(min(cv_results$RMSE), 2), ")\n")
  
  # 3. Create final map
  cat("→ Creating rainfall map...\n")
  rain_map <- create_rainfall_map(df, best_p, yr)
  print(rain_map)
}

cat("\nANALYSIS COMPLETE! All maps generated successfully.\n"