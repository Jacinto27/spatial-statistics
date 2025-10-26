#! a) Build an R function to choose the parameter p in the inverse distance weighting estimator. Hint: use the idw implemented in gstat.
library(gstat)
library(sf)
library(geoR)
library(ggplot2)
data(wolfcamp)

df <- wolfcamp
#structure(df)
# remove old yucky geodata structure and turn into a dataframe NOTE: MBA package has wolfcamp as dataframe, step not needed if using
coords <- df$coords 
vals <- df$data
df2 <- data.frame(x = coords[, 1], y = coords[, 2], z = vals)
#Define as sf
pts <- st_as_sf(df2, coords = c("x", "y"), crs = NA)


#leave one out for loop with rmse validation
p_vals = c(1, 1.5, 2, 2.5, 3, 4, 5) #I assume theres greater value in miniscule differences between 1 and 3 than 4 and 5

  #score
  rmse <- function(e) {
      sqrt(mean(e^2, na.rm = TRUE)) 
  }

#function that loops for leave-one-out evaluation
loovcf <- function(pts, p_vals = c(1, 2, 3, 4), z_col = "z", func = c("idw", "mba.surf")) { #maybe sane default for p? check with prof

z_vec <- st_drop_geometry(pts)$z
scores <- numeric(length(p_vals))  # one RMSE per p

if (func == "idw") {
  for (j in seq_along(p_vals)) {
    p <- p_vals[j]
    cat("for p =", p, "...\n")
    errs <- numeric(nrow(pts))

    for (i in seq_len(nrow(pts))) {
      train <- pts[-i, ]     # all except i
      test  <- pts[i, ]    
      pred  <- idw(
        formula   = z ~ 1,
        locations = train,
        newdata   = test,
        idp       = p
      )$var1.pred
      errs[i] <- pred - z_vec[i]
    }

    scores[j] <- rmse(errs)
  }

}else {

}

cv_table <- data.frame(p = p_vals, RMSE = scores)
cv_table <- cv_table[order(cv_table$RMSE), ] # sort results by lowest error

return(cv_table)
} 

results <- loovcf(pts, p_vals, func = "idw")
cat("Lowest value: ", results[1,]$p,results[1,]$RMSE)

#! b) Work on the Wolfcamp dataset. Interpolate Wolfcamp data on a 20x20 grid (built from the coordinates) using IDW and multilevel bi-splines (function mba.surf from MBA) and compare the results using a proper score
# Hint: divide data into training and testing sets

# library(geoR)
# data(wolfcamp)
# df <- wolfcamp
# coords <- df$coords 
# vals <- df$data
# df2 <- data.frame(x = coords[, 1], y = coords[, 2], z = vals)
# #Define as sf
# library(sf)
# pts <- st_as_sf(df2, coords = c("x", "y"), crs = NA)

library(gstat)
library(sf)
library(MBA)
library(fields)


make_grid <- function(pts, n = 20) {
  bb <- st_bbox(pts)
  gx <- seq(bb["xmin"], bb["xmax"], length.out = n)
  gy <- seq(bb["ymin"], bb["ymax"], length.out = n)
  expand.grid(x = gx, y = gy)
}

grid20 <- make_grid(pts, n = 20)

# Ensure grid20 is a SpatialPointsDataFrame or sf object
grid20_sf <- st_as_sf(grid20, coords = c("x", "y"), crs = st_crs(pts))


idw_res <- idw(z ~ 1, pts, newdata = grid20_sf, idp = results[1,]$p)
grid20$z_idw <- idw_res$var1.pred


p1<-ggplot(grid20, aes(x = x, y = y, fill = z_idw)) +
  geom_tile() +
  scale_fill_viridis_c() +  # Use a nice color scale
  theme_minimal() +
  labs(title = "IDW Interpolation", fill = "Z-Value IDW") +
  theme(axis.title = element_blank(), axis.text = element_blank())

loovcfmod <- function(pts,
                      p_vals = c(1, 2, 3, 4),
                      n_vals = c(5, 10, 15, 20),
                      m_vals = c(2, 3, 4),
                      h_vals = c(4, 6, 8),
                      extend = "extend",
                      z_col = "z",
                      func = "mba.surf" ) {#c("idw", "mba.surf") #maybe sane default for p? check with prof

z_vec <- st_drop_geometry(pts)$z
scores <- numeric(length(p_vals))  # one RMSE per p

if (func == "idw") {
  for (j in seq_along(p_vals)) {
    p <- p_vals[j]
    cat("for p =", p, "...\n")
    errs <- numeric(nrow(pts))

    for (i in seq_len(nrow(pts))) {
      train <- pts[-i, ]     # all except i
      test  <- pts[i, ]    
      pred  <- idw(
        formula   = z ~ 1,
        locations = train,
        newdata   = test,
        idp       = p
      )$var1.pred
      errs[i] <- pred - z_vec[i]
    }

    scores[j] <- rmse(errs)
  }

  cv_table <- data.frame(p = p_vals, RMSE = scores)
  cv_table <- cv_table[order(cv_table$RMSE), ] # sort results by lowest error

  return(cv_table)

  } else {
    # Parameter grid for MBA
    param_grid <- expand.grid(n = n_vals, m = m_vals, h = h_vals, stringsAsFactors = FALSE)
    param_grid$RMSE <- NA_real_
    param_grid$fold_rmse <- vector("list", nrow(param_grid))

    for (j in seq_len(nrow(param_grid))) {
      n <- param_grid$n[j]
      m <- param_grid$m[j]
      h <- param_grid$h[j]

      errs <- numeric(nrow(pts))  # one error per point

      for (i in seq_len(nrow(pts))) {
        train <- pts[-i, ]     # all except i
        test  <- pts[i, , drop = FALSE]  # just the i-th point

        # Fit the MBA model
        mba_model <- mba.surf(
          xyz = cbind(st_coordinates(train), st_drop_geometry(train)$z),
          no.X = n, no.Y = n, m = m, h = h, extend = TRUE
        )

        # Predict for the test point using the fitted model
        pred_res <- mba.points(convert_mba_grid_to_xyz(mba_model), data.frame(x = st_coordinates(test)[,1], y = st_coordinates(test)[,2]))

        # Extract the prediction (if it exists)
        # if (is.data.frame(pred_res) && "z" %in% names(pred_res)) {
        pred <- pred_res$xyz.est[3] # this corresponds to the value z
        # } else {
        #   pred <- NA_real_
        # }

        # Calculate the error for the i-th point
        errs[i] <- pred - z_vec[i]
      }

      # Store the RMSE for this combination of n, m, h
      param_grid$RMSE[j] <- rmse(errs)
      # param_grid$fold_rmse[[j]] <- errs
      param_grid <- param_grid[order(param_grid$RMSE), ]
    }
    return(param_grid)
  }

}
# 4) Extra: mba_model-to-dataframe helper, why doesnt MBA come with this already?
convert_mba_grid_to_xyz <- function(mba_model) {
  est <- mba_model$xyz.est
  if (is.null(est$x) || is.null(est$y) || is.null(est$z)) {
    stop("mba_model$xyz.est missing")
  }
  grid <- expand.grid(x = est$x, y = est$y, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  zvec <- as.vector(est$z)
  if (length(zvec) != nrow(grid)) stop("z length and x,y rows differ")
  xyz_df <- data.frame(x = grid$x, y = grid$y, z = zvec)
  # drop NA predictions
  xyz_df <- xyz_df[!is.na(xyz_df$z), , drop = FALSE]
  rownames(xyz_df) <- NULL
  xyz_df
}


results_mba <- loovcfmod(pts, n_vals = c(20), func = "mba.surf")
 

mba_res <- mba.surf(xyz = cbind(st_coordinates(pts), st_drop_geometry(pts)$z),
                    no.X = 20, no.Y = 20, m =  results_mba$m[1], h =  results_mba$h[1], extend = TRUE)

mba_df <- expand.grid(x = mba_res$xyz.est$x, y = mba_res$xyz.est$y)
mba_df$z_mba <- as.vector(mba_res$xyz.est$z)

p2<-ggplot(mba_df, aes(x = x, y = y, fill = z_mba)) +
  geom_tile() +
  scale_fill_viridis_c() +  # Use a nice color scale
  theme_minimal() +
  labs(title = "MBA Surface", fill = "Z-Value MBA") +
  theme(axis.title = element_blank(), axis.text = element_blank())

library(cowplot)
plot_grid(p1, p2, ncol = 2)


#! scratchpad


# quick time complexity calc
# source('timeestimation.r', chdir = TRUE)
# res <- estimate_runtime_2loops(n = 4, m = 85, steps_outer = 4, steps_inner = 85)
# big_o_2loops()


# source("gd_to_sf.r")
# gd_to_sf(df, value_name = df$data)
# Required libs (you already loaded them in your session)

# best_p <- best_idw$p
# best_mba_combo <- mba_res[1, ]
# # per-fold RMSE vectors:
# idw_best_fold_rmse <- idw_res$fold_rmse[as.character(best_p), ]   # row from matrix
# mba_best_fold_rmse <- unlist(best_mba_combo$fold_rmse)
# cat("Per-fold RMSEs (IDW vs MBA):\n")
# print(data.frame(fold = 1:max(fold_ids),
#                  idw = idw_best_fold_rmse,
#                  mba = mba_best_fold_rmse))

#  paired t-test on fold RMSEs 
# if (all(!is.na(idw_best_fold_rmse)) && all(!is.na(mba_best_fold_rmse))) {
#   cat("Paired t-test on fold RMSEs:\n")
#   print(t.test(idw_best_fold_rmse, mba_best_fold_rmse, paired = TRUE))
# }
