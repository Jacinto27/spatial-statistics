# a) Build an R function to choose the parameter p in the inverse distance weighting estimator. Hint: use the idw implemented in gstat.
library(gstat)
library(sf)
library(geoR) # wolfcamp data, also in MBA package

data(wolfcamp)
df <- wolfcamp
#structure(df)
# remove old geodata structure and turn into a dataframe NOTE: MBA package has wolfcamp as dataframe, step not needed if using
coords <- df$coords 
vals <- df$data
df2 <- data.frame(x = coords[, 1], y = coords[, 2], z = vals)
#Define as sf with the nomenclature from class: x & y are coordinates, z are values
pts <- st_as_sf(df2, coords = c("x", "y"), crs = NA)


#leave-one-out loop with rmse validation
p_vals = c(1, 1.5, 2, 2.5, 3, 4, 5) #I assume theres greater value in miniscule differences between 1 and 3 than 4 and 5, ask the prof.

#error score
rmse <- function(e) {
    sqrt(mean(e^2, na.rm = TRUE)) 
}

#function that loops
loovcf <- function(pts, p_vals = c(1, 2, 3, 4), z_col = "z") { #generic defaults defalt for p? check with prof

z_vec <- st_drop_geometry(pts)$z
scores <- numeric(length(p_vals))  # one RMSE per p


for (j in seq_along(p_vals)) { # for each value in the list P
  p <- p_vals[j]
  cat("for p =", p, "...\n")
  errs <- numeric(nrow(pts))

  for (i in seq_len(nrow(pts))) {# run the regression for all values in the data except for 1, and do this until every single row has been removed and replaced 
    train <- pts[-i, ]     # all except i
    test  <- pts[i, ]    
    pred  <- idw(
      formula   = z ~ 1, # is this formula correct?
      locations = train,
      newdata   = test,
      idp       = p
    )$var1.pred
    errs[i] <- pred - z_vec[i] #store each error term into a list for each regression,
  }

  scores[j] <- rmse(errs) # calculate the mean square error of all regression for each p value
}

cv_table <- data.frame(p = p_vals, RMSE = scores) #store in table and 
cv_table <- cv_table[order(cv_table$RMSE), ] # sort results by lowest error

return(cv_table)
} 

results <- loovcf(pts, p_vals)
## lowest rmse for p:4 for wolfcamp data

# b) Work on the Wolfcamp dataset. Interpolate Wolfcamp data on a 20x20 grid (built from the coordinates) using IDW and multilevel bi-splines (function mba.surf from MBA) and compare the results using a proper score
#  Hint: divide data into training and testing sets

library(geoR)


# scatchpad


# quick time complexity calc
# source('timeestimation.r', chdir = TRUE)
# res <- estimate_runtime_2loops(n = 4, m = 85, steps_outer = 4, steps_inner = 85)
# big_o_2loops()


# source("gd_to_sf.r")
# gd_to_sf(df, value_name = df$data)
