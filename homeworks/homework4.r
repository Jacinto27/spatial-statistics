# Homework 4
## Using the available data and the interpolation grid:
### assumed to be raindata
load("class files/data/RainData.RData")
library(geoR)
library(sf)
##ls()
##[1] "Bas"             "Basquota32"      "Cal"             "Calabriaquota32"
##[5] "dati"            "grid.int"
head(dati)
assigned_years <- c(1976, 1979, 1982)
colnames(Calabriaquota32)=colnames(Basquota32)
grid=rbind(Basquota32[,c(2:4)],Calabriaquota32[,c(2:4)])

## a) Interpolate the dependent variable on the estimation grid using kriging in a maximum likelihood framework. (choose variogram and trend, use the available covariates)
dati_1976 <- as.geodata(dati[dati$anno==1976,],coords.col=6:7,data.col=10,covar.col=8,covar.names="quota")
dati_1979 <- as.geodata(dati[dati$anno==1979,],coords.col=6:7,data.col=10,covar.col=8,covar.names="quota")
dati_1982 <- as.geodata(dati[dati$anno==1982,],coords.col=6:7,data.col=10,covar.col=8,covar.names="quota")

years_list <- list("1976" = dati_1976, "1979" = dati_1979, "1982" = dati_1982)

for (year in names(years_list)) {
  dev.new()  # Opens a new graphics device
  par(mfrow = c(2, 2))
  plot(years_list[[year]], name=paste("Data for", year))
}

## There seems to be NO correlation with the x and y coordinates, what then?
# should i test normality fo residuals even with no correlation? doesn't make sense 

for (year in names(years_list)) {
cat("year=", year, "p value=", format(shapiro.test(years_list[[year]]$data)$p.value, scientific = FALSE), "\n")
}

# Output
## year= 1976 p value= 0.000001527713 
## year= 1979 p value= 0.0000008326988 
## year= 1982 p value= 0.0003342081 

## values are non normal, exponetnial? trend? not sure what to do

for (year in names(years_list)) {
cat("year=", year, "p value (residuals) =", 
format(
  shapiro.test(
    residuals(
      lm(years_list[[year]]$data ~ years_list[[year]]$coords[, 1] + years_list[[year]]$coords[, 2]))
      )$p.value, scientific = FALSE), "\n")
}
#doesn't fix anything

for (year in names(years_list)) {
cat("year=", year, "p value (residuals) =", 
format(
  shapiro.test(
    residuals(
      lm(years_list[[year]]$data ~ years_list[[year]]$covariate[, 1] ))
      )$p.value, scientific = FALSE), "\n")
}
#doesn't fix anything either perhaps boxcox?

for (year in names(years_list)) {
cat("year=", year, "p value (residuals) =", 
format(
  shapiro.test(
      log(years_list[[year]]$data)
      )$p.value, scientific = FALSE), "\n")
}
## Output:
## year= 1976 p value (residuals) = 0.4267471 is gaussian when log 
## year= 1979 p value (residuals) = 0.1677097 is gaussian when log 
## year= 1982 p value (residuals) = 0.0117682 is not even when transformed

## for the last year maybe i should try sqrt o algo 

for (year in names(years_list)) {
cat("year=", year, "p value (residuals) =", 
format(
  shapiro.test(
      sqrt(years_list[[year]]$data)
      )$p.value, scientific = FALSE), "\n")
}
# conclusion: data isn't gaussian? only when transformed for first and second year
#when transforming to sqrt to the last dataset becomes gaussian 
## year= 1976 p value (residuals) = 0.004546831 
## year= 1979 p value (residuals) = 0.04664353 
## year= 1982 p value (residuals) = 0.4461693 

lambda_mapping <- c("1976" = 0, "1979" = 0, "1982" = 0.5)

for (year in names(years_list)) {
  dev.new() # Opens a new graphics device
  par(mfrow = c(1, 1))
  vv <- variog(years_list[[year]], max.dist = 250000, trend = "1st", lambda = lambda_mapping[[year]])
  plot(vv, main = paste("Data for", year))
}
# Data looks good enough for first few years with log transfor, for the last
# year the log doesn't make it gaussian but hte sqrt does, however, the data for
# the last year still needs to be adjusted for the quota covariate to make it
# stationary, what should i do? trend all years? keep them as separate
# estimations?

# fuction that does cross validation 
run_cv_for_year <- function(geodata_year, max_distance = 250000, kappa_values = c(0.5, 1.0, 1.5)) {
  # Step 1: Compute the variogram
  vv <- variog(geodata_year, max.dist = max_distance, trend = "1st")
  
  # Step 2: Fit the Exponential and Matérn models
  vv_fit_exp <- variofit(vv, cov.model = "exponential")
  best_model <- "Exponential"
  best_error <- Inf
  best_params <- NULL
  
  # Step 3: Cross-validation for Exponential model
  geostat_exp <- likfit(geodata_year, model = vv_fit_exp, trend = "cte", lik.method = "ML")
  xvalid_exp <- xvalid(geodata = geodata_year, model = geostat_exp, locations.xvalid = "all")
  exp_error <- mean((xvalid_exp$observed - xvalid_exp$predicted)^2)  # Calculate RMSE or other error metric
  
  # If Exponential is the best so far, save the result
  if (exp_error < best_error) {
    best_error <- exp_error
    best_model <- "Exponential"
    best_params <- vv_fit_exp
  }
  
  # Step 4: Cross-validation for Matérn model (with different kappa values)
  for (kappa in kappa_values) {
    vv_fit_matern <- variofit(vv, cov.model = "matern", kappa = kappa)
    
    geostat_matern <- likfit(geodata_year, model = vv_fit_matern, trend = "cte", lik.method = "ML")
    xvalid_matern <- xvalid(geodata = geodata_year, model = geostat_matern, locations.xvalid = "all")
    matern_error <- mean((xvalid_matern$observed - xvalid_matern$predicted)^2)  # RMSE for Matérn
    
    # If Matérn is the best, update
    if (matern_error < best_error) {
      best_error <- matern_error
      best_model <- "Matérn"
      best_params <- vv_fit_matern
    }
  }
   # Step 4: Cross-validation for Spherical model
  vv_fit_sph <- variofit(vv, cov.model = "spherical")
  geostat_sph <- likfit(geodata_year, model = vv_fit_sph, trend = "cte", lik.method = "ML")
  xvalid_sph <- xvalid(geodata = geodata_year, model = geostat_sph, locations.xvalid = "all")
  sph_error <- mean((xvalid_sph$observed - xvalid_sph$predicted)^2)  # RMSE for Spherical
  
  # If Spherical is the best, update
  if (sph_error < best_error) {
    best_error <- sph_error
    best_model <- "Spherical"
    best_params <- vv_fit_sph
  } 
  # Return the best model, its parameters, and the error
  return(list(
    best_model = best_model,
    best_params = best_params,
    best_error = best_error
  ))
}


# Main loop for processing each year
cv_results <- list()  # Store results for each year
##years_list <- list()  # Your list of data for each year (geodata)

for (year in names(years_list)) {
  geodata_year <- years_list[[year]]
  
  # Run the cross-validation for the current year
  result <- run_cv_for_year(geodata_year)
  
  # Store the results for this year
  cv_results[[year]] <- result
  
  # Output the best model and error for this year
  cat("Best model for", year, "is", result$best_model, "with error", result$best_error, "\n")
}
## b) given a) build pointwise conf  
## An Rmarkdown with code and a few comments will be optimal (send compiled version too).

##scratchpad
## data(wolfcamp)
## head(wolfcamp)
## 
## coords <- df$coords 
## vals <- df$data
## df2 <- data.frame(x = coords[, 1], y = coords[, 2], z = vals)
## pts <- st_as_sf(df2, coords = c("x", "y"), crs = NA)