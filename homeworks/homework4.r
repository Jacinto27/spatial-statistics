# Homework 4
## Using the available data and the interpolation grid:
### assumed to be raindata
load("class files/data/RainData.RData")
library(geoR)
library(sf)
## ls()
## [1] "Bas"             "Basquota32"      "Cal"             "Calabriaquota32"
## [5] "dati"            "grid.int"
head(dati)
assigned_years <- c(1976, 1979, 1982)
colnames(Calabriaquota32) <- colnames(Basquota32)
grid <- rbind(Basquota32[, c(2:4)], Calabriaquota32[, c(2:4)])
grid[,1:2]<-grid[,1:2]/1000
dati[,6:7]<-dati[,6:7]/1000
## a) Interpolate the dependent variable on the estimation grid using kriging in a maximum likelihood framework. (choose variogram and trend, use the available covariates)
dati_1976 <- as.geodata(dati[dati$anno == 1976, ], coords.col = 7:6, data.col = 10, covar.col = 8, covar.names = "quota")
dati_1979 <- as.geodata(dati[dati$anno == 1979, ], coords.col = 6:7, data.col = 10, covar.col = 8, covar.names = "quota")
dati_1982 <- as.geodata(dati[dati$anno == 1982, ], coords.col = 6:7, data.col = 10, covar.col = 8, covar.names = "quota")
grid_geodata <- as.geodata(grid,coords.col=1:2,covar.col=3)
years_list <- list("1976" = dati_1976, "1979" = dati_1979, "1982" = dati_1982)
#lambda_mapping <- c("1976" = 0, "1979" = 0, "1982" = 0.5)
#Replicating kriging with geoR
plot(dati_1976)
vv <-  variog(years_list[['1976']], lambda = lambda_mapping[['1976']],trend=~ quota,max.dist = 100)
plot(vv, type='b')
vv_fit_matern <- variofit(vv, cov.model = "exponential", max.dist = 100)
vv.exp <- likfit(dati_1976,cov.model = "exponential",ini.cov.pars = vv_fit_matern$cov.pars,lambda = lambda_mapping[['1976']])
lines(vv.exp,col=2)
    xvalid_matern <- xvalid(geodata = geodata_year, model = geostat_matern, locations.xvalid = "all")
    matern_error <- mean((xvalid_matern$error)^2) # RMSE for Matérn
kk.c<-krige.control(obj.model = vv.exp,lambda = lambda_mapping[['1976']])
krig1<-krige.conv(dati_1976,locations = as.matrix(grid),krige = kk.c)
#Plotting datka
 # Convert kriging result to a SpatialGridDataFre
 #plot(grid)
 library(ggplot2)
 ggplot(grid, aes(x = XUTM, y = YUTM, color = krig1$predict)) +
   geom_point() +
   scale_color_viridis_c() +
   theme_minimal() +
   labs(title = "Kriging Estimates", x = "XUTM", y = "YUTM", color = "Predicted Value")
#todo: use spplot by changing it to sp
for (year in names(years_list)) {
  dev.new() # Opens a new graphics device
  par(mfrow = c(2, 2))
  plot(years_list[[year]], name = paste("Data for", year))
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

## values are non normal

for (year in names(years_list)) {
  cat(
    "year=", year, "p value (residuals) =",
    format(
      shapiro.test(
        residuals(
          lm(years_list[[year]]$data ~ years_list[[year]]$coords[, 1] + years_list[[year]]$coords[, 2])
        )
      )$p.value,
      scientific = FALSE
    ), "\n"
  )
}
# residuals
# doesn't fix anything

for (year in names(years_list)) {
  cat(
    "year=", year, "p value (residuals) =",
    format(
      shapiro.test(
        residuals(
          lm(years_list[[year]]$data ~ years_list[[year]]$covariate[, 1])
        )
      )$p.value,
      scientific = FALSE
    ), "\n"
  )
}
# doesn't fix anything either perhaps boxcox?

for (year in names(years_list)) {
  cat(
    "year=", year, "p value (residuals) =",
    format(
      shapiro.test(
        residuals(
          lm(years_list[[year]]$data ~ years_list[[year]]$covariate[, 1])
        )
      )$p.value,
      scientific = FALSE
    ), "\n"
  )
}
# residuals adjusting for x and y and for height
# doesn't fix anything

for (year in names(years_list)) {
  cat(
    "year=", year, "p value (residuals) =",
    format(
      shapiro.test(
        log(years_list[[year]]$data)
      )$p.value,
      scientific = FALSE
    ), "\n"
  )
}
# log transform:
## Output:
## year= 1976 p value (residuals) = 0.4267471 is gaussian when log
## year= 1979 p value (residuals) = 0.1677097 is gaussian when log
## year= 1982 p value (residuals) = 0.0117682 is not even when transformed

## for the last year maybe i should try sqrt or something

for (year in names(years_list)) {
  cat(
    "year=", year, "p value (residuals) =",
    format(
      shapiro.test(
        sqrt(years_list[[year]]$data)
      )$p.value,
      scientific = FALSE
    ), "\n"
  )
}
# conclusion: data is guassian only when log-transformed for first and second year
# when transforming to sqrt to the last dataset becomes gaussian
## year= 1976 p value (residuals) = 0.004546831
## year= 1979 p value (residuals) = 0.04664353
## year= 1982 p value (residuals) = 0.4461693

for (year in names(years_list)) {
  dev.new() # Opens a new graphics device
  par(mfrow = c(1, 1))
  vv <- variog(years_list[[year]],  lambda = lambda_mapping[[year]])
  plot(vv, main = paste("Data for", year))
}
# Data looks good enough for first few years with log transfor, for the last
# year the log doesn't make it gaussian but hte sqrt does, however, the data for
# the last year still needs to be adjusted for the quota covariate to make it
# stationary, what should i do?
# Data looks good enough for first few years with log transfor, for the last
# year the log doesn't make it gaussian but hte sqrt does, however, the data for
# the last year still needs to be adjusted for the quota covariate to make it
# stationary, what should i do? trend all years? keep them as separate
# estimations?
lambda_mapping <- list()
for (year in names(years_list)) {
  res <- boxcoxfit(years_list[[year]]$data)
  lambda_mapping[year] <- res$lambda  
}
# Data looks good enough for first few years with log transfor, for the last
for (year in names(years_list)) {
  dev.new() # Opens a new graphics device
  par(mfrow = c(1, 1))
  vv <- variog4(years_list[[year]], lambda = lambda_mapping[[year]])
  plot(vv, main = paste("Data for", year))
}
for (year in names(years_list)) {
  dev.new() # Opens a new graphics device
  par(mfrow = c(1, 1))
  vv <- variog4(years_list[[year]], lambda = lambda_mapping[[year]], trend='1st')
  plot(vv, main = paste("Data for", year))
}
for (year in names(years_list)) {
  dev.new() # Opens a new graphics device
  par(mfrow = c(1, 1))
  vv <- variog4(years_list[[year]], lambda = lambda_mapping[[year]],trend=~ quota)
  plot(vv, main = paste("Data for", year))
}

#todo: wrap in for loop later
  vv_1976 <- variog(years_list[['1976']], lambda = lambda_mapping[['1976']],trend=~ quota,max.dist = 100)
  vv_1979 <- variog(years_list[['1979']], lambda = lambda_mapping[['1979']],trend='1st',max.dist = 100)
  vv_1982 <- variog(years_list[['1982']], lambda = lambda_mapping[['1982']],trend=~ quota,max.dist = 100)
variog_mapping <- list()

variog_mapping[["1976"]] = vv_1976 
 variog_mapping[["1979"]] = vv_1979 
 variog_mapping[["1982"]] = vv_1982
# year the log doesn't make it gaussian but hte sqrt does, however, the data for
# the last year still needs to be adjusted for the quota covariate to make it
# stationary, what should i do?
# fuction that does cross validation
run_cv_for_year <- function(variog,geodata_year, max_distance = 100, kappa_values = c(0.1, 0.3, 0.5, 0.7, 1.0, 1.5), lambda_mapping) {
  # Step 1: Compute the variogram
  ##  # Step 2: Fit the Exponential and Matérn models
  trend<- as.formula('~ quota')
  ##  vv_fit_exp <- variofit(vv, cov.model = "exponential")
  ##  best_model <- "Exponential"
  best_error <- Inf
  best_params <- NULL
  ##
  ##  # Step 3: Cross-validation for Exponential model
  ##  geostat_exp <- likfit(geodata_year, model = vv_fit_exp, trend = "cte", lik.method = "ML")
  ##  xvalid_exp <- xvalid(geodata = geodata_year, model = geostat_exp, locations.xvalid = "all")
  ##  exp_error <- mean((xvalid_exp$observed - xvalid_exp$predicted)^2)  # Calculate RMSE or other error metric
  ##
  ##  # If Exponential is the best so far, save the result
  ##  if (exp_error < best_error) {
  ##    best_error <- exp_error
  ##    best_model <- "Exponential"
  ##    best_params <- vv_fit_exp
  ##  }
  # Step 4: Cross-validation for Matérn model (with different kappa values)
  for (kappa in kappa_values) { # kappa of 0.5 for mattern is just exponential
    vv_fit_matern <- variofit(variog, cov.model = "matern", max.dist = max_distance, kappa = kappa)

    geostat_matern <- likfit(geodata_year,
      ## model = vv_fit_matern,
       trend = trend,
      #lik.method = "ML",
      ini.cov.pars = c(vv_fit_matern$cov.pars[1], vv_fit_matern$cov.pars[2]),
      nugget = vv_fit_matern$nugget,
      lambda = lambda_mapping
    )

    xvalid_matern <- xvalid(geodata = geodata_year, model = geostat_matern, locations.xvalid = "all")
    matern_error <- mean((xvalid_matern$error)^2) # RMSE for Matérn
    matern_cv <-sd(xvalid_matern $std.error) / mean(xvalid_matern$std.error) 
    # If Matérn is the best, update
    if (matern_error < best_error) {
      best_error <- matern_error
      best_model <- "matern"
      best_params <- geostat_matern
    }
  }
  # Step 4: Cross-validation for Spherical model
  vv_fit_sph <- variofit(vv, cov.model = "spherical", max.dist = max_distance)

  geostat_sph <- likfit(geodata_year,
    ## model = vv_fit_sph,
     trend = trend,
    #lik.method = "ML",
    ini.cov.pars = c(vv_fit_sph$cov.pars[1], vv_fit_sph$cov.pars[2]),
    nugget = vv_fit_sph$nugget, lambda = lambda_mapping
  )

  xvalid_sph <- xvalid(geodata = geodata_year, model = geostat_sph, locations.xvalid = "all")
  sph_error <- mean((xvalid_sph$error)^2) # RMSE for Spherical

  # If Spherical is the best, update
  if (sph_error < best_error) {
    best_error <- sph_error
    best_model <- "spherical"
    best_params <- geostat_sph 
  }
  # Return the best model, its parameters, and the error
  return(list(
    best_model = best_model,
    best_params = best_params,
    best_error = best_error
  ))
  # Step 5: Cross-validation for Gaussian model
  vv_fit_gaus <- variofit(vv, cov.model = "gaussian", max.dist = max_distance)

  geostat_gaus <- likfit(geodata_year,
    ## model = vv_fit_sph,
    trend = trend,
    #lik.method = "ML",
    ini.cov.pars = c(vv_fit_gaus$cov.pars[1], vv_fit_gaus$cov.pars[2]),
    nugget = vv_fit_gaus$nugget, lambda = lambda_mapping
  )

  xvalid_gaus <- xvalid(geodata = geodata_year, model = geostat_gaus, locations.xvalid = "all")
  gaus_error <- mean((xvalid_gaus$error)^2) # RMSE for gaus

  # If Spherical is the best, update
  if (gaus_error < best_error) {
    best_error <- gaus_error
    best_model <- "gaussian"
    best_params <- geostat_gaus
  }
  # Return the best model, its parameters, and the error
  return(list(
    best_model = best_model,
    best_params = best_params,
    best_error = best_error
  ))
}
# Main loop for processing each year
cv_results <- list() # Store results for each year
## years_list <- list()  # Your list of data for each year (geodata)

for (year in names(years_list)) {
  geodata_year <- years_list[[year]]

  # Run the cross-validation for the current year
  result <- run_cv_for_year(geodata_year,variog=variog_mapping[[year]], lambda_mapping = lambda_mapping[[year]])

  # Store the results for this year
  cv_results[[year]] <- result

  # Output the best model and error for this year
  cat("Best model for", year, "is", result$best_model, "with error", result$best_error, "\n")
}

library(gstat)
# estimating the grid
krig_res <- list()
kriging_dataframe <- list()
for (year in names(years_list)) {
   trend_d <- trend.spatial('cte',geodata = years_list[[year]] )
   trend_l <- trend.spatial('cte',geodata=grid_geodata)

  kk_control <- krige.control(
    #type.krige = "OK",
    trend.d = trend_d,
    trend.l =trend_l,
    obj.model = cv_results[[year]]$best_params,
    cov.model = cv_results[[year]]$best_params$cov.model,
    cov.pars = cv_results[[year]]$best_params$cov.pars,
    nugget = cv_results[[year]]$best_params$nugget,
    lambda = lambda_mapping[[year]]
  )

  krig_res[[year]] <- krige.conv(years_list[[year]], locations = grid[,1:2], krige = kk_control)
  ## pred_krig <- krige(formula = "data~1", years_list[[year]]$data, model = cv_results[[year]]$best_params$cov.model,newdata=grid)
  kriging_dataframe[[year]] <-data.frame(
    year=year,
    XUTM=grid$XUTM,
    YUTM=grid$YUTM,
    Quota=grid$Quota,
    kriging_pred=krig_res[[year]]$predict,
    kiring_variance=krig_res[[year]]$krige.var,
    model=cv_results[[year]]$best_model
  )
  plotting_df <- data.frame(
    XUTM <- grid$XUTM,
    YUTM <- grid$YUTM,
    kriging_pred<- krig_res[[year]]$predict
  )
 
}
##for (year in names(years_list)) {
##  dev.new() # Opens a new graphics device
##  par(mfrow = c(1, 1))
##ggplot(grid, aes(x = XUTM, y = YUTM, color = krig_res[[year]]$predict)) +
##  geom_point() +
##  scale_color_viridis_c() +
##  theme_minimal() +
##  labs(title = "Kriging Estimates", x = "XUTM", y = "YUTM", color = "Predicted Value")
##}
library(ggplot2)
ggplot(grid, aes(x = XUTM, y = YUTM, color = krig_res[["1982"]]$predict)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Kriging Estimates", x = "XUTM", y = "YUTM", color = "Predicted Value")
# Data looks good enough for first few years with log transfor, for the last
# year the log doesn't make it gaussian but hte sqrt does, however, the data for
# the last year still needs to be adjusted for the quota covariate to make it
# stationary, what should i do?

# c) map results, including confidence intervals, and comment on the output.


## scratchpad
## data(wolfcamp)
## head(wolfcamp)
##
## coords <- df$coords
## vals <- df$data
## df2 <- data.frame(x = coords[, 1], y = coords[, 2], z = vals)
## pts <- st_as_sf(df2, coords = c("x", "y"), crs = NA)

## for (year in names(years_list)) {
vv <- variog(years_list[["1979"]], max.dist = 250000, trend = "1st", lambda = 0)

plot(vv,
  type = "l", lwd = 2, main = paste("Variogram (Log-Rain) - 1979"),
  xlab = "Distance (m)"
)

lines(cv_results[["1979"]]$best_param, col = "blue", lwd = 2)

# Convert kriging result to a SpatialGridDataFrame
plot(grid)
library(ggplot2)
ggplot(grid, aes(x = XUTM, y = YUTM, color = krig_res[["1982"]]$predict)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Kriging Estimates", x = "XUTM", y = "YUTM", color = "Predicted Value")

new_sf <- st_as_sf(dati[dati$anno == 1976, ],coords=c('XUTM','YUTM'))
plot(new_sf)

## library(interp)

## xy <- interp(years_list[['1976']]$XUTM, years_list[['1976']]$YUTM, years_list[['1976']]$totanno)
## surface(xy)
## lines(Cal[, 1] / 1000, Cal[, 2] / 1000, lwd = 2, col = "white")
## lines(Bas[, 1] / 1000, Bas[, 2] / 1000, lwd = 2, col = "white")

dati.sp <- dati
sp::coordinates(dati.sp) <- ~XUTM + YUTM
sp::gridded(grid)=~XUTM+YUTM
geodat_grid <- as.geodata(grid,coords.col=c("XUTM","YUTM"),covar.col="Quota")
pred_krig <- gstat::krige( formula = totanno ~ 1, dati[dati$anno==1982] ,locations=geodat_grid$coords, model = cv_results[['1982']]$best_params$cov.model, newdata=dati.sp)

dati_1976$coords<-dati_1976$coords/1000
grid[,1:2]<-grid[,1:2]/1000
vv<-variog(dati_1976,trend='1st',max.dist = 250,lambda=0)
plot(vv)
vv.exp<-likfit(dati_1976,trend='1st',cov.model='exponential',ini.cov.pars = c(0.08,50),lambda=0)
lines(vv.exp,col=2)
kk.c<-krige.control(trend.d = "1st",trend.l = "1st",obj.model = vv.exp,lambda=0)
krig1<-krige.conv(dati_1976,locations = grid[,1:2],krige = kk.c)
