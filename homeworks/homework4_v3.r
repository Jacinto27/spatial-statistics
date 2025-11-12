# homework 4 fixed

load("class files/data/RainData.RData")
library(geoR)
library(sf)
## ls()
## [1] "Bas"             "Basquota32"      "Cal"             "Calabriaquota32"
## [5] "dati"            "grid.int"
#head(dati)
assigned_years <- c(1975, 1979, 1982)
colnames(Calabriaquota32) <- colnames(Basquota32)
grid <- rbind(Basquota32[, c(2:4)], Calabriaquota32[, c(2:4)])
grid[, 1:2] <- grid[, 1:2] / 1000
dati[, 6:7] <- dati[, 6:7] / 1000
## a) Interpolate the dependent variable on the estimation grid using kriging in a maximum likelihood framework. (choose variogram and trend, use the available covariates)
dati_1975 <- as.geodata(dati[dati$anno == 1975, ], coords.col = 7:6, data.col = 10, covar.col = 8, covar.names = "quota")
dati_1979 <- as.geodata(dati[dati$anno == 1979, ], coords.col = 7:6, data.col = 10, covar.col = 8, covar.names = "quota")
dati_1982 <- as.geodata(dati[dati$anno == 1982, ], coords.col = 7:6, data.col = 10, covar.col = 8, covar.names = "quota")
grid_geodata <- as.geodata(grid, coords.col = 1:2, covar.col = 3)
years_list <- list("1975" = dati_1975, "1979" = dati_1979, "1982" = dati_1982)
# lambda_mapping <- c("1975" = 0, "1979" = 0, "1982" = 0.5)
# Replicating kriging with geoR
lambda_mapping <- c()
for (year in names(years_list)) {
  res <- boxcoxfit(years_list[[year]]$data)
  lambda_mapping[year] <- res$lambda
}
variog_mapping <- c()
# for (year in names(years_list)) {
# ifelse(year=='1979',trend<-'1st',trend <- as.formula("~ quota"))
# vv <-  variog(years_list[[year]], lambda = lambda_mapping[year],trend=tren,max.dist = 100)
# variog_mapping[[year]]<-vv
# }

vv1_raw <- variog(years_list[["1975"]], lambda = lambda_mapping["1975"], trend = ~quota)
vv2_raw <- variog(years_list[["1979"]], lambda = lambda_mapping["1979"], trend = ~quota)
vv3_raw <- variog(years_list[["1982"]], lambda = lambda_mapping["1982"], trend = ~quota)

plot(vv1_raw, main = "Variog 1975", type = "b")
plot(vv2_raw, main = "Variog 1979", type = "b")
plot(vv3_raw, main = "Variog 1982", type = "b")

vv1_4 <- variog4(years_list[["1975"]], lambda = lambda_mapping["1975"], trend = ~quota)
vv2_4 <- variog4(years_list[["1979"]], lambda = lambda_mapping["1979"], trend = ~quota)
vv3_4 <- variog4(years_list[["1982"]], lambda = lambda_mapping["1982"], trend = ~quota)

plot(vv1_4, main = "Variog 1975 Anisotropy", type = "b")
plot(vv2_4, main = "Variog 1979 Anisotropy", type = "b")
plot(vv3_4, main = "Variog 1982 Anisotropy", type = "b")

vv1 <- variog(years_list[["1975"]], lambda = lambda_mapping["1975"], trend = ~quota, max.dist = 100)
vv2 <- variog(years_list[["1979"]], lambda = lambda_mapping["1979"], trend = ~quota, max.dist = 100)
vv3 <- variog(years_list[["1982"]], lambda = lambda_mapping["1982"], trend = ~quota, max.dist = 100)
#variog_mapping[["1975"]] <- vv1
#variog_mapping[["1979"]] <- vv2
#variog_mapping[["1982"]] <- vv3

plot(vv1, main = "Variog 1975 Truncated: 100km", type = "b")
plot(vv2, main = "Variog 1979 Truncated: 100km", type = "b")
plot(vv3, main = "Variog 1982 Truncated: 100km", type = "b")

xvalid_matern <- list()
matern_error <- c()
matern_cv <- c()

for (year in names(years_list)) {
  kappa=3
  vv <- variog(years_list[[year]], lambda = lambda_mapping[year], trend = ~quota, max.dist = 100)
  vv_fit_matern <- variofit(vv, cov.model = "matern", max.dist = 100, kappa = kappa)
  geostat_matern <- likfit(years_list[[year]],
    lambda = lambda_mapping[year],
    cov.model = "matern",
    ini.cov.pars = vv_fit_matern$cov.pars
  )
  xvalid_matern[[year]] <- xvalid(geodata = years_list[[year]], model = geostat_matern, locations.xvalid = "all")
  matern_error[year] <- sqrt(mean((xvalid_matern[[year]]$error^2))) # RMSE for Matérn
  matern_cv[year] <-mean(xvalid_matern[[year]]$std.error^2) 
}

matern_error # RMSE for Matérn
matern_cv 

xvalid_exp <- list()
exp_error <- c()
exp_cv <- c()
kk.a <- list()
krig0 <- list()

for (year in names(years_list)) {
  vv <- variog(years_list[[year]], lambda = lambda_mapping[year], trend = ~quota, max.dist = 100)
  vv_fit_exp <- variofit(vv, cov.model = "exponential", max.dist = 100)
  geostat_exp <- likfit(years_list[[year]],
    lambda = lambda_mapping[year],
    cov.model = "exponential",
    ini.cov.pars = vv_fit_exp$cov.pars
  )
  xvalid_exp[[year]] <- xvalid(geodata = years_list[[year]], model = geostat_exp, locations.xvalid = "all")
  exp_error[year] <- sqrt(mean((xvalid_exp[[year]]$error^2))) # rmse for matérn
  exp_cv[year] <-mean(xvalid_exp[[year]]$std.error^2)  
  kk.a[[year]] <- krige.control(obj.model = geostat_exp, lambda = lambda_mapping[year])
  krig0[[year]] <- krige.conv(years_list[[year]], locations = as.matrix(grid), krige = kk.a[[year]])
  if (year %in% c("1979", "1982")) {
    plot(vv, main = paste("Empirical vs Fitted Variogram -", year)) # base variogram
    lines(vv_fit_exp, col = "red", lwd = 2) # variofit (black)
  }
}

xvalid_sph <- list()
sph_error <- c()
sph_cv <- c()
kk.c <- list()
krig1 <- list()
for (year in names(years_list)) {
  vv <- variog(years_list[[year]], lambda = lambda_mapping[year], trend = ~quota, max.dist = 100)
  vv_fit_sph <- variofit(vv, cov.model = "spherical", max.dist = 100)
  geostat_sph <- likfit(years_list[[year]],
    lambda = lambda_mapping[year],
    cov.model = "spherical",
    ini.cov.pars = vv_fit_sph$cov.pars
  )
  xvalid_sph[[year]] <- xvalid(geodata = years_list[[year]], model = geostat_sph, locations.xvalid = "all")
  sph_error[year] <- sqrt(mean((xvalid_sph[[year]]$error^2))) # RMSE for Matérn
  sph_cv[year] <-mean(xvalid_sph[[year]]$std.error^2) 

  kk.c[[year]] <- krige.control(obj.model = geostat_sph, lambda = lambda_mapping[year])
  krig1[[year]] <- krige.conv(years_list[[year]], locations = as.matrix(grid), krige = kk.c[[year]])
  if (!(year %in% c("1979", "1982"))) {
    plot(vv, main = paste("Empirical vs Fitted Variogram -", year)) # base variogram
    lines(vv_fit_sph, col = "red", lwd = 2) # variofit (black)
  }
}


xvalid_gaus <- list()
gaus_error <- c()
gaus_cv <- c()

for (year in names(years_list)) {
  vv <- variog(years_list[[year]], lambda = lambda_mapping[year], trend = ~quota, max.dist = 100)
  vv_fit_gaus <- variofit(vv, cov.model = "gaussian", max.dist = 100)
  geostat_gaus <- likfit(years_list[[year]],
    lambda = lambda_mapping[year],
    cov.model = "gaussian",
    ini.cov.pars = vv_fit_gaus$cov.pars
  )
  xvalid_gaus[[year]] <- xvalid(geodata = years_list[[year]], model = geostat_gaus, locations.xvalid = "all")
  gaus_error[year] <- sqrt(mean((xvalid_exp[[year]]$error^2))) # rmse for matérn
  gaus_cv[year] <- mean(xvalid_gaus[[year]]$std.error^2) 
}
## kappa_values <- seq(1, 5, length.out = 5)

## for(year in names(years_list)) {
## vv <- variog(years_list[[year]],
## lambda = lambda_mapping[year],
## trend = ~ quota,
## max.dist = 100)

### Initialize a vector to store errors for each kappa
## kappa_rmse <- numeric(length(kappa_values))

## for(i in seq_along(kappa_values)) {
## kappa_val <- kappa_values[i]

### Fit variogram with current kappa
## vv_fit_matern <- variofit(vv, cov.model = "matern",
## max.dist = 100, kappa = kappa_val)

### Fit model
## geostat_matern <- likfit(years_list[[year]],
## lambda = lambda_mapping[year],
## cov.model = "matern",
## ini.cov.pars = vv_fit_matern$cov.pars,
## kappa = kappa_val)

### Cross-validation
## xval <- xvalid(geodata = years_list[[year]],
## model = geostat_matern,
## locations.xvalid = "all")

### Compute RMSE for this kappa
## kappa_rmse[i] <- sqrt(mean((xval$error^2)))
## }

### Find the best kappa for this year (min RMSE)
## best_idx <- which.min(kappa_rmse)
## best_kappa[[year]] <- kappa_values[best_idx]

### Fit again using the best kappa
## vv_fit_matern_best <- variofit(vv, cov.model = "matern",
## max.dist = 100, kappa = best_kappa[[year]])
## geostat_matern_best <- likfit(years_list[[year]],
## lambda = lambda_mapping[year],
## cov.model = "matern",
## ini.cov.pars = vv_fit_matern_best$cov.pars,
## kappa = best_kappa[[year]])

### Store results for the best model
## xvalid_matern[[year]] <- xvalid(geodata = years_list[[year]],
## model = geostat_matern_best,
## locations.xvalid = "all")
## matern_error[year] <- sqrt(mean((xvalid_matern[[year]]$error^2))) # RMSE
## matern_cv[year] <- mean(xvalid_matern[[year]]$std.error / years_list[[year]]$data)
## }

## Fill with na the missing values of gaus 
gaus_error <- gaus_error[match(names(years_list), names(gaus_error))]
gaus_cv    <- gaus_cv[match(names(years_list), names(gaus_cv))]

results <- data.frame(
  year = names(exp_error),
  sph_error = unlist(sph_error),
  sph_cv = unlist(sph_cv),
  exp_error = unlist(exp_error),
  exp_cv = unlist(exp_cv),
  matern_error = unlist(matern_error),
  matern_cv = unlist(matern_cv),
  gaus_error = unlist(gaus_error),
  gaus_cv = unlist(gaus_cv)
)
rownames(results) <- results$year
results$year <- NULL
results
cat("kappa =",kappa)

for (year in rownames(results)) {
  # extract the row of errors for this year
  errors <- results[year, grep("_error$", names(results))]

  # find which model has the lowest error
  best_model <- names(errors)[which.min(errors)]

  # print formatted result
  cat(year, "lowest value:", best_model, "\n")
}
library(sp)
#spplot(krig1[["1975"]])

# confidence interval
alpha <- qnorm(0.975)
kkl1975 <- krig1[["1975"]]$predict - alpha * sqrt(krig1[["1975"]]$krige.var)
kku1975 <- krig1[["1975"]]$predict + alpha * sqrt(krig1[["1975"]]$krige.var)


kkl1979 <- krig0[["1979"]]$predict - alpha * sqrt(krig1[["1979"]]$krige.var)
kku1979 <- krig0[["1979"]]$predict + alpha * sqrt(krig1[["1979"]]$krige.var)


kkl1982 <- krig0[["1982"]]$predict - alpha * sqrt(krig1[["1982"]]$krige.var)
kku1982 <- krig0[["1982"]]$predict + alpha * sqrt(krig1[["1982"]]$krige.var)

results_1975 <- data.frame(
  XUTM = grid[, "XUTM"],
  YUTM = grid[, "YUTM"],
  pred = krig1[["1975"]]$predict,
  var = krig1[["1975"]]$krige.var,
  lower95 = kkl1975,
  upper95 = kku1975
)
coordinates(results_1975) <- ~ XUTM + YUTM


results_1979 <- data.frame(
  XUTM = grid[, "XUTM"],
  YUTM = grid[, "YUTM"],
  pred = krig0[["1979"]]$predict,
  var = krig0[["1979"]]$krige.var,
  lower95 = kkl1979,
  upper95 = kku1979
)
coordinates(results_1979) <- ~ XUTM + YUTM



results_1982 <- data.frame(
  XUTM = grid[, "XUTM"],
  YUTM = grid[, "YUTM"],
  pred = krig0[["1982"]]$predict,
  var = krig0[["1982"]]$krige.var,
  lower95 = kkl1982,
  upper95 = kku1982
)
coordinates(results_1982) <- ~ XUTM + YUTM

library(sp)
# 1975
p1<-spplot(results_1975,pretty=TRUE,scales=list(draw=TRUE), zcol = "pred", main = "Rain data - 1975")
p2<-spplot(results_1975,pretty=TRUE,scales=list(draw=TRUE), zcol = "lower95", main = "Rain data - 1975 (LB)")
p3<-spplot(results_1975,pretty=TRUE,scales=list(draw=TRUE), zcol = "upper95", main = "Rain data - 1975 (UP)")
  gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
# 1979
p1<-spplot(results_1979,pretty=TRUE,scales=list(draw=TRUE), zcol = "pred", main = "Rain data - 1979", col.regions = topo.colors(100))
p2<-spplot(results_1979,pretty=TRUE,scales=list(draw=TRUE), zcol = "lower95", main = "Rain data - 1979 (LB)", col.regions = topo.colors(100))
p3<-spplot(results_1979,pretty=TRUE,scales=list(draw=TRUE), zcol = "upper95", main = "Rain data - 1979 (UP)", col.regions = topo.colors(100))
  gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
# 1982
p1<-spplot(results_1982,pretty=TRUE,scales=list(draw=TRUE), zcol = "pred", main = "Rain data - 1982", col.regions = topo.colors(100))
p2<-spplot(results_1982,pretty=TRUE,scales=list(draw=TRUE), zcol = "lower95", main = "Rain data - 1982 (LB)", col.regions = viridis(100))
p3<-spplot(results_1982,pretty=TRUE,scales=list(draw=TRUE), zcol = "upper95", main = "Rain data - 1982 (UP)", col.regions = viridis(100))
  gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
#1975
mean(results_1975$var, na.rm=TRUE) / var(results_1975$pred, na.rm=TRUE)

results_1975$ci_width_rel <- (results_1975$upper95 - results_1975$lower95) / results_1975$pred
summary(results_1975$ci_width_rel)
#
#library(ggplot2)
#ggplot(results_1975@data, aes(x = pred, y = sqrt(var))) +
#  geom_point(alpha = 0.6) +
#  labs(x = "Predicted value", y = "Prediction Std. Dev.",
#       title = "Uncertainty vs. Predicted value (1975)") +
#  theme_minimal()

#1979
mean(results_1979$var, na.rm=TRUE) / var(results_1979$pred, na.rm=TRUE)

results_1979$ci_width_rel <- (results_1979$upper95 - results_1979$lower95) / results_1979$pred
summary(results_1979$ci_width_rel)
##
##library(ggplot2)
##ggplot(results_1979@data, aes(x = pred, y = sqrt(var))) +
##  geom_point(alpha = 0.6) +
##  labs(x = "Predicted value", y = "Prediction Std. Dev.",
##       title = "Uncertainty vs. Predicted value (1979)") +
##  theme_minimal()

#1982
mean(results_1982$var, na.rm=TRUE) / var(results_1982$pred, na.rm=TRUE)

results_1982$ci_width_rel <- (results_1982$upper95 - results_1982$lower95) / results_1982$pred
summary(results_1982$ci_width_rel)
#
#library(ggplot2)
#ggplot(results_1982@data, aes(x = pred, y = sqrt(var))) +
#  geom_point(alpha = 0.6) +
#  labs(x = "Predicted value", y = "Prediction Std. Dev.",
#       title = "Uncertainty vs. Predicted value (1982)") +
#  theme_minimal()




# xvalid_gaus <- list()


# gaus_error  <- list()
# gaus_cv     <- list()

# for(year in names(years_list)){
# vv <-  variog(years_lIst[[year]], lambda = lambda_mapping[year],trend=~ quota,max.dist = 100)
# vv_fit_gaus <- variofit(vv, cov.model = "gaussian", max.dist = 100)
# geostat_gaus <- likfit(years_list[[year]],
# lambda=lambda_mapping[year],
# cov.model = "gaussian",
# ini.cov.pars=vv_fit_gaus$cov.pars
# )
# xvalid_gaus[[year]] <- xvalid(geodata = years_list[[year]], model = geostat_exp, locations.xvalid = "all")
# gaus_error[[year]] <- mean((xvalid_sph[[year]]$error)^2) # RMSE for Matérn
# gaus_cv[[year]] <-mean(xvalid_sph[[year]]$std.error/years_list[[year]]$data)
# }


