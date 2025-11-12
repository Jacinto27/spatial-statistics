load("class files/data/RainData.RData")
library(geoR)
library(interp)
library(sf)
colnames(Calabriaquota32) <- colnames(Basquota32)
grid <- rbind(Basquota32[, c(2:4)], Calabriaquota32[, c(2:4)])
grid[, 1:2] <- grid[, 1:2] / 1000
dati[, 6:7] <- dati[, 6:7] / 1000

dati_1975 <- as.geodata(dati[dati$anno == 1975, ], coords.col = 7:6, data.col = 10, covar.col = 8, covar.names = "quota")
dati_1979 <- as.geodata(dati[dati$anno == 1979, ], coords.col = 7:6, data.col = 10, covar.col = 8, covar.names = "quota")
dati_1982 <- as.geodata(dati[dati$anno == 1982, ], coords.col = 7:6, data.col = 10, covar.col = 8, covar.names = "quota")
grid_geodata <- as.geodata(grid, coords.col = 1:2, covar.col = 3)
years_list <- list("1975" = dati_1975, "1979" = dati_1979, "1982" = dati_1982)
# Replicating kriging with geoR

lambda_mapping <- c()
trend.d <- list()
model.list <- list()
variofit <- list()
reltao <- c()

# This trend is only defined once, taken out of the loop for runtime optimization
trend.l <- trend.spatial(~Quota, geodata = grid)

for (year in names(years_list)) {
  # lambda
  res <- boxcoxfit(years_list[[year]]$data)
  lambda_mapping[year] <- res$lambda

  # spatial trend for each year
  trend.d[[year]] <- trend.spatial(~quota, geodata = years_list[[year]])

  # model settings
  model <- ifelse(year %in% c("1979", "1982"), "exponential", "spherical")
  model.list[[year]] <- model.control(trend.d = trend.d[[year]], trend.l = trend.l, cov.model = model, lambda = lambda_mapping[year])

  vv <- variog(years_list[[year]], lambda = lambda_mapping[[year]], trend = ~quota, max.dist = 100)
  variofit[[year]] <- variofit(vv, cov.model = model, max.dist = 100)
  reltao[year] <- variofit[[year]]$nugget / variofit[[year]]$cov.pars[1]
}

# split because the warnings were crowding the terminal
for (year in names(years_list)) {
  cat("Relative Tao=", reltao[year], "Phi=", variofit[[year]]$cov.pars[2], "\n")
}

PC <- prior.control(
  phi.prior = "uniform", phi.disc = seq(10, 90, 1),
  tausq.rel.prior = "uniform", tausq.rel.discrete = seq(0, 1.4, 0.01)
)

## takes a while for results
#xy.bayes_try <- list()
#for (year in names(years_list)) {
#xy.bayes_try[[year]] <- krige.bayes(years_list[[year]], model = model.list[[year]], prior = PC)
#}
#save(xy.bayes_try,file='xy.bayes_try.Rdata')
load('xy.bayes_try.Rdata')

par(mfrow = c(1, 2))
plot(xy.bayes_try[["1975"]], type = "bars", name = "1975")
plot(xy.bayes_try[["1979"]], type = "bars", name = "1979")
plot(xy.bayes_try[["1982"]], type = "bars", name = "1982")

# 1000 simulations per location
output <- output.control(
  n.posterior = 1000, n.predictive = 1000,
  quantile = c(0.025, 0.975)
)

# per-year model control
PC75 <- prior.control(
  phi.prior = "uniform", phi.disc = seq(68, 70, 1),
  tausq.rel.prior = "uniform", tausq.rel.discrete = seq(0.85, 0.91, 0.01)
)

PC79 <- prior.control(
  phi.prior = "uniform", phi.disc = seq(64, 70, 1),
  tausq.rel.prior = "uniform", tausq.rel.discrete = seq(0.25, 0.29, 0.01)
)

PC82 <- prior.control(
  phi.prior = "uniform", phi.disc = seq(58, 65, 1),
  tausq.rel.prior = "uniform", tausq.rel.discrete = seq(0.24, 0.30, 0.01)
)

# takes a long time, run only once

# 1975
#xy75.bayes <- krige.bayes(years_list[[year]], model = model.list[[year]], loc = grid[, 1:2], prior = PC75, output = output)
load("xy75.bayes.Rdata") #20 parameters
interp75 <- interp(grid[, 1], grid[, 2], xy75.bayes$predictive$mean)

interp75.var <- interp(grid[, 1], grid[, 2], xy75.bayes$predictive$variance)

par(mfrow = c(1, 1))
image(interp75)
title("Predicted values for year 1982")

par(mfrow = c(1, 1))
image(interp75.var)
title("Prediction variance for year 1982")
points(dati_1975, add = T, col = "red")

CI75 <- xy75.bayes$predictive$quantiles.simulations
head(CI75)
save(xy75.bayes, file = "xy75.bayes.Rdata")

# 1979
#xy79.bayes <- krige.bayes(years_list[[year]], model = model.list[[year]], loc = grid[, 1:2], prior = PC79, output = output)
load("xy79.bayes.Rdata") #35 paramters
interp79 <- interp(grid[, 1], grid[, 2], xy79.bayes$predictive$mean)

interp79.var <- interp(grid[, 1], grid[, 2], xy79.bayes$predictive$variance)

par(mfrow = c(1, 1))
image(interp79)
title("Predicted values for year 1982")

par(mfrow = c(1, 1))
image(interp79.var)
title("Prediction variance for year 1982")
points(dati_1979, add = T, col = "red")

CI79 <- xy79.bayes$predictive$quantiles.simulations
head(CI79)
save(xy79.bayes, file = "xy79.bayes.Rdata")

# 1982
 xy82.bayes <- krige.bayes(years_list[[year]], model = model.list[[year]], loc = grid[, 1:2], prior = PC82, output = output)
#load("xy82.bayes.Rdata") #56 paramters
interp82 <- interp(grid[, 1], grid[, 2], xy82.bayes$predictive$mean)

interp82.var <- interp(grid[, 1], grid[, 2], xy82.bayes$predictive$variance)

par(mfrow = c(1, 1))
image(interp82)
title("Predicted values for year 1982")

par(mfrow = c(1, 1))
image(interp82.var)
title("Prediction variance for year 1982")
points(dati_1982, add = T, col = "red")

CI82 <- xy82.bayes$predictive$quantiles.simulations
head(CI82)
save(xy82.bayes, file = "xy82.bayes.Rdata")

# CV calculation
library(scoringutils)
rows <- length(dati_1982$data)
ing82 <- matrix(NA, rows, 2)
data82 <- dati[dati$anno == 1982, ]
for (i in 1:rows) {
  LOO_1982_train <- data82[-i, ]
  LOO_1982_test <- data82[i, ]
  LOO82_train <- as.geodata(LOO_1982_train,
    coords.col = c(7, 6), data.col = 10, covar.col = 8
  )
  LOO82_test <- as.geodata(LOO_1982_test,
    coords.col = c(7, 6), data.col = 10, covar.col = 8
  )
  trend.d <- trend.spatial(~Quota, geodata = LOO82_train)
  trend.l <- trend.spatial(~Quota, geodata = LOO82_test)
  model.list <- model.control(
    trend.d = trend.d, trend.l = trend.l,
    cov.model = "exponential", lambda = lambda_mapping["1982"]
  )
  interp <- krige.bayes(LOO82_train,
    model = model.list,
    loc = LOO82_test$coords, prior = PC82, output = output
  )
  ing82[i, 1:2] <- as.numeric(interp$predictive$quantiles.simulations)
  rm(list = c("LOO_1982_train", "LOO_1982_test", "LOO82_train", "LOO82_test", "trend.d", "trend.l"))
}
# Function interval_score
interval_score <- function(observed, lower, upper, interval_range, weigh = TRUE, separate_results = FALSE) {
  alpha <- (100 - interval_range) / 100

  dispersion <- (upper - lower)
  underprediction <- (2 / alpha) * (lower - observed) * as.numeric(observed < lower)
  overprediction <- (2 / alpha) * (observed - upper) * as.numeric(observed > upper)

  score <- dispersion + underprediction + overprediction

  if (weigh) {
    dispersion <- dispersion * alpha / 2
    underprediction <- underprediction * alpha / 2
    overprediction <- overprediction * alpha / 2
  }
  return(score)
}

IS.82 <-
  interval_score(
    observed = log(data82$totanno),
    lower = log(ing82[, 1]),
    upper = log(ing82[, 2]),
    interval_range = 95
  )

mean(IS.82)
