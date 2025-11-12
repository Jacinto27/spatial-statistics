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
# lambda_mapping <- c("1975" = 0, "1979" = 0, "1982" = 0.5)
# Replicating kriging with geoR
lambda_mapping <- c()
for (year in names(years_list)) {
  res <- boxcoxfit(years_list[[year]]$data)
  lambda_mapping[year] <- res$lambda
}

# trend on the dataset
#trend.d<-trend.spatial("1st", geodata=dati_1975, add = ~Quota)
trend.d<-trend.spatial(~quota, geodata=dati_1975)

# trend on the locations where to apply the interpolation
trend.l<-trend.spatial(~Quota, geodata=grid)

# model settings
model.list=model.control(trend.d = trend.d, trend.l = trend.l, cov.model = "spherical",lambda = lambda_mapping['1975'])

# Estimate the likfit for each year to get the aproximate phi value, can also be read by running homework 4 and getting the results

  #vv <- variog(years_list[['1975']], lambda = lambda_mapping[['1975']], trend = ~quota, max.dist = 100)
  #vv_fit_sph <- variofit(vv, cov.model = "spherical", max.dist = 100)
  #geostat_sph <- likfit(years_list[['1975']],
    #lambda = lambda_mapping[['1975']],
    #cov.model = "spherical",
    #ini.cov.pars = vv_fit_sph$cov.pars
  #)

PC75 <- prior.control(phi.prior = "uniform", phi.disc=seq(50,150,1),
tausq.rel.prior= "uniform", tausq.rel.discrete=seq(0,1.0,0.01))

# testing the behavior of the priors
xy75.bayes_try <- krige.bayes(dati_1975, model = model.list, prior = PC75) 

par(mfrow=c(1,2)) 
plot(xy75.bayes_try, type="bars")

# 1000 simulations per location
output=output.control(n.posterior=1000, n.predictive=1000,
quantile=c(0.025,0.975))
PC75<- prior.control(phi.prior = "uniform", phi.disc=seq(68,71,1),
tausq.rel.prior= "uniform", tausq.rel.discrete=seq(0.74,0.84,0.01))
# takes a long time, run only once 
xy75.bayes <- krige.bayes(dati_1975,model = model.list,loc=grid[,1:2],prior = PC75, output = output)
interp75<-interp(grid[,1],grid[,2],xy75.bayes$predictive$mean)

interp75.var<-interp(grid[,1],grid[,2],xy75.bayes$predictive$variance)

par(mfrow=c(1,1)) 
image(interp75)
title("Predicted values for year 1975")

par(mfrow=c(1,1)) 
image(interp75.var)
title("Prediction variance for year 1975")
points(dati_1975,add=T,col="red")

CI75 = xy75.bayes$predictive$quantiles.simulations
head(CI75)

library(scoringutils)
rows <- length(dati_1975$data)
ing75<-matrix(NA,rows,2)
data75 <- dati[dati$anno == 1975, ]
for (i in 1:rows){
LOO_1975_train<-data75[-i,]
LOO_1975_test<-data75[i,]
LOO75_train<-as.geodata(LOO_1975_train,
                        coords.col = c(7,6),data.col= 10, covar.col= 8)
LOO75_test<-as.geodata(LOO_1975_test,
                       coords.col = c(7,6),data.col= 10, covar.col= 8)
trend.d<-trend.spatial(~Quota, geodata=LOO75_train)
trend.l<-trend.spatial(~Quota, geodata=LOO75_test)
model.list=model.control(trend.d = trend.d, trend.l = trend.l,
                         cov.model = "spherical",lambda = lambda_mapping['1975'])
interp<-krige.bayes(LOO75_train, model = model.list,
                    loc=LOO75_test$coords,prior = PC75, output = output)
ing75[i,1:2]<-as.numeric(interp$predictive$quantiles.simulations)
rm(list=c("LOO_1975_train","LOO_1975_test","LOO75_train","LOO75_test","trend.d","trend.l"))
}  
#Function interval_score
interval_score <- function(observed, lower, upper, interval_range, weigh=TRUE, separate_results = FALSE) {
  
  alpha <- (100 - interval_range) / 100
  
  dispersion <- (upper-lower)
  underprediction <- (2/alpha)*(lower-observed)*as.numeric(observed < lower)
  overprediction <- (2/alpha)*(observed-upper)*as.numeric(observed > upper)
  
  score = dispersion + underprediction + overprediction
  
  if (weigh) {
    dispersion <- dispersion * alpha / 2
    underprediction <- underprediction * alpha / 2
    overprediction <- overprediction * alpha / 2
  }
 return(score)
}

IS.75 <-
  interval_score(
  observed = log(data75$totanno), 
  lower = log(ing75[,1]),
  upper = log(ing75[,2]),
  interval_range = 95
)

mean(IS.75)
