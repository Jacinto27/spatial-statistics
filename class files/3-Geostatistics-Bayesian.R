## Bayesian kriging: libraries -----
library(bmstdr)
# For bmstr check this vignette https://cran.r-project.org/web/packages/bmstdr/vignettes/bmstdr-vig_bookdown.html
# Bayesian kriging requires less data exploration as the variogram model becomes less impacting in this estimation approach
# also read properscoringrules.pdf and Zamo-Naveau2018_Article_EstimationOfTheContinuousRanke.pdf
#from the further reading folder on the elearning to understand the scores provided for model validation
#Bspatial function fits regression models to point referenced spatial data. ----
#It is a wrapper for several available functions. The option package allows to recall spBayes, inla, stan and the bmstdr internal function. The latter is called when package = "none" and model = "spat". It fits a spatial model with no nugget. We can use the function to fit a linear regression model with independent error too. This choice may be very useful for model choice. 
# The option validrows allows to define a validation set 

### Model choice (mchoice = TRUE) ----
load("class files/data/AllGrids.Rdata")
N <- 5000
burn.in<- N/2 #useful for stan and spBayes
n.report<- 2 #for SpBayes
set.seed(12345)#to allow replication of results fix the simulation seed
val<-sort(sample(c(1:nrow(wolf)),trunc(0.1*nrow(wolf)))) #10% of sites used for validation
# independent regression 1 ----
simple <- Bspatial(piezo~x+y,data=wolf,validrows = val, package = "none", model = "lm", N=N,mchoice = T)

# model choice statistics: 
# goodenss of fit
# pdic   pdicalt   dicorig    dicalt    pwaic1    pwaic2     waic1     waic2       gof   penalty      pmcc 
# 4.02      4.07    861.46    861.57      4.42      4.80    861.87    862.62 293424.22 296329.66 589753.89 
#predictive
# rmse     mae    crps     cvg 
# 57.369  53.081  40.077 100.000 

plot(simple)
# Spatial regression with no nugget ----
spat.nonug <- Bspatial(piezo~x+y,data=wolf,validrows = val, package = "none", model = "spat", N=N, coords = 1:2,coordtype = "plain", phi = 0.1, mchoice = T )

# model choice statistics: 
# goodenss of fit
#  pdic   pdicalt   dicorig    dicalt    pwaic1    pwaic2     waic1     waic2       gof   penalty      pmcc 
#4.02      4.07    845.79    845.89      3.66      4.05    832.42    833.19 301497.29 301305.46 602802.76
# predictive
#   rmse     mae    crps     cvg 
#54.389  46.054  24.177 100.000 

plot(spat.nonug)
# Spatial regression with nugget: spBayes ----

spat.spB <- Bspatial(piezo~x+y,data=wolf,validrows = val, package = "spBayes", N=N, coords = 1:2,coordtype = "plain", prior.phi.param = c(0.0005,1), burn.in = N/2, mchoice = T, n.report = n.report )

# model choice statistics: 
# goodenss of fit
#  pdic   pdicalt   dicorig    dicalt    pwaic1    pwaic2     waic1     waic2       gof   penalty      pmcc 
#   4.18      4.65    862.05    862.99      4.47      4.90    862.34    863.21 293575.16 308284.04 601859.20 
# predictive
# rmse     mae    crps     cvg 
#57.536 52.679 30.24    100 

plot(spat.spB)

# is N large enough for convergence of the spatial model?

N1 <- 10000
burn.in1 <- N1/2
spat.spB <- Bspatial(piezo~x+y,data=wolf,validrows = val, package = "spBayes", N=N1, coords = 1:2,coordtype = "plain", prior.phi.param = c(0.0005,1), burn.in = burn.in1, mchoice = T, n.report = n.report )
# model choice statistics: 
# goodenss of fit
# pdic    pdicalt    dicorig     dicalt     pwaic1     pwaic2      waic1      waic2        gof    penalty pmcc
#3.88      4.33    861.43    862.33      4.18      4.53    861.72    862.44 295009.63 309979.83 604989.46  

# predictive
# rmse     mae    crps     cvg 
#57.788 52.91  17.70       100

#some improvement in terms of GoF and much better crps --> improvement in predictive performance

# # Spatial regression with nugget: stan----
spat.stan <- Bspatial(piezo~x+y,data=wolf,validrows = val, package = "stan",  N=N, coords = 1:2,coordtype = "plain", phi = 0.1, burn.in = N/2, mchoice = T, n.report = n.report )

# model choice statistics: 
# goodenss of fit
# pdic   pdicalt   dicorig    dicalt    pwaic1    pwaic2     waic1     waic2       gof   penalty      pmcc 
# 3.97      4.02    861.34    861.45      4.33      4.66    861.70    862.35 293769.38 307512.09 601281.47 
# predictive
# rmse     mae    crps     cvg 
#57.91  53.68  37.36       100

# the only drw