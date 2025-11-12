### generate some data ----
require(geoR) 
##### 
set.seed(285) 

# using goR we can simulate GP using grf, by default it implements a Matern covariance with k=0.5 that is an exponential covariance
ex.data <- grf(50, cov.pars=c(10, .25)) 
######## plot results
op <- par(no.readonly = TRUE) #
##some graphical adjustments
 par(mfrow=c(1,2))
 par(mar=c(1,1,2,2)) 
par(mgp = c(2,1,1)) 

######### bubble plot ----
points(ex.data,xlab="",ylab="", main="Bubble graph",sub="diameter proportional \n to the size of the value")
####### variogram ----
plot(variog(ex.data,max.dist=0.4),main="empirical variogram",type="l")
par<-op
#Bayesian kriging ----
## ## use default in the Bayesian kriging 
ex.post <- krige.bayes(ex.data)
## 
names(ex.post)
### visualize the default priors
 ex.post$prior 
 
ex.post$prior$beta # lol
## visualize results
plot(ex.data) 
lines(ex.post, sum = mean,col=2) 
plot(ex.post)

### fix the range 
ex1 <- krige.bayes(ex.data, prior = list(phi.prior = "fixed", phi = 0.3))
## change covariance model
ex1 <- krige.bayes(ex.data, model = list(cov.model="spherical"))
## choose the number of posterior samples
ex1 <- krige.bayes(ex.data, output = list(n.posterior = 100)) 

### some predictions on a 10x10 grid
ex.grid <- as.matrix(expand.grid(seq(0,1,l=10), seq(0,1,l=10)))
## discrete priors for the range phi and the tau_rel parameters
ex.bayes <- krige.bayes(ex.data, loc=ex.grid, 
prior = prior.control(phi.discrete=c(0, 1,2), 
tausq.rel.discrete=c(0,1,2)), output=output.control(n.post=100)) 
#
names(ex.bayes$predictive) 
## Graphical representations ----.. 
## .posterior distributions
par(mfrow=c(1,2))
plot(ex.bayes) #only phi appears as tau_rel is set to zero
plot(ex.data,lty=1) 
lines(ex.bayes, sum = mean, col=2) 
legend(0.2,4,c("estimated on the data", "estimated on the prediction"),lty=1,col=c(1,2),cex=0.8)
## predictive distribution results ----
op <- par(no.readonly = TRUE)
 par(mfrow=c(2,2))
 
image(ex.bayes)
title("predicted values") 
image(ex.bayes, val="variance")
title("prediction variance")
image(ex.bayes, val= "simulation", 
number.col=1)
title("a simulation from the \n predictive distribution") 
image(ex.bayes, val= "simulation", number.col=2)
title("another simulation from \n the predictive distribution") 

par(mfrow=c(1,1))
image(ex.bayes)
title("predicted values") 
points(ex.data,add=T)
image(ex.bayes, val="variance")
title("prediction variance")
points(ex.data,add=T)
## Change the setting: add nugget 
## lpriors for phi and tau_rel are always discrete
PC <- prior.control(phi.prior = c(.1, .2, .3, .2, .1, .1), 
phi.disc=seq(0.1, 0.6, l=6), 
tausq.rel.prior=c(.1, .4, .3, .2), 
tausq.rel.discrete=c(0,.1,.2,.3)) 
ex.user <- krige.bayes(ex.data,loc=ex.grid, prior = PC) 


######### plot results -----
par(mfrow=c(1,2))
plot(ex.user) 
par(mfrow=c(1,1))
plot(ex.data,lty=1) 
lines(ex.user, sum = mean, col=2) 
legend(0.2,4,c("stima sui dati", "stima a posteriori"),lty=1,col=c(1,2),cex=0.8)

 par(mfrow=c(1,2))
 
image(ex.user)
title("predicted values") 
points(ex.data,add=T)
image(ex.user, val="variance")
title("prediction variance")
points(ex.data,add=T)

