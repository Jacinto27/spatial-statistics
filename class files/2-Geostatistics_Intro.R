####################################################################################
######################## Variograms and Covariances -----

require(MASS)
require(MBA)
require(fields)
require(sf)
require(gstat)
require(geoR)
set.seed(1234)
# we build a function that simulates a Gaussian spatial process with exponential coveriance ----
#First we randomly generated coordinates in the interval of the plan [0,1]x[0,1]
dat<-data.frame(x=runif(100),y=runif(100))
grf.exp<-function(coord,m=NULL,sigma=1,rr=0.25,nugget=0){
  require(MASS)
  require(sp)
  require(sf)
  require(fields)
  #mean =0 by default
  #sigma=1 and range=0.25, nugget=0
  if(is.null(m)){m<-rep(0,nrow(coord))}
  D<-as.matrix(dist(coord))
  SigmaI<-diag(nugget,nrow=nrow(coord))+sigma*Exponential(D,range=rr)
  value<-(mvrnorm(n=1,mu=m,Sigma = SigmaI))
  dat<-data.frame(coord,value)
  
  dat<-SpatialPointsDataFrame(data=dat,coords=dat[,1:2])
  #coordinates(dat)<-~x+y
  return(dat)
}
## generate an example ----
sim1 <- grf.exp(dat)
## plot values and locations

bubble(sim1)
spplot(sim1)
   
# empirical and theoretical variograms ----
require(gstat)
sim1.v<-variogram(value~1,data=sim1) #empirical variogram ----
plot(sim1.v,sub=expression(paste(sigma^2,"=1,",phi,"=0.25")),main="exponential",ylab="semivariogram",type="b")
## theoretical variogram ----
m <- vgm(1, "Exp",range= 0.25, nugget=0.12)
plot(sim1.v,model=m)
#variogram in 4 directions ----
plot(variogram(value~1,sim1, alpha=c(0,45,135, 270)), type="b")
#on the same plot 
plot(variogram(value~1,sim1, alpha=c(0,45,135,270)), multipanel=F)
#we can use this plot to verify if the directional variograms are really different

#### surface representation
sim1.sur<-mba.surf(cbind(coordinates(sim1),sim1$value),no.X=100,no.Y=100)$xyz.est ### we only need the estimated surface
fields::surface(sim1.sur,main=expression(paste("exponential ",sigma^2,"=1,",phi,"=0.25")))

## Exercise 1 ----
# build a function that simulates from a spherical Gaussian process 

## We build a function that simulates from a matern Gaussian process ----

dat<-data.frame(x=runif(100),y=runif(100))
grf.mat<-function(coord,m=NULL,sigma=1,rr=0.25,nugget=0,kappa=0.5){
  require(MASS)
  require(sp)
  require(sf)
  require(fields)
  #mean =0 by default
  #sigma=1, range=0.25, nugget=0, smoothness kappa=0.5
  if(is.null(m)){m<-rep(0,nrow(coord))}
  D<-as.matrix(dist(coord))
    SigmaI<-as.matrix(diag(nugget,nrow=nrow(coord))+sigma*Matern(D, range=rr,nu=kappa))
  value<-(mvrnorm(n=1,mu=m,Sigma = SigmaI))
  dat<-data.frame(coord,value)
  dat<-SpatialPointsDataFrame(data=dat,coords=dat[,1:2])
  #coordinates(dat)<-~x+y
  return(dat)
}

sim2<-grf.mat(dat,nugget=0.1, kappa=2.5)

bubble(sim2)
spplot(sim2)
sim2.surf<-mba.surf(cbind(coordinates(sim2),sim2$value),no.X=100,no.Y=100)$xyz.est ### we only need the estimated surface
surface(sim2.surf,main=expression(paste("matern ",sigma^2,"=1,",phi,"=0.25,",kappa,"=2.5")))
plot(variogram(value~1,data=sim2),sub=expression(paste(sigma^2,"=1,",phi,"=0.25,",kappa,"=2.5" )),main="matern",ylab="semivariogram",type="b")

#let's see a different smoothing

sim3<- grf.mat(data.frame(x=runif(100),y=runif(100)),kappa=0.1, nugget = 0.2,sigma=0.5)
bubble(sim3)
sim3.surf<-mba.surf(cbind(coordinates(sim3),sim3$value),no.X=100,no.Y=100)$xyz.est ### we only need the estimated surface
#surface and bvariogram
surface(sim3.surf,main=expression(paste("matern ",sigma^2,"=0.5,",phi,"=0.25,",kappa,"=0.1,"," nugget=0.2")))
plot(variogram(value~1,data=sim3),sub=expression(paste(sigma^2,"=0.5,",phi,"=0.25,",kappa,"=0.1",  ", nugget=0.2" )),main="matern",ylab="semivariogram",type="b")


par(mfrow=c(1,1))

####### more simulations to define the variogram ----

sim5<-grf.mat(data.frame(x=runif(100),y=runif(100)),sigma=1.2,rr=0.5,nugget=0.5,kappa=0.1)
sim6<-grf.exp(data.frame(x=runif(100),y=runif(100)),rr=0.5,nugget=0.5,sigma=1.2)
sim5.v<-variogram(value~1,sim5)
sim6.v<-variogram(value~1,sim6)
#pdf("../introductionAndGeostat/Images/sim5sim6Variog.pdf")
par(mfrow=c(1,2))
plot(sim5.v$dist,sim5.v$gamma,ylim=c(0,max(sim5.v$gamma+.05)), xlim=c(0,max(sim5.v$dist)),ylab=expression(gamma(h)),type="b",main="matern")
arrows(0,0,0,0.5,code=3)
text(0.05,0.25,"nugget")
plot(sim6.v$dist,sim6.v$gamma,ylim=c(0,max(sim6.v$gamma+.05)), xlim=c(0,max(sim6.v$dist)),ylab=expression(gamma(h)),type="b",main="exponential")
arrows(0,0,0,0.5,code=3)
text(0.01,0.25,"nugget")
#dev.off()
par(mfrow=c(1,1))

#plot theoretical variogram----
show.vgms(models = c("Exp", "Mat", "Gau"), nugget = 0.1,)

# show a set of Matern models with different smoothness:
show.vgms(kappa.range = c(.1, .2, .5, 1, 2, 5, 10), max = 10,as.groups = T)

######################### Anisotropy ----
# see the notes here https://r-spatial.github.io/gstat/reference/vgm.html
require(sp)
data(meuse)
meuse.sp <- meuse
coordinates(meuse.sp) <- ~x + y
vgm.aniso <- variogram(log(zinc) ~ 1, meuse.sp, alpha = c(0, 45, 90, 135))
plot(vgm.aniso,type="b")
plot(vgm.aniso,multipanel=FALSE)

# Redo with ggplot
require(ggplot2)
zinc.aniso.plot <- ggplot(aes(x = dist, y = gamma), data = vgm.aniso)  #Initialize
zinc.aniso.plot <- zinc.aniso.plot + geom_point()  # Tell it to plot points
zinc.aniso.plot <- zinc.aniso.plot + facet_wrap(~dir.hor)  # tell it to separate by dir.hor
zinc.aniso.plot
# clean up the axes a bit.
zinc.aniso.plot <- zinc.aniso.plot + scale_x_continuous(limits = c(0, 1500), 
                                                        expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1.2), expand = c(0, 
                                                                                                                              0))
zinc.aniso.plot
### Spatial correlation seems to be strongest in the 45 degree direction, and weakest in the 135 degress direction (things stay similar for longer if you move in the 45 degree direction). 

# We fit the model in the 45 degree direction. ----
vgm.45 <- subset(vgm.aniso, vgm.aniso$dir.hor == 45)  # Subset the variogram data for just this direction
plot(vgm.45,type="b")
#from the plot we can try a graphical estimation of the model picking the parameters from the plot
eye.sph <- vgm(0.42, "Sph", 1100, 0.05)
plot(vgm.45, eye.sph)
#We might try an automatic regression procedure to fit it. This might be more objective.

fit.sph <- fit.variogram(vgm.45, vgm(0.42, "Sph", 1200, 0.05))
fit.sph
plot(vgm.45, fit.sph)

#gaussian covariance sigma^2 exp(-h^2/range^2)

fit.exp <- fit.variogram(vgm.45, vgm(0.42, "Exp", 1100/3, 0.05))
fit.mat.025 <- fit.variogram(vgm.45, vgm(0.42, "Mat", 1100, 0.05, kappa = 0.25))
fit.mat.100 <- fit.variogram(vgm.45, vgm(0.42, "Mat", 1100/4, 0.05, kappa = 1))
fit.gau <- fit.variogram(vgm.45, vgm(0.42, "Gau", 1100/sqrt(3), 0.05))

plot(vgm.45, fit.exp)
plot(vgm.45,fit.mat.025)
### Plot the other variograms and choose one 
#Kriging with trend  - Regression Kriging ----

#We will now fit a regression kriging model to the meuse data. The dist from the river is appropriate as a covariate. Also, we have this as a raster layer for the entire study region.

#First, it is good to verify that distance is an appropriate predictor of log zinc.

#We'll also check for anisotropy. The zinc is anisotropic, but maybe that's completely explained by elevation differences. We should check.

summary(lm(log(zinc) ~ dist, data = meuse))
# Then we should plot the variogram of the residuals, and see if there is
# spatial structure.
vgm <- variogram(log(zinc) ~ dist, locations = ~x + y, data = meuse)
plot(vgm)
# Plot anisotropic variogram
vgm.aniso <- variogram(log(zinc) ~ dist, meuse.sp, alpha = c(0, 45, 90, 135))
plot(vgm.aniso)

##Any difference could just be noise. When in doubt, keep it simple.
fit.reg.sph <- fit.variogram(vgm, vgm(0.25, "Sph", 1500, 0.02))
# note, it's important that the regressor - dist - is present in both the
# sample data (meuse) and in the prediction data (meuse.grid).
data (meuse.grid)
pred.reg <- krige(log(zinc) ~ dist, ~x + y, meuse, meuse.grid, model = fit.reg.sph)
## plot kriging ----
pred.plot <- ggplot(aes(x = x, y = y), data = pred.reg)
pred.plot <- pred.plot + geom_tile(aes(fill = var1.pred))
pred.plot <- pred.plot + scale_fill_gradient("log(zinc)",low = "yellow", high = "blue")
pred.plot + coord_equal()+theme_light()

#The Smoothing Effect of Kriging: ----
#The kriging map is the best prediction (when the data are Gaussian! otherwise, it's just a good map, but not the best). But don't think the kriging map is a fully realistic interpretation of the true, underlying process. In particular, kriging is a smoother, it loses texture from the original process. Remember that the kriging prediction is the mean of the underlying Gaussian process. 

#the variogram of the fitted data is not close to the observed one
pred.reg$dist <- meuse.grid$dist
plot(variogram(var1.pred ~ dist, ~x + y, data = pred.reg), fit.reg.sph, ylim = c(0,0.3))
# a possible correction is obtained using conditional simulations ----
# we simulate from a Gaussian process with exactly the same covariance structure as the one observed
sim.reg <- krige(log(zinc) ~ dist, ~x + y, meuse, meuse.grid, model = fit.reg.sph,nsim = 9, nmax = 20)

#to draw a nice plot Convert from 'wide' format to 'long' format
library(reshape)
sim.gdf <- melt(sim.reg, id.vars = c("x", "y"), variable_name = "sim")
names(sim.gdf)

sim.plot <- ggplot(aes(x = x, y = y), data = sim.gdf) + facet_wrap(~sim)
sim.plot <- sim.plot + geom_tile(aes(fill = value))
sim.plot <- sim.plot + scale_fill_gradient(low = "blue", high = "yellow")
sim.plot + coord_equal()+theme_light()

sim.reg$dist <- meuse.grid$dist
plot(variogram(sim1 ~ dist, ~x + y, data = sim.reg), fit.reg.sph, ylim = c(0, 
                                                                           0.3))


require(geoR)
### GeoR exploratory analysis----
data("wolfcamp")
class(wolfcamp)
summary(wolfcamp)
plot(wolfcamp)
vv<-variog(wolfcamp,trend="1st", max.dist = 250)
plot(vv,type="b")
eyefit(vv) #graphical estimation

#least square estimation
vv.exp<-variofit(vv,fix.nugget=F,nugget=1000)
plot(vv,type="l")
lines(vv.exp,col=2)
## likelihood estimation
vv.exp.l<-likfit(wolfcamp,trend="1st",cov.model = "exp",ini.cov.pars = c(2000,200))
lines(vv.exp.l,col=4)
#cross validation
vv.l.v<-xvalid(wolfcamp,model=vv.exp.l)
names(vv.l.v)
### crossvalidation coefficient
mean(vv.l.v$std.error^2)

vv.wls.v<-xvalid(wolfcamp,model=vv.exp)
mean(vv.wls.v$std.error^2)

## geoR spatial interpolation-----

x<-seq(min(wolfcamp$coord[,1]), max(wolfcamp$coord[,1]))
y<-seq(min(wolfcamp$coord[,2]), max(wolfcamp$coord[,2]))
grid<-expand.grid(x,y)
kk<-krige.control(trend.d="1st",trend.l="1st", obj.model=vv.exp.l,cov.model="exp", cov.pars=c(3715.87359,34.12164))
krig.1<-krige.conv(wolfcamp,locations = grid,krige=kk)
names(krig.1)
