###################################################
### chunk number 1: 
###################################################
rm(list=ls())
source('myfunction.R')
library(spatstat)
library(spdep)
library(scatterplot3d)
library(sgeostat)
library(akima)
load("datiAll.RData")


###################################################
### chunk number 2: 
###################################################
wolf.interp<-interp(wolf$x,wolf$y,wolf$piezo)
persp(wolf.interp$x,wolf.interp$y,wolf.interp$z,xlab="x",ylab="y",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Perspective plot",box=TRUE)


###################################################
### chunk number 3: linear trend
###################################################

aq.fit.1<-lm(piezo~x+y,data=wolf)
xgrid<-seq(min(wolf$x),max(wolf$x),length=20)
ygrid<-seq(min(wolf$y),max(wolf$y),length=20)
aq.grid<- expand.grid(x=xgrid,y=ygrid)
aq.surf<-predict(aq.fit.1,newdata=aq.grid)
persp(xgrid,ygrid,matrix(aq.surf,20,20),xlab="x",ylab="y",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Linear trend",box=TRUE)


###################################################
### chunk number 4: 
###################################################
library(akima)
wolf.interp<-interp(wolf$x,wolf$y,wolf$piezo)
persp(wolf.interp$x,wolf.interp$y,wolf.interp$z,xlab="x",ylab="y",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Perspective plot",box=TRUE)


###################################################
### chunk number 5: Parabolic trend
###################################################

aq.fit.2<-lm(piezo~x+y+I(x^2)+I(y^2)+x*y,data=wolf)
xgrid<-seq(min(wolf$x),max(wolf$x),length=20)
ygrid<-seq(min(wolf$y),max(wolf$y),length=20)
aq.grid<- expand.grid(x=xgrid,y=ygrid)
aq.surf<-predict(aq.fit.2,newdata=aq.grid)
persp(xgrid,ygrid,matrix(aq.surf,20,20),xlab="x",ylab="y",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Parabolic trend",box=TRUE)

###################################################
### chunk number 5b: Cubic trend
###################################################

aq.fit.3<-lm(piezo~x+y+I(x^2)+I(y^2)+x*y+I(x^3)+I(y^3)+I(x^2)*y+I(y^2)*x,data=wolf)
xgrid<-seq(min(wolf$x),max(wolf$x),length=20)
ygrid<-seq(min(wolf$y),max(wolf$y),length=20)
aq.grid<- expand.grid(x=xgrid,y=ygrid)
aq.surf<-predict(aq.fit.3,newdata=aq.grid)
persp(xgrid,ygrid,matrix(aq.surf,20,20),xlab="x",ylab="y",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Cubic trend",box=TRUE)


###################################################
### chunk number 6: 
###################################################
wolf.interp<-interp(wolf$x,wolf$y,wolf$piezo)
persp(wolf.interp$x,wolf.interp$y,wolf.interp$z,xlab="x",ylab="y",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Perspective plot",box=TRUE)



###################################################
### chunk number 7: residuals from trends
###################################################
res.interp<-interp(wolf$x,wolf$y,aq.fit.1$res)
persp(res.interp$x,res.interp$y,res.interp$z,xlab="x",ylab="y",zlab='residuals',theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Residui",box=TRUE)

### there still a lot of variability, but the spatial gradient has been removed
aq.surf<-predict(aq.fit.1,newdata=aq.grid)

par(mfrow=c(1,3))

persp(wolf.interp$x,wolf.interp$y,wolf.interp$z,xlab="x",ylab="y",zlab="",theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Observed",box=TRUE)

persp(xgrid,ygrid,matrix(aq.surf,20,20),xlab="x",ylab="y",zlab='',theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Linear Trend",box=TRUE)

persp(res.interp$x,res.interp$y,res.interp$z,xlab="x",ylab="y",zlab='',theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Residui",box=TRUE)
