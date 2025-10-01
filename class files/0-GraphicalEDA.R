###################################################
### chunk number 1:  libraries ----
###################################################
rm(list=ls())
source('myfunction.R') #aggiornare il percorso
library(spatstat)
library(spdep)
library(scatterplot3d)
library(sgeostat)
load("datiAll.RData")



###################################################
### chunk number 2:  histogram ----
###################################################
summary(wolf$piezo)
par(mfrow=c(1,1),pty='s',cex.lab=0.7,cex.main=1.0)
hist(wolf$piezo,  freq=FALSE, breaks=20,
     main = "Piezometric head measurements \n at the wolf Aquifer (Texas, USA)",xlab='Piezometric head')
lines(density(wolf$piezo))


###################################################
### chunk number 3: scatter plot of the coordinates ----
###################################################
par(mfrow=c(1,1),pty='s')
plot(wolf$x,wolf$y,pch=20,xlab="x",ylab="y")
title(main="wolfcamp data") 


###################################################
### chunk number 4:  3D scatterplot ----
###################################################
par(mfrow=c(1,1),pty='s')
scatterplot3d(wolf$x,wolf$y,wolf$piezo,
              type = "h", highlight.3d = TRUE,pch = 20, 
              main = "Piezometric head measurements \n at the wolf Aquifer (Texas, USA)",xlab='x',ylab='y',zlab='Piezometric head' )


###################################################
### chunk number 5: 3D surface ----
###################################################
par(mfrow=c(1,1),pty='s')
library(akima)
wolf.interp<-interp(wolf$x,wolf$y,wolf$piezo)
persp(wolf.interp$x,wolf.interp$y,wolf.interp$z,xlab="x",ylab="y",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Perspective plot",box=TRUE)

#or

require(fields)
surface(wolf.interp,xlab="x",ylab="y",zlab='z',theta=60,
        phi=30,expand=0.9,ltheta = 120, shade = 0.75, type="p",
        main="Observed values",box=TRUE)

###################################################
### chunk number 6: simbol map four classes ----
###################################################
  require(ggplot2)
  wolf1<-wolf
  wolf1$val<-cut(wolf[,3], breaks=quantile(wolf[,3]),include.lowest = T)
    ggplot(wolf1)+geom_point(aes(x=x,y=y,fill=val,shape=val,col=val),show.legend = T,size=2.2)+theme_light()+xlab("x")+ylab("y")+
      scale_color_manual("Piezoelectric height", values = c("blue", "green", "orange","darkred"))+
      scale_fill_manual("Piezoelectric height", values = c("blue", "green", "orange","darkred"))+
      scale_shape_manual("Piezoelectric height", values = c(15, 20, 3,17))
      


###################################################
### chunk number 9: Contour plots ----
###################################################
surface(wolf.interp)


#########################################
######### Example using meuse data gstat ----
#############################################

par(mfrow=c(1,1))
require(gstat)
require(sp)
    data(meuse)
s3d <- scatterplot3d(meuse[,c(1,2,6)], type = "h", highlight.3d = TRUE, pch = 20, main = "Zinc concentration")


#We now build a square grid using the meuse coordinates and verify if the horizontal and vertical directions have different means
xgrid<-seq(min(meuse$x),max(meuse$x),length=29)
ygrid<-seq(min(meuse$y),max(meuse$y),length=41)
meuse.grid<- expand.grid(x=xgrid,y=ygrid)
fx<-factor(meuse$x,levels=xgrid)
fy<-factor(meuse$y,levels=ygrid)
meuse.aov<-aov(meuse$zinc~meuse$x+meuse$y)
summary(meuse.aov)  

#only the vertical direction seems to have a significant cintribution, we may expect some anisotropy

par(mfrow=c(1,1))
plot(meuse.grid,type="n",main="Meuse data",frame.plot=F)
for (i in 1:dim(meuse.grid)[1]){
  abline(h=meuse.grid$y[i], col="lightgray")
  abline(v=meuse.grid$x[i],col="lightgray")
}
points(meuse$x,meuse$y,col=2,pch=18)

### Notice data samples are not aligned to the regular grid we built, if we want to analyze them useing a regular grid we need to align them to the grid itself
# a possible choice is to assign the samples to the nearest grid point
# for example we may compute the absolute difference between observed coordinates and grid coordinates 
########################################################################

mmx<-rep(0,155)
mmy<-rep(0,155)
for (i in 1: 155){
  a<-xgrid-meuse$x[i]
  b<-ygrid-meuse$y[i]
  mmx[i]<-match(min(abs(a)),abs(a))
  mmy[i]<-match(min(abs(b)),abs(b))}
cc.meuse<-matrix(0,155,2)
for (i in 1:155){
  cc.meuse[i,]<-c(xgrid[mmx[i]],ygrid[mmy[i]])
}
points(cc.meuse,col="blue",pch=20)

####repeat the anova 
fx<-as.factor(cc.meuse[,1])
fy<-as.factor(cc.meuse[,2])
################ ANOVA meuse ##########
meuse.aov<-aov(meuse$zinc~fx+fy)
summary(meuse.aov)
#no changes


###################################################
### chunk number 10:  exploring spatial variation ----
###################################################
par(mfrow=c(1,1),pty='s')
wolf.sp<-wolf
coordinates(wolf.sp)<-~x+y
wolf.cloud<-variogram(piezo~1,cloud = T,data = wolf.sp)
plot(wolf.cloud,main='variogram cloud',pch=20,cex=0.8)


###################################################
### chunk number 11:  sample variogram ----
###################################################
wolf.bin<-variogram(piezo~1,cloud = F,data = wolf.sp)
ggplot(wolf.cloud)+geom_point(aes(x=dist,y=gamma),col="blue", alpha=0.3)+
  geom_point(data=wolf.bin, aes(x=dist,y=gamma),col="red",size=3)+
  theme_light()+xlab("distance")+ylab("semivariogram")+ggtitle("empirical variogram")
###################################################
### chunk number 12:  variogram in 4 directions ----
###################################################
wolf.aniso <- variogram(piezo ~ 1,wolf.sp, alpha = c(0, 45, 90, 135))
ggplot(wolf.aniso)+geom_line(aes(x = dist, y= gamma, group=factor(dir.hor),col=factor(dir.hor)))+
  theme_light()+xlab("distance")+ylab("semivariogram")+
  scale_color_manual("directions",values = c("black","green","blue","red"))+
  ggtitle("variogram in 4 directions")

