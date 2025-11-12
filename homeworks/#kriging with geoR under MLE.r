#kriging with geoR under MLE
library(geoR)
data("wolfcamp")
#2 step procedure---
#Step1: choose variogram model---
vv<-variog(wolfcamp,trend="1st",max.dist = 250)
plot(vv,type="b")
vv.exp<-likfit(wolfcamp,trend = "1st",cov.model="exponential", ini.cov.pars = c(2000,100))
lines(vv.exp,col=2)
summary(vv.exp)
##Step2: run kriging with fixed covariance ----
grid<- expand.grid(seq(min(wolfcamp$coords[,1]),max(wolfcamp$coords[,1]),length=20),
                   seq(min(wolfcamp$coords[,2]),max(wolfcamp$coords[,2]),length=20)          )

kk.c<-krige.control(trend.d = "1st",trend.l = "1st",obj.model = vv.exp)

krig1<-krige.conv(wolfcamp,locations = grid,krige = kk.c)

image(krig1)
contour(krig1,add=F,nlev=20)

names(krig1)
#confidence intervals 95%----
low<-krig1$predict-1.96*sqrt(krig1$krige.var)
up<-krig1$predict+1.96*sqrt(krig1$krige.var)
