library(gstat)
library(sp)
library(raster)
library(gridExtra)
library(gstat)
load("class files/data/RainData.RData")

names(Basquota32)
names(Calabriaquota32)


#names(Calabriaquota32)[names(Calabriaquota32) == "X"] <- "XUTM"
#names(Calabriaquota32)[names(Calabriaquota32) == "Y"] <- "YUTM"
#names(Calabriaquota32)[names(Calabriaquota32) == "quota"] <- "Quota"


#create grid
#grid <- rbind(
#  Basquota32[, c("XUTM","YUTM","Quota")],
#  Calabriaquota32[, c("XUTM","YUTM","Quota")]
#)


colnames(Calabriaquota32)=colnames(Basquota32)
grid=rbind(Basquota32[,c(2:4)],Calabriaquota32[,c(2:4)])

# convert in km
grid$x <- grid$XUTM / 1000
grid$y <- grid$YUTM / 1000

coordinates(grid_prof) <- ~x + y
#idw interpolation

fit_70 <- idw(totanno ~ 1,
              locations = coord_1970,
              newdata = grid,
              idp = idp_70)

fit_77 <- idw(totanno ~ 1,
              locations = coord_1977,
              newdata = grid,
              idp = idp_77)

fit_79 <- idw(totanno ~ 1,
              locations = coord_1979,
              newdata = grid,
              idp = idp_79)

#plot

sp_70 <- spplot(
  fit_70["var1.pred"],
  scales = list(draw = TRUE),
  col.regions = viridis::plasma(6),
  main = "1970",
  xlab = "X UTM (km)",
  ylab = "Y UTM (km)"
)

sp_77 <- spplot(
  fit_77["var1.pred"],
  scales = list(draw = TRUE),
  col.regions = viridis::plasma(6),
  main = "1977",
  xlab = "X UTM (km)",
  ylab = "Y UTM (km)"
)

sp_79 <- spplot(
  fit_79["var1.pred"],
  scales = list(draw = TRUE),
  col.regions = viridis::plasma(6),
  main = "1979",
  xlab = "X UTM (km)",
  ylab = "Y UTM (km)"
)

#every map shows the precipitation interpolated in one year

grid.arrange(sp_70, sp_77, sp_79, ncol = 3)
