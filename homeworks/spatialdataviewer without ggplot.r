#spatialdataviewer without ggplot

library(geoR)
data(wolfcamp)

# Look at first few rows
head(wolfcamp)

# Structure of the dataset
str(wolfcamp)

# Summary statistics
summary(wolfcamp)
summary(wolfcamp$data)

wolfcamp$x <- wolfcamp$coords[, 1]
wolfcamp$y <- wolfcamp$coords[, 2]
wolfcamp$z <- wolfcamp$data

# Histogram of z values
hist(wolfcamp$z,
     breaks = 30,
     col = "lightblue",
     main = "Distribution of z",
     xlab = "z")

# Simple scatterplot of coordinates

plot(wolfcamp$x, wolfcamp$y,
     pch = 20,
     col = "black",
     xlab = "x",
     ylab = "y",
     main = "wolfcamp coordinates")

# Scatterplot colored by z (using a color ramp)

# Bounding box (extent of coordinates)
range(wolfcamp$x)
range(wolfcamp$y)

#values and stuff
library(lattice)
library(MASS)
dens <- kde2d(wolfcamp$x, wolfcamp$y, n = 100)
filled.contour(dens$x, dens$y, dens$z,
               xlab = "x", ylab = "y",
               main = "Wolfcamp point density (kde2d)")

image(dens$x, dens$y, dens$z)
contour(dens$x, dens$y, dens$z, add = TRUE)
