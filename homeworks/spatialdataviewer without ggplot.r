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

# Histogram of z values
hist(wolfcamp$data,
     breaks = 30,
     col = "lightblue",
     main = "Distribution of z",
     xlab = "z")

# Simple scatterplot of coordinates
wolfcamp$x <- wolfcamp$coords[, 1]
wolfcamp$y <- wolfcamp$coords[, 2]

plot(wolfcamp$x, wolfcamp$y,
     pch = 20,
     col = "black",
     xlab = "x",
     ylab = "y",
     main = "wolfcamp coordinates")

# Scatterplot colored by z (using a color ramp)
cols <- colorRampPalette(c("blue", "green", "yellow", "red"))(100)
z_col <- cols[cut(wolfcamp$z, breaks = 100)]
plot(wolfcamp$x, wolfcamp$y,
     col = z_col,
     pch = 20,
     xlab = "x", ylab = "y",
     main = "wolfcamp points colored by z")
legend("topright", legend = c("low z", "high z"),
       fill = c("blue", "red"), bty = "n")

# Bounding box (extent of coordinates)
range(wolfcamp$x)
range(wolfcamp$y)

#values and stuff
library(lattice)
dens <- kde2d(wolfcamp$x, wolfcamp$y, n = 100)
levelplot(dens$z ~ dens$x * dens$y,
          xlab = "x", ylab = "y",
          main = "wolfcamp density (lattice)")