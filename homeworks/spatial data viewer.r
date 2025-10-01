#spatial data viewer
# Packages
library(ggplot2)
library(sf)

# 1) Load the data
#data("x", package = "x")
data("Wolfcamp", package = "geoR")
# Take a look at the first rows
head(Wolfcamp)

# Check dimensions and structure
str(Wolfcamp)
summary(Wolfcamp)

# 2) Basic summaries of the z variable (the thing you’ll interpolate)
hist(Wolfcamp$z, breaks = 30, col = "skyblue", main = "Distribution of z",
     xlab = "z value")

# 3) Scatterplot of coordinates, colored by z
ggplot(Wolfcamp, aes(x = x, y = y, color = z)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  coord_equal() +
  theme_minimal() +
  labs(title = "Wolfcamp data", color = "z")

# 4) Quick sf conversion (if you want to treat it spatially)
wc_sf <- st_as_sf(Wolfcamp, coords = c("x", "y"), crs = NA)

# Bounding box
st_bbox(wc_sf)

# 5) Optional: check spatial density
ggplot(Wolfcamp, aes(x, y)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  geom_point(alpha = 0.4) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Sampling density of Wolfcamp points")
