# a) Build an R function to choose the parameter p in the inverse distance weighting estimator. Hint: use the idw implemented in gstat.
library(gstat)
library(sf)
 
 p_vals = c(1, 1.5, 2, 2.5, 3, 4, 5) #I assume theres greater value in miniscule differences between 1 and 3 than 4 and 5

# b) Work on the Wolfcamp dataset. Interpolate Wolfcamp data on a 20x20 grid (built from the coordinates) using IDW and multilevel bi-splines (function mba.surf from MBA) and compare the results using a proper score
#  Hint: divide data into training and testing sets

library(geoR)