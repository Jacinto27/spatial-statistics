# Notes on the first homework  

## a) Build an R function to choose the parameter p in the inverse distance weighting estimator. Hint: use the idw implemented in gstat.

For this first part we must find the optimal value for p for the IDW (inverse distance weighted) interpolation model. The IDW model (as the name suggests) is an inverse function with a p parameter that serves as it's power.The weight has this form:

$\displaystyle\sum_{k=1}^n \frac{1}{d^p}

Where d is a distance value. We must guess our values for p until we find one that minimizes our meassurement error. We do this by dividing our data into training and testing datapoints, training our data with a certain p value, estimating the missing values in our test data and meassuring their estimation error. The value that has the smallest error is our "optimal" p value. Our professor susggested we first check the values between 1 and 10 on intervals of 1, and then "zoom in" on the values that minimize our meassurement error and test out those ranges on intervals of 0.2. 

(Because our true optimal value for the wolfcamp dataset is 3.8, the range that we tested at first was ok since it went from 0 to 5)

There are two ways to split our data to find our optimal values:
- K-fold cros validation: fast but imprecise 
- Leave-one-out cross validation: thorough but slow

K-fold cross validation involves splitting the data into k groups for training and testing and doing our estimations on all of the groups one by one. For example, if we choose k=5 we split the data into 5 groups, we train our data in 4 of these groups and test it on the 5th, then we repeat the process but replacing our training and testing groups so that we train on a different set of 4 groups and a different 5th test group, we repeat the process until all of our groups have been a tested at least once. Our groups can be split by randomly choosing points, or by choosing cutoff points in our data (as in, ordering our data by longitude and dividing it by 5, for example).

Leave-one-out cross validation is a much thorough validation process where we only leave out 1 data point out of our dataset, creat our model just to estimate that one data point, and then reinsert the data point back and choose another data point to leave out, we repeat this until we have tested our entire dataset This is the algorithm that we're gonna be using on the entire homework. Leave-one-out cross validation is essentially a k-fold validation where k=1. Notice how in loocv we don't need to choose our training and testing datasets as our whole data is both training and testing at the same time.

For our homework, we simply make a loocv algorithm that finds the best p value for the wolfcamp dataset.  

## b) Work on the Wolfcamp dataset. Interpolate Wolfcamp data on a 20x20 grid (built from the coordinates) using IDW and multilevel bi-splines (function mba.surf from MBA) and compare the results using a proper score

For this part we must interpolate our data into a grid. In our case we want a very rough grid of 20x20 squares, for this we must essentially choose a model and predict every grid (so, 400 estimations). Luckily the packages already have functions that do this for us, we must only divide our territory into the grids and then interpolate them using those packages.

The professor wants us to test two interpolation methods to construct this grid:
- idw (what we did on the first part)
- mba.surf (surface estimation)

IDW is essentially taking random point and comparing it to a sample value (rain data, height, biomass, etc) and then (inversly) weighting that value by the distance our point has to the sample, essentially saying that the closer we are to a value the more "similar" our reading is. We do this for all points in our field until we have a surface.

Bi-spline interpolation is the opposite effect, it's finding the smoothest path betwene the sample points using 2dimensional splines. We can imagine this as each sample point being a stick that grows upwards depending on it's value, and we just "lay a blanket over them", the resulting shape that the blanket takes is the result of our bi-spline interpolation. Bi-spline interpolation takes a few more parameters, (namely number of intervals on the x axis (n), number of intervals on the y axis (m), number of levels (h) ). To find the optimal value for our bi-spline, we perform a leave one out CV algorithm like we did for the IDW, but this time we test a few random combinations between n, m and h.

After we do that, we simply choose which interpolation model between IDW or mba.surf has the best results by picking the one with the lowest meassurement error. Since our algorithm meassures the errors using the Rooted Mean Square Error (RMSE) function this serves as the "proper score" the professor talks about.