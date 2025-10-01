# notes from class

install fields for output

usefl function ?expand.grid

idw location, formula argument for the power P (idp)

## OCTOBER 1

To implemente ipd interpolator the distance msut be omuted based on the ype of distinace
- long lat: greater circle distance:    checking in R if this distance is implemented
- if not, must onvert into (metric) coordinates to then use euclidean distance

best beggingi model is a lienar model:
z(s)=\mu(s)+error

however note the assumptions:
- dnormal distrubution of error
- mean 0
- variance constant
- independent errors betwwen 

errors are not independent they might be linked across space, and the variance might not be constant,
instead of assuming independence we can assume dependence,

if we create a full gaussean covaraince matrix we can identify the dependence in the gaussian setting.

since were doing a gaussiean normal with covariancce

z~N(\mu,\sigma (covariance))

we could probalby reduce the number of elements by grouping the elements together, perhaps by distances.
first hypotehsis we can make to simplify the covariance matrix

our covariance is the relationship of two distances, because of tobblers law we know closer objects are more positively correlated, therefore our covariane matrix can be populate dby vaues of a function we can define such as:
h() where h is a distance function between two points, and c is our function. we want to define c such as it has a higher value when the distance is short, and decreases exponentially? (why?)

covariance(s_1,s_2)= C(h())

sine stationarity also helps us in reduce values of covariance we can also assume (or induce) stationrity in our data to have a simpler covariance matrix

the best stationary model to apply is second order stationarity.

this stationary assumption in stationary in space is as such:
we have a process X() with a fixed value u, our distribution (underlying process) doesnt change if you move your observation by u units to any direction

non stationarity i thhe mean is not important becaude we can simply check the properties for z(s)-B_0+bx(s))=E(s)

we can resutrctrute it to have a process with mean 0 this is less of a problem than non stationary covariance

returnng to the covarinnce function

small note sometimes if we move accorss the same route on oppiste direction the values mayb be different:
west to east :c1
east to west: c2

think of lakes with currents, moving in one direction will require less distance than the other
this property is called aniisotropic property, we want our process ot be isotropic, so oour distance is jsut hte model of our distance.

so we want:

weaakly stationary:
    - constant mean (this is not a serious problem because i can fit a spatial trend)
    - covariance that depends on the distance
    - and that it is (isotropic)
in the slides X(s) is the same as Z(s)
theres a process for brownian motion that allows us to assert this (for stationary and isotropic models):
\gamna(u)=E(X(s+u)-X(s))^2=Var(X(s+u)=X(s))=C(0)-C(u) where c is our covariance function
