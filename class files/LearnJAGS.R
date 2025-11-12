#Learning JAGS----
require(R2jags) #this is the library linking R and Jags
require(geoR)
data(wolfcamp)
#the first step will be to prepare a script with the model we want to estimate
## data list ----
## in R we prepare the data we'll pass to JAGS
tobeloaded=list(Z=wolfcamp$data,
                X=wolfcamp$coords[,1],
                Y=wolfcamp$coords[,2],
                D=as.matrix(dist(wolfcamp$coords)),
                N=length(wolfcamp$data),
                az=2,
                bz=4000,
                aw=2,
                bw=4000
)
# parameters in the output----
#we have to define which parameters we want in the output
partosave=c("alfa",
            "betaX",
            "betaY",
            "tauz",
            "tauw",
            "phi")
# define the iterations, burnin and thinning ----
n.iter=10000
n.burn=n.iter/4
n.thin=10
#call jags and run the script (it takes time) ----
model.sim=jags(tobeloaded,
               #inits,
               parameters.to.save=partosave,
               model.file="spatial.txt",
               n.chains=2,
               n.iter=n.iter,
               n.burnin=n.burn,
               n.thin=n.thin
)
#rjags
#Azure
traceplot(model.sim, ask=F,mfrow=c(3,3))
# if are in the same session of the first call to the jags function we can run an update that is we can restart the algorithm at the very last iteration in model.sim and we can obtain longer chains
# first we further check converge
round(model.sim$BUGSoutput$summary,5)
model.sim<-update(model.sim,n.iter=10000,n.burnin=1000,n.thin=10)

# output's post processing ----
require(coda)
names(model.sim$BUGSoutput)
#extract the chains and make them as mcmc obkect in coda
cc<-mcmc(data=model.sim$BUGSoutput$sims.matrix )
plot(cc)
# from this function we can obtain point estimates and credible intervals
summary(cc)
# to obtain HPD that are more reliable then credible intervals:
HPDinterval(cc)
## solution of the exercise----
#phi proportional to 1/range
tobeloaded=list(Z=wolfcamp$data,
                X=wolfcamp$coords[,1],
                Y=wolfcamp$coords[,2],
                invR=solve(exp(-0.0001*as.matrix(dist(wolfcamp$coords)))),
                N=length(wolfcamp$data),
                az=2,
                bz=4000,
                aw=2,
                bw=4000
                
)

partosave=c("alfa",
            "betaX",
            "betaY",
            "tauz",
            "tauw")#,
#            "phi")
# define the iterations, burnin and thinning ----
n.iter=10000
n.burn=n.iter/4
n.thin=10
#call jags and run the script (it takes time) ----
model.sim=jags(tobeloaded,
               #inits,
               parameters.to.save=partosave,
               model.file="spatial_1.txt",
               n.chains=2,
               n.iter=n.iter,
               n.burnin=n.burn,
               n.thin=n.thin
)
traceplot(model.sim, ask=F,mfrow=c(2,3))
# if are in the same session of the first call to the jags function we can run an update that is we can restart the algorithm at the very last iteration in model.sim and we can obtain longer chains
# first we further check converge
round(model.sim$BUGSoutput$summary,5)
model.sim$BUGSoutput$DIC
modelsim0<-update(model.sim,n.iter = 10000,n.thin=n.thin,n.burnin=5000) #phi=0.0001 DIC=1111.619
modelsim1<-model.sim #phi=0.001 DIC=1125.824
modelsim2<-model.sim #phi = 0.01 DIC=1197.067
modelsim3<-update(model.sim,n.iter = 10000,n.thin=n.thin,n.burnin=5000)#phi=0.1 DIC 2224.117
round(modelsim3$BUGSoutput$summary,5)
traceplot(modelsim3, ask=F,mfrow=c(2,3))
modelsim3$BUGSoutput$DIC

###prediction ----
grid<-expand.grid(seq(min(wolfcamp$coords[,1]),max(wolfcamp$coords[,1]),length=10),seq(min(wolfcamp$coords[,2]),max(wolfcamp$coords[,2]),length=10))

tobeloaded=list(Z=wolfcamp$data,
                X=wolfcamp$coords[,1],
                Y=wolfcamp$coords[,2],
                X1=grid[,1],
                Y1=grid[,2],
                invR1=solve(exp(-0.0001*as.matrix(dist(grid)))),
                invR=solve(exp(-0.0001*as.matrix(dist(wolfcamp$coords)))),
                N=length(wolfcamp$data),
                N1=nrow(grid),
                az=2,
                bz=4000,
                aw=2,
                bw=4000
                )

partosave=c("alfa",
            "betaX",
            "betaY",
            "tauz",
            "tauw","Zpred")#,
#            "phi")
# define the iterations, burnin and thinning ----
n.iter=10000
n.burn=n.iter/4
n.thin=10
#call jags and run the script (it takes time) ----
model.sim=jags(tobeloaded,
               #inits,
               parameters.to.save=partosave,
               model.file="spatial_2.txt",
               n.chains=2,
               n.iter=n.iter,
               n.burnin=n.burn,
               n.thin=n.thin
)
pdf("traces.pdf")
traceplot(model.sim, ask=F)
dev.off()
round(model.sim$BUGSoutput$summary,5)
Zpred<-model.sim$BUGSoutput$sims.matrix[,grep("Zpred",colnames(model.sim$BUGSoutput$sims.matrix))]
Zpredhat<-apply(Zpred,2,mean)
ZpredCI<-apply(Zpred,2,quantile,prob=c(0.025,0.975))
