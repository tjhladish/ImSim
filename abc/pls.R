#metrics : matrix of metrics computed for each simulation (for ex, 8 metrix and 50,000 simuls will give a 50000x8 matrix)
#vobs: vector of observed metrics (ex, vector of size 8 if you have 8 metrics)

raw_data = read.table("particles.1", sep=',', header=TRUE)
weights = read.table("weights.1", sep=',', header=F)
params = raw_data[,1:4]
metrics = raw_data[,5:12] + 1

vobs = c(0.0483716190476, 0.044, 0.0255449427004, 0.00437161904762, 0.371031746032, 0.232142857143, 0.19246031746, 0.204365079365)
vobs = vobs + 1

lam <- numeric(8)
pdf("/home/tjhladish/Dropbox/Eugene/pls_plots/pls_plots.pdf", width=11, height=8.5)
par(mfrow=c(4,2));for (i in 1:length(lam)) hist(metrics[,i])

library(MASS) ###computes the lambda parameters for the box-cox transformation
#for(i in 1:8){
for(i in (1:length(lam))){
    #sam.i<-1:nrow(metrics)
    print(summary(metrics[,i]))
    b <- boxcox(metrics[,i]~R0+Ih+h+P0, data=params, lambda=seq(-50,50,0.1), plotit=F)
    #b <- boxcox(metrics[,i]~params, lambda=seq(-2,2,0.01), plotit=F)
    lam[i] <- b$x[rev(order(b$y))[1]]
    plot(b, type="l")
}

library(car) ###does the box cox transformation
metsim3 <- as.matrix(metrics)
obsmet <- numeric(8)
for (i in 1:length(lam)) metsim3[,i] <- box.cox(metrics[,i], lam[i])  #box-cox transforms metrics

for (i in 1:length(lam)) obsmet[i] <- box.cox(vobs[i], lam[i])  #box-cox transforms vobs

par(mfrow=c(4,2));for (i in 1:length(lam)) hist(metsim3[,i])


#PLS
library(pls)
samp<-sample(1:nrow(metrics),2000) ##you can't run the PLS with 50000 data points. Too many. So pick 2000 of them.

y <- as.matrix(params[samp,]) ##params is the matrix of parameter values (=particles)
x <- metsim3[samp,]
res <- plsr(y~x, ncomp=8, scale=T,validation="LOO") 
summary(res)
plot(RMSEP(res), legendpos = "topright") #use this to find out how many components you should keep (See PLS reference manual).
                                         #It's a manual step. Can't easily be made automatic. Here, for ex, assume 6 components are sufficient

COMP<- 5
#visually inpsect and pick number of components
plot(res, ncomp = COMP, asp = 1, line = TRUE) 

#Fit the PLS model again with 6 components
res <- plsr(y~x, ncomp=COMP, scale=T,validation="LOO") 

sim.metr <- predict(res, comps=1:COMP, newdata = metsim3, type="scores")  ##8 transformed, orthogonal metrics
obs.metr <- predict(res, comps=1:COMP, newdata = t(obsmet), type="scores")

#euclidean distances, poorly coded but it does the job
d <- numeric(length(sim.metr[,1]))
for (i in 1:COMP){
d <- d+ (sim.metr[,i]-obs.metr[i])^2
}

predictive.prior <- params[order(d)[1:1000],]#simple rejection. Try small numbers
met.reject  <- sim.metr[order(d)[1:1000],]
weights.filtered <- weights[order(d)[1:1000],]

dev.off()
##here you are.
write(t(as.matrix(predictive.prior)), file="predictive_prior.out", append=F, ncol = ncol(predictive.prior))
write(t(as.matrix(met.reject)), file="best_transformed_metrics.out", append=F, ncol = ncol(met.reject))
write(t(as.matrix(weights.filtered)), file="predictive_prior_weights.out", append=F, ncol = 1)
