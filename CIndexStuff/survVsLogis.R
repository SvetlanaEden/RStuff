### comparing c-index of suvival vs logistic regression
source("/Users/svetlanaeden/stuff/RStuff/skeMisc.R")


### this function gives a probability of concordance for each pair
### based on their values x1 and x2 (which could be ranks after subjects are ordered)
concProb = function(x1, x2, normalizeBy){
	### probability of concordance increases with distance
	#interpInASetOfPoints(abs(x1 - x2)/normalizeBy, Xs=c(0, (normalizeBy-1)/normalizeBy), Ys=c(.5, .9))
	### probability of concordance increases with distance ever so slightly
	interpInASetOfPoints(abs(x1 - x2)/normalizeBy, Xs=c(0, (normalizeBy-1)/normalizeBy), Ys=c(.5, .6))
	### probability of concordance decreases with distance
	#interpInASetOfPoints(abs(x1 - x2)/normalizeBy, Xs=c(0, (normalizeBy-1)/normalizeBy), Ys=c(.9, .5))
}

### this matrix represents probabilities of concordance b/w pairs of subjects.
### row and column indices are subject ranks in the order of their outcome times

nSubj = 5
probMatr = matrix(0, nrow=nSubj, ncol=nSubj)
concMatr = matrix(0, nrow=nSubj, ncol=nSubj)

for (i in 1:(nrow(probMatr)-1)){
	for (j in (i+1):ncol(probMatr)){
		probMatr[i, j] = concProb(i, j, nSubj)
		concMatr[i, j] = sample(x=c(0, 1), size=1, replace = TRUE, prob = c(1-probMatr[i, j], probMatr[i, j]))
	}
}
#probMatr
#concMatr

survCInd = sum(concMatr)/((nSubj^2 - nSubj)/2)
survCInd

logiCInd = rep(0, nSubj-1)
nPairs = (1:(nSubj-1))*((nSubj-1):1)
for(i in 1:(nSubj-1)){
	logiCInd[i] = sum(concMatr[1:i, (i+1):nSubj])/nPairs[i]
}
logiCInd

plot(1:length(logiCInd), logiCInd, ylim=c(0, 1), type="l")
abline(h=survCInd, col="red")

### Now, let's simmulate concordance with normal distr. 
### First, let's assume that the means of X's are distributed uniformly.
set.seed(20150827)
nSubj = 3
nSim = 500

######### setting the means
######### setting the means
######### setting the means
  #xMeans = runif(nSubj, min=0, max=10)
  xMeans = rnorm(nSubj, mean=10, sd=10)
  #xMeans = c(1,2,3)
  xMeans = xMeans[order(xMeans)]
	
### Now let's simulated each X according to its mean
xValue = matrix(NA, ncol=nSubj, nrow=nSim)
for (i in 1:ncol(xValue)){
	xValue[,i] = rnorm(nSim, mean = xMeans[i], sd=5)
}

# Now let's build a "concordance" matrix for each realization of Xs' (for each row of xValue)
cIndMatr = array(data = 0, dim = c(nSubj, nSubj, nSim))
for (k in 1:nSim){
	for (j in 1:(nSubj-1)){
	  cIndMatr[j, (j+1):nSubj, k] = xValue[k, j] < xValue[k, (j+1):nSubj]
  }
}
# Now let's compute a C-ind for each simulation
cIndForNormSurv = apply(cIndMatr, 3, function(x){sum(x)/((nSubj^2 - nSubj)/2)})

# Now let's compute a C-ind for each simulation for each threashold
nPairs = (1:(nSubj-1))*((nSubj-1):1)
cIndForNormLogi = matrix(NA, nrow=nSim, ncol=nSubj-1)
for (k in 1:nSim){
	### identify order in which a researcher would impose threasholds
	orderedXInd = order(xValue[k,])
	for (j in 1:(nSubj-1)){
		### get rid of indices in the group small than a threashold and larger than a threashold
	  cIndForNormLogi[k, j] = (sum(cIndMatr[,,k]) - sum(cIndMatr[orderedXInd[1:j],orderedXInd[1:j],k]) - sum(cIndMatr[orderedXInd[(j+1):(nSubj)],orderedXInd[(j+1):(nSubj)],k]))/nPairs[j]
  }
	cat("k=", k, "\n")
}
#cIndForNormLogi

plot(cIndForNormSurv, type="n", ylim=c(0, 1), xlim=c(1, nSubj-1))
for (i in 1:nSim){
	lines(cIndForNormLogi[i,], col="gray")
}
abline(h=mean(cIndForNormSurv))
abline(h=(cIndForNormSurv))
lines(apply(cIndForNormLogi, 2, mean), col="red")


############ explore the case with three subjects:
############ explore the case with three subjects:
############ explore the case with three subjects:

set.seed(20150902)
nSubj = 3
nSim = 50000

######### setting the means
  xMeans = (1:nSubj))*3
  xMeans = xMeans[order(xMeans)]
	
### Now let's simulate each X according to its mean
xValue = matrix(NA, ncol=nSubj, nrow=nSim)
for (i in 1:ncol(xValue)){
	xValue[,i] = rnorm(nSim, mean = xMeans[i], sd=5)
}

# Now let's build a "concordance" matrix for each realization of Xs' (for each row of xValue)
cIndMatr = array(data = 0, dim = c(nSubj, nSubj, nSim))
for (k in 1:nSim){
	for (j in 1:(nSubj-1)){
	  cIndMatr[j, (j+1):nSubj, k] = xValue[k, j] < xValue[k, (j+1):nSubj]
  }
}
# Now let's compute a C-ind for each simulation
cIndForNormSurv = apply(cIndMatr, 3, function(x){sum(x)/((nSubj^2 - nSubj)/2)})
statesOfPairs = apply(t(apply(cIndMatr, 3, function(x){as.vector(x)})), 1, paste, collapse = "")
length(table(statesOfPairs))
100*table(statesOfPairs)/nSim

# Now let's compute a C-ind for each simulation for each threashold
nPairs = (1:(nSubj-1))*((nSubj-1):1)
cIndForNormLogi = matrix(NA, nrow=nSim, ncol=nSubj-1)
for (k in 1:nSim){
	### identify order in which a researcher would impose threasholds
	orderedXInd = order(xValue[k,])
	for (j in 1:(nSubj-1)){
		### get rid of indices in the group small than a threashold and larger than a threashold
	  cIndForNormLogi[k, j] = (sum(cIndMatr[,,k]) - sum(cIndMatr[orderedXInd[1:j],orderedXInd[1:j],k]) - sum(cIndMatr[orderedXInd[(j+1):(nSubj)],orderedXInd[(j+1):(nSubj)],k]))/nPairs[j]
  }
	cat("k=", k, "\n")
}
#cIndForNormLogi

plot(cIndForNormSurv, type="n", ylim=c(0, 1), xlim=c(1, nSubj-1))
for (i in 1:nSim){
	lines(cIndForNormLogi[i,], col="gray")
}
abline(h=mean(cIndForNormSurv))
abline(h=(cIndForNormSurv))
lines(apply(cIndForNormLogi, 2, mean), col="red")

