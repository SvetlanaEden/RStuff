### exploring c-index for a case with 
### three subjects with exponential survival time with age x1 > x2 > x3
### and beta. We assume that the x's are spaced evenly: x1 - x2 = x2 - x3 = delta

abline

beta  = 1
delta = seq(.01, 10, .1)
x3 = 0

### a case when compute a c-index of the survival model:
### a case when compute a c-index of the survival model:
### a case when compute a c-index of the survival model:
numer = exp(3*beta*x3 + 2*beta*delta)*(3*exp(3*beta*delta) + 2 + 2*exp(beta*delta) + 2*exp(2*beta*delta))
denom = exp(3*beta*x3)*(exp(2*beta*delta)+exp(beta*delta)) * (exp(beta*delta) + 1) * (exp(2*beta*delta) + 1)
cIndexSurv = (1/3)*numer/denom
plot(delta, cIndexSurv)


### a case when compute a c-index of the logistic regression model:
### a case when compute a c-index of the logistic regression model:
### a case when compute a c-index of the logistic regression model:
numerLogist = exp(3*beta*x3+2*beta*delta)*(2*exp(3*beta*delta) + 3*exp(2*beta*delta) + 1)
cIndexLogist = (1/2)*numerLogist/denom
plot(delta, cIndexLogist)

par(mar=c(7,4,3,2))
plot(delta, cIndexLogist, type="l", col="red", ylim=c(0, 1))
lines(delta, cIndexSurv)
lines(delta, cIndexLogist - cIndexSurv, col="green")
