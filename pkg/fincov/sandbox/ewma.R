# http://faculty.chicagobooth.edu/ruey.tsay/teaching/bs41202/sp2011/EWMAvol.R
## Compute exponentially weighted moving average covariance matrix.
## If lambda <= 0, likelihood function is used to estimate lambda.
if(!is.matrix(rtn))rtn=as.matrix(rtn)
Mean=apply(rtn,2,mean)
T=dim(rtn)[1]; k=dim(rtn)[2]
x=rtn
for (i in 1:k){
  x[,i]=rtn[,i]-Mean[i]
}
XX <<- x
#
##
if(lambda > 0){
  h1=1-lambda; Sigt=cov(XX); V1=c(Sigt)
  for (t in 2:T){
    xx=as.numeric(x[t-1,])
    for (i in 1:k){
      Sigt[i,]= h1*xx*xx[i]+lambda*Sigt[i,]
    }
    V1=rbind(V1,c(Sigt))
  }
}
  
x <- matrix(edhec[1, 1:2], nrow=1)
t(x) %*% x

x[1] * x[1]
x[1] * x[2]
x[2] * x[2]


R <- edhec[, 1:4]

# http://faculty.washington.edu/ezivot/econ589/multivariategarch.pdf

# This version stores each predicted covariance matrix in a list 
# Compute the initial covariance matrix
# Should Sigma0 be the sample covariance?
# Or maybe it could be a covariance matrix estimated via ledoit wolf shrinkage?
Sigma0 <- cov(R)
lambda <- 0.94
tmpMat <- matrix(0, nrow=ncol(R), ncol=ncol(R))
listEWMA <- list()
# populate the listEWMA[[1]] with Sigma0
listEWMA[[1]] <- Sigma0
for(i in 2:3){
  tmpMat <- (1 - lambda) * matrix(R[(i-1),], ncol=1) %*% R[(i-1),] + lambda * listEWMA[[i-1]]
  listEWMA[[i]] <- tmpMat
}
listEWMA

# This version continually overwrites tmpMat
Sigma0 <- cov(R)
lambda <- 0.94
tmpMat <- matrix(0, nrow=ncol(R), ncol=ncol(R))
tmpMat <- Sigma0
for(i in 2:3){
  tmpMat <- (1 - lambda) * matrix(R[(i-1),], ncol=1) %*% R[(i-1),] + lambda * tmpMat
}
all.equal(listEWMA[[3]], tmpMat)

# http://r.789695.n4.nabble.com/EWMA-covariance-matrix-td929109.html
set.seed(123)
testData = matrix(rnorm(100), 50, 2)
lam = 0.9
i = 0:49
ewma.wt = lam^i
ewma.wt = ewma.wt/sum(ewma.wt)
cov.ewma = cov.wt(testData, wt=rev(ewma.wt))
cov.ewma

meanx <- colMeans(testData)
x <- testData - matrix(rep(meanx, nrow(testData)), ncol=ncol(testData), byrow=TRUE)
Sigma0 <- cov(x)
lambda <- 0.9
tmpMat <- matrix(0, nrow=ncol(x), ncol=ncol(x))
tmpMat <- Sigma0
for(i in 2:nrow(x)){
  tmpMat <- (1 - lambda) * matrix(x[(i-1),], ncol=1) %*% x[(i-1),] + lambda * tmpMat
}
tmpMat
