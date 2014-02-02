
library(PerformanceAnalytics)
data(edhec)
R <- edhec[, 1:4]

##### Starting Values #####
# Startup period length
m <- 10

# Fit a multivariate regression model using the robust affine equivariant 
# method of Rousseeuw et al. (2004)

# This will get me the starting values for the startup period and I can use the
# fitted values from the regression model as the starting value.

# Ignore this step for now and assume the mth period is the vector of starting
# values


##### Smoothed Values #####
smoothedEWMA <- function(Lambda, yt, yhat_lag1){
  p <- length(yt)
  yt <- matrix(yt, nrow=p, ncol=1)
  yhat_lag1 <- matrix(yhat_lag1, nrow=p, ncol=1)
  out <- Lambda %*% yt + (diag(p) - Lambda) %*% yhat_lag1
  return(as.vector(out))
}

# Compute the forecast errors
# \Lambda_{opt} = argmin det \hat{cov}(R)
# \hat{Cov}(R) = 1 \ (T - m) \sum_{t=m+1}^T r_t %*% t(r_t)

# objective function

objEWMA <- function(sol, R, m){
  # Dimension of R
  n <- nrow(R)
  p <- ncol(R)
  
  # Index vector for training period
  idx <- (m+1):n
  
  # yhat for initial value
  yhat <- R[m,]
  
  # Initialize tmp matrix to store the smoothed values
  tmp <- matrix(0, nrow=(n-m), ncol=ncol(R))
  
  # Smoothing matrix
  # Lambda <- matrix(sol, nrow=p, ncol=p)
  Lambda <- diag(sol)
  
  # Initialize counter for tmp
  i <- 1
  for(t in idx){
    tmp[i, ] <- smoothedEWMA(Lambda=Lambda, yt=R[t,], yhat_lag1=yhat)
    yhat <- tmp[i,]
    i <- i + 1
  }
  
  # Compute the forecast error
  tmpR <- coredata(R[(m+1):n, ])
  forecast_error <- tmpR - tmp
  
  # Compute the covariance matrix of the forecast errors
  tmp_sigma <- rep(0, p)
  for(j in 1:nrow(forecast_error)){
    tmp_sigma <- tmp_sigma + (forecast_error[j,] %*% t(forecast_error[j,]))
  }
  tmp_sigma
  eig_vals <- eigen(sm)$values 
  
  if(all(eig_vals >= 0 & eig_vals <= 1)){
    out <- det(tmp_sigma)
  } else {
    out <- 10000
  }
  return(out)
}

# sm <- objEWMA(sol=rep(0.25, 4*4), R, 10)
# sm

# Optimization to find values for the smoothing matrix
# What are the correct upper and lower bounds?
# What is a good starting value?
# Is optim good or should I use DEoptim for optimization
opt <- optim(runif(4), objEWMA, R=R, m=10, method="L-BFGS-B", lower=rep(0, 4^2), upper=rep(1, 4^2))
if(opt$convergence != 0) message("function did not converge")
opt

opt_de <- DEoptim(objEWMA, R=R, m=10, lower=rep(0, 4), upper=rep(1, 4), control=DEoptim.control(itermax=50))
opt_de
