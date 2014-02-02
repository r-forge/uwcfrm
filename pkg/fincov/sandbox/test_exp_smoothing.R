# Classic Exponential Smoothing

library(PerformanceAnalytics)
data(edhec)
R <- edhec[, 1:2]

# Estimating the smoothing matrix

# Step 1: select a training period
# break the data into a training set and a test set

# must satisfy training_period <= nrow(R)

# Step 2: select a startup period
# This is how we obtain the starting values
# The paper uses a robust multivariate linear regression model to compute 
# starting values. Work on this part later

# m = startup period
# must satisfy m > p

# calculate:
# y_m (fitted values from model)
# Sigma_m (covariance matrix of model residuals)

# Need to calculate the one-step-ahead forecast errors
# r_t = y_t - \hat{y}_{t | t-1}
# is y_t the smoothed value or the actual value? I think it is the actual value
# is \hat{y}_{t | t-1} the smoothed value or the forecast value?

# It makes sense to use the forecast (or predicted) values, but the notation 
# makes me think it is the smoothed values

# R := {r_{m+1}, ..., r_T}
# r_{m+1} indicates that I do not use the starting values in the forecast.
# If m=10, then the first forecasted value is at t=11 (i.e. R[11,])

#####


# DEoptim
set.seed(123)
opt1 <- DEoptim(objLambda, R=R, startup_period=10, training_period=36,
                lower=rep(0, p^2), upper=rep(1, p^2), control=DEoptim.control(itermax=50))
opt1$optim$bestmem
opt1$optim$bestval

# optim
# not sure if optim is a good method for this problem

# init_vals <- as.numeric(diag(p)) * 0.8
# lb <- rep(0, p^2)
# ub <- rep(1, p^2)
# opt2 <- optim(init_vals, objLambda, 
#               R=R, startup_period=10, training_period=36, 
#               method="L-BFGS-B", lower=lb, upper=ub,
#               control=list(factr=.Machine$double.eps))
# opt2

set.seed(123)
tmpL <- solveLambda(R=R, startup_period=10, training_period=36)

# The paper uses a grid search for optimization
# I could also logic from random portfolios for a 'random' optimization method
#opt3 <- gridSearch(objLambda, R=R, startup_period=10, training_period=36, lower=lb, upper=ub)
#opt3$minfun
#opt3$minlevels

tmp <- classicExponentialSmoothing(R=R, startup_period=12, training_period=48)

head(tmp$smoothed_values)
head(R[49:152])

plot(tmp$smoothed_values[,1], type="l")
lines(R[49:152,1], col="red")

plot(tmp$smoothed_values[,2], type="l")
lines(R[49:152,2], col="red")

