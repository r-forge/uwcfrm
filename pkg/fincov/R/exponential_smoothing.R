
#' Exponential Smoothing Constructor
#' 
#' If the smoothing matrix is NULL, it will be solved for through optimization.
#' 
#' @param R xts object of asset returns
#' @param smoothing_matrix smoothing matrix
#' @param startup_period number of periods to use for startup data to estimate starting values
#' @param training_period number of periods to use for training data to estimating smoothing matrix
#' @param type
#' @author Ross Bennett
#' @export
ExponentialSmoothing <- function(R, smoothing_matrix=NULL, startup_period=10, training_period=36, type=c("classic", "robust")){
  
  # Match the argument for type
  type <- match.arg(type)
  
  # Get dimensions of R
  n <- nrow(R)
  p <- ncol(R)
  
  # Checks
  if(startup_period <= p) stop("startup_period must be greater than number of assets")
  if(startup_period > n) stop("startup_period must be less than number of observations")
  if(training_period > n) stop("training_period must be less than number of observations")
  if(training_period < startup_period) stop("training_period must be greater than startup_period")
  
  # Initialize list to return
  out <- list()
  call <- match.call()
  out$call <- call
  
  # List of initial model specifications
  init <- list()
  init$R <- R
  init$startup_period <- startup_period
  init$training_period <- training_period
  
  # List of initial data
  data <- list()
  data$startup_data <- R[1:startup_period,]
  # should the training data include startup data?
  data$training_data <- R[(startup_period+1):training_period,]
  data$test_data <- R[(training_period+1):n,]
  
  out$data <- data
  
  # List of starting values
  starting_values <- list()
  # TODO: Fit a model to compute the starting values
  # Use the startup data to fit a model for the starting values
  # The Croux paper uses a robust mcd regression to compute starting values
  # If the returns are the response variable, what should be used as the
  # predictor variable(s)?
  # Other methods?
  # GARCH, copula, ...?
  # Use the actual data for now
  starting_values$fitted <- data$startup_data
  starting_values$sigma <- cov(data$startup_data)
  
  out$starting_values <- starting_values
  
  if(!is.null(smoothing_matrix)){
    # Checks for smoothing matrix
    # p x p
    if(!all.equal(dim(smoothing_matrix), c(p, p))){
      stop("smoothing matrix must have dimensions p x p")
    }
    
    # eigen values
    eig_vals <- eigen(smoothing_matrix)$values
    if(any(eig_vals < 0) | any(eig_vals > 1)){
      warning("eigenvalues are not in [0, 1]")
    }
    Lambda <- smoothing_matrix
  } else {
    if(type == "robust"){
      robust <- TRUE
    } else {
      robust <- FALSE
    }
    # Optimization to compute optimal smoothing matrix
    Lambda <- solveLambda(data$training_data, starting_values$fitted, robust)
  }
  
  init$smoothing_matrix <- Lambda
  
  out$init <- init
  
  tmp_data <- R[(startup_period+1):n,]
  
  # Calculate the smoothed series
  smoothed_data <- smoothedSeries(tmp_data, Lambda, starting_values$fitted)
  out$smoothed_data <- smoothed_data
  
  # Calculate the forecast series
  # forecast errors
  
  # Calculate the error
  # smoothed errors
  error_smoothed <- tmp_data - smoothed_data
  out$error <- error_smoothed
  
  if(robust){
    robust_estimates <- list()
    # Calculate the scale estimate
    scale_estimate <- scaleEstimate(error_smoothed, starting_values$sigma, lambda_0=0.2)
    out$robust_scale_estimate <- scale_estimate
    
    # Calculate the cleaned series
    cleaned_series <- cleanedSeries(error_smoothed, scale_estimate, smoothed_data)
    out$robust_cleaned_series <- cleaned_series
    
    # calculate the robust smoothed series
    robust_smoothed <- smoothedValuesRobust(Lambda, cleaned_series, smoothed_data, starting_values$fitted)
    out$robust_smoothed_series <- robust_smoothed
  }
  # Structure and return
  return(structure(out, class=paste("exponential_smoothing", type, sep="_")))
}

# \hat{y}_t = \Lambda y_t + (I - \Lambda) \hat{y}_{t-1})
# Compute the smoothed values of a single time period (i.e. \hat{y}_{t|t-1})
.smoothed_values <- function(Lambda, y, y_lag1){
  # Lambda should be a p x p matrix
  if(length(y) != length(y_lag1)) stop("y and y_lag1 must be the same length")
  p <- length(y)
  y <- matrix(as.numeric(y), ncol=1)
  y_lag1 <- matrix(as.numeric(y_lag1), ncol=1)
  out <- Lambda %*% y + (diag(p) - Lambda) %*% y_lag1
  return(as.numeric(out))
}

#' Exponential Smoothing Smoothed Values
#' 
#' Compute the exponential smoothing smoothed values. The smoothed values are
#' computed for the same timeframe as R. The starting values should be immediately
#' before R.
#' 
#' @param R xts object of asset returns
#' @param Lambda smoothing matrix
#' @param starting_values xts object of 
#' @return xts object of smoothed series in the same time index of R
#' @author Ross Bennett
#' @export
smoothedSeries <- function(R, Lambda, starting_values){
  
  # Check R
  # Check starting_values
  
  # Get dimensions of asset returns
  n <- nrow(R)
  p <- ncol(R)
  
  # Checks for smoothing matrix
  # p x p
  if(!all.equal(dim(Lambda), c(p, p))){
    stop("smoothing matrix must have dimensions p x p")
  }
  
  # eigenvalues
  eig_vals <- eigen(Lambda)$values
  if(any(eig_vals < 0) | any(eig_vals > 1)){
    warning("eigenvalues are not in [0, 1]")
  }
  
  # Get the starting values
  y_m <- as.numeric(last(starting_values))
  y_lag1 <- y_m
  
  # Use the rest of the data to calculate the smoothed values
  smoothed_mat <- matrix(0, nrow=nrow(R), ncol=ncol(R))
  colnames(smoothed_mat) <- colnames(R)
  for(i in 1:nrow(smoothed_mat)){
    tmp_y <- R[i,]
    smoothed_mat[i,] <- .smoothed_values(Lambda, tmp_y, y_lag1)
    y_lag1 <- smoothed_mat[i,]
  }
  smoothed_mat <- xts(smoothed_mat, index(R))
  return(smoothed_mat)
}

# \hat{y}_{T+1|T} = \hat{y} = \Lambda \sum_{k=0}^{T-1} (I - \Lambda)^k y_{T-k}
# Forecast the next value y_{T+1} at time T
# Single period one-step-ahead forecast values
.forecast_values <- function(R, Lambda){
  # Get dimensions of asset returns
  n <- nrow(R)
  p <- ncol(R)
  
  # Checks for smoothing matrix
  # p x p
  if(!all.equal(dim(Lambda), c(p, p))){
    stop("smoothing matrix must have dimensions p x p")
  }
  
  # eigenvalues
  eig_vals <- eigen(Lambda)$values
  if(any(eig_vals < 0) | any(eig_vals > 1)){
    warning("eigenvalues are not in [0, 1]")
  }
  
  # p x p identity matrix
  identity <- diag(p)
  
  tmp <- 0
  for(k in 0:(n-1)){
    tmp_y <- matrix(as.numeric(R[(n-k),]), ncol=1)
    tmp <- tmp + (identity - Lambda)^k %*% tmp_y
  }
  out <- Lambda %*% tmp
  return(as.numeric(out))
}

#' Exponential Smoothing Forecast
#' 
#' Compute the exponential smoothing one-step-ahead forecast
#' 
#' @param R xts object of asset returns
#' @param Lambda smoothing matrix
#' @param startup_period periods to use for starting values
#' @return forecasted values
#' @export
forecastSeries <- function(R, Lambda, startup_period=10){
  # Get dimensions of asset returns
  n <- nrow(R)
  p <- ncol(R)
  
  # Checks for smoothing matrix
  # p x p
  if(!all.equal(dim(Lambda), c(p, p))){
    stop("smoothing matrix must have dimensions p x p")
  }
  
  # eigenvalues
  eig_vals <- eigen(Lambda)$values
  if(any(eig_vals < 0) | any(eig_vals > 1)){
    warning("eigenvalues are not in [0, 1]")
  }
  
  # Check startup_period
  if(startup_period <= p) stop("startup_period must be greater than number of assets")
  if(startup_period > n) stop("startup_period must be less than number of observations")
  
  # Total number of forecast periods
  m <- n - startup_period
  
  # Get the date index of the asset returns
  idx <- index(R)
  
  # Compute the one-step-ahead forecast
  tmp_forecast <- matrix(0, nrow=m, ncol=p)
  colnames(tmp_forecast) <- colnames(R)
  for(i in 1:nrow(tmp_forecast)){
    tmp_forecast[i,] <- .forecast_values(R[1:(startup_period+i-1),], Lambda)
  }
  forecast <- xts(tmp_forecast, idx[(startup_period+1):n])
  # Alternate way
  #tmp_forecast <- matrix(0, nrow=n, ncol=p)
  #colnames(tmp_forecast) <- colnames(R)
  #for(i in 1:m){
  #  tmp_forecast[(i+startup_period),] <- .forecast_values(R[1:(startup_period+i-1),], Lambda)
  #}
  #tmp_forecast[1:startup_period,] <- NA
  #forecast <- xts(tmp_forecast, idx)
  return(forecast)
}

#' Objective function to find optimal smoothing matrix using classic methods
#' 
#' Objective function used in \code{\link{solveLambda}} for computing the 
#' optimal smoothing matrix via optimization. The objective value to be 
#' minimized is the determinant of the covariance matrix of the errors.
#' 
#' The time index of R, smoothed values, and error matrix are all equal.
#' 
#' @param params vector of parameters for smoothing matrix
#' @param R xts object of asset returns. This should be the training data to estimate the smoothing matrix.
#' @param starting_values starting values for computing the smoothed values.
#' @return determinant of covariance of error matrix
#' @author Ross Bennett
#' @export
objLambdaClassic <- function(params, R, starting_values){
  
  # Dimensions of asset returns
  n <- nrow(R)
  p <- ncol(R)
  
  # Construct the Lambda matrix
  Lambda <- matrix(params, nrow=p, ncol=p)
  
  # All eigen values of smoothing matrix must be in [0, 1]
  eig_vals <- eigen(Lambda)$values
  if(any(eig_vals < 0) | any(eig_vals > 1)){
    out <- 10000
  } else {
    
    # Compute the smoothed values given the training data
    train_smoothed <- smoothedSeries(R, Lambda, starting_values)
    
    # Compute the forecast values given the training data
    # The forecast is very slow
    #train_forecast <- forecastValues(training_data, Lambda, startup_period)
    
    # Compute the error matrix
    error_mat <- coredata(R) - coredata(train_smoothed)
    #error_mat <- coredata(tmp_data) - coredata(train_forecast)
    
    # Compute the sample covariance of the errors
    tmp_cov <- 0
    for(i in 1:nrow(error_mat)){
      tmp_e <- error_mat[i,]
      tmp_cov <- tmp_cov + tmp_e %*% t(tmp_e)
    }
    cov_error <- (1 / nrow(error_mat)) * tmp_cov
    
    # Compute the determinant of the covariance of the error matrix. This is the
    # objective value in the function to minimize.
    out <- det(cov_error)
    
    # The determinant is very small, will I have tolerance issues with the 
    # optimization?
    
    # sum of squared errors
    # out <- sum(apply(error_mat, 1, function(x) sum(x^2)))
  }
  return(out)
}

#' Objective function to find optimal smoothing matrix using robust methods
#' 
#' Objective function used in \code{\link{solveLambda}} for computing the 
#' optimal smoothing matrix via optimization. The objective value to be 
#' minimized is the determinant of the covariance matrix of the errors.
#' 
#' The time index of R, smoothed values, and error matrix are all equal.
#' 
#' @param params vector of parameters for smoothing matrix
#' @param R xts object of asset returns. This should be the training data to 
#' estimate the smoothing matrix.
#' @param starting_values starting values for computing the smoothed values.
#' @return determinant of covariance of error matrix
#' @author Ross Bennett
#' @export
objLambdaRobust <- function(params, R, starting_values){
  
  # Dimensions of asset returns
  n <- nrow(R)
  p <- ncol(R)
  
  # Construct the Lambda matrix
  Lambda <- matrix(params, nrow=p, ncol=p)
  
  # All eigen values of smoothing matrix must be in [0, 1]
  eig_vals <- eigen(Lambda)$values
  if(any(eig_vals < 0) | any(eig_vals > 1)){
    out <- 10000
  } else {
    
    # Compute the smoothed values given the training data
    train_smoothed <- smoothedSeries(R, Lambda, starting_values)
    
    # Compute the forecast values given the training data
    # The forecast is very slow
    #train_forecast <- forecastValues(training_data, Lambda, startup_period)
    
    # Compute the error matrix
    error_mat <- coredata(R) - coredata(train_smoothed)
    #error_mat <- coredata(tmp_data) - coredata(train_forecast)
    
    #h <- floor((nrow(error_mat) + ncol(error_mat) + 1) / 2)
    h <- floor(0.75 * nrow(error_mat))
    cov_error <- cov.rob(error_mat, quantile.used=h, method="mcd", nsamp="sample")$cov
    
    # Compute the determinant of the covariance of the error matrix. This is the
    # objective value in the function to minimize.
    out <- det(cov_error)
    
    # The determinant is very small, will I have tolerance issues with the 
    # optimization?
    
    # sum of squared errors
    # out <- sum(apply(error_mat, 1, function(x) sum(x^2)))
  }
  return(out)
}

#' Solve for the optimal smoothing matrix
#' 
#' Solve for the optimal smoothing matrix by minimizing the determinant of the
#' covariance of the error matrix
#' 
#' @param R xts object of asset returns. This should be the training data used
#' to solve for the optimal smoothing matrix.
#' @param starting_values starting values used as initial condition to estimate 
#' smoothed values.
#' @param robust TRUE/FALSE (default=TRUE). If TRUE, use robust methods to solve
#' for Lambda, the optimal smoothing matrix.
#' @param method optimization method
#' @return optimal smoothing matrix
#' @export
solveLambda <- function(R, starting_values, robust=TRUE, method=c("DEoptim")){
  require(DEoptim)
  
  # Match the argument for the optimization method
  method <- match.arg(method)
  
  # Dimensions of returns
  n <- nrow(R)
  p <- ncol(R)
  
  # DEoptim method
  if(method == "DEoptim"){
    # Is lower value of 0 and upper value of 1 the correct bounds?
    if(robust){
      opt <- DEoptim(objLambdaRobust, 
                     R=R, 
                     starting_values=starting_values,
                     lower=rep(0, p^2), 
                     upper=rep(1, p^2),
                     control=list(itermax=50))
    } else {
    opt <- DEoptim(objLambdaClassic, 
                   R=R, 
                   starting_values=starting_values, 
                   lower=rep(0, p^2), 
                   upper=rep(1, p^2),
                   control=list(itermax=50))
    }
    opt_par <- opt$optim$bestmem
  }
  # Construct the Lambda matrix
  Lambda <- matrix(opt_par, nrow=p, ncol=p)
  return(Lambda)
}

#' Robust local estimate of scale
#' 
#' @param error_matrix matrix of forecast or smoothed values erros
#' @param starting_sigma initial condition covariance matrix of errors
#' @param lambda_0 priori chosen smoothing constant
#' @return list of local estimates of scale for each time period coinciding to \code{error_matrix}
#' @author Ross Bennett
#' @export
scaleEstimate <- function(error_matrix, starting_sigma, lambda_0=0.2){
  
  # \Sigma_{t-1}
  sigma_lag1 <- starting_sigma
  
  # Compute \Sigma_t for each time step
  sigma_list <- vector("list", nrow(error_matrix))
  for(i in 1:length(sigma_list)){
    tmp_et <- as.numeric(error_matrix[i,])
    sigma_list[[i]] <- .scale_estimate(tmp_et, sigma_lag1, lambda_0=lambda_0)
    sigma_lag1 <- sigma_list[[i]]
  }
  names(sigma_list) <- index(error_matrix)
  return(sigma_list)
}

# helper function used in scaleEstimate
.scale_estimate <- function(e_t, sigma_lag1, lambda_0=0.2){
  tuning_constant <- qchisq(0.95, length(e_t))
  e_t <- matrix(e_t, ncol=1)
  tmp <- as.numeric(t(e_t) %*% solve(sigma_lag1) %*% e_t)
  out <- lambda_0 * (.biweight(sqrt(tmp), tuning_constant) / tmp) * e_t %*% t(e_t) + (1 - lambda_0) * sigma_lag1
  return(out)
}

# biweight function used in .scale_estimate
.biweight <- function(x, c){
  # the constant gamma_cp is selected such that E[\rho_{c,p}(||X||)] = p
  # where X is a p-variate normal distribution
  # ||X|| is the norm of the element x of a normed vector space (i.e. sqrt(x %*% x))
  # How to calculate gamma_cp?
  gamma_cp <- 1
  if(abs(x) <= c){
    out <- gamma_cp * (1 - (1 - (x / c)^2)^3)
  } else {
    out <- gamma_cp
  }
  return(out)
}

# huber function used in .cleaned_series
.huber <- function(x, k){
  min(k, max(x, -k))
}

# Single period cleaned series y^*_t 
.cleaned_series <- function(sigma_t, e_t, smoothed_y_t){
  boundary_value <- qchisq(0.95, length(e_t))
  e_t <- matrix(e_t, ncol=1)
  sigma_inv <- solve(sigma_t)
  tmp <- sqrt(as.numeric(t(e_t) %*% sigma_inv %*% e_t))
  out <- (.huber(tmp, boundary_value) / tmp) * e_t + smoothed_y_t
  return(out)
}

#' Exponential Smoothing Robust Cleaned Series
#' 
#' Compute the robust cleaned series of the exponential smoothing model
#' 
#' \code{sigma_t} must be a list with length equal to the number of observations
#' in \code{error_matrix} and \code{smoothed_values}.
#' 
#' @param error_matrix xts object of errors
#' @param sigma_t list of covariance matrix of errors at time t
#' @param smoothed_values xts object of smoothed values
cleanedSeries <- function(error_matrix, sigma_t, smoothed_values){
  if(nrow(error_matrix) != nrow(smoothed_values)) stop("error_matrix and smoothed_values must have the same number of rows")
  if(!is.list(sigma_t)) stop("sigma_t must be a list")
  if(nrow(error_matrix) != length(sigma_t)) stop("The length of sigma_t must equal the number of observations in error_matrix and smoothed_values")
  
  if(!all.equal(as.character(index(error_matrix)), names(sigma_t), check.attributes=FALSE)){
    stop("index of error_matrix must match names of sigma_t")
  }
  
  cleaned_series <- matrix(0, nrow=nrow(error_matrix), ncol=ncol(error_matrix))
  colnames(cleaned_series) <- colnames(R)
  for(i in 1:nrow(cleaned_series)){
    tmp_et <- as.numeric(error_matrix[i,])
    tmp_smoothed <- as.numeric(smoothed_values[i,])
    tmp_sigma <- sigma_t[[i]]
    cleaned_series[i,] <- .cleaned_series(tmp_sigma, tmp_et, tmp_smoothed)
  }
  cleaned_series <- xts(cleaned_series, index(error_matrix))
  return(cleaned_series)
}

#' Smoothed Values 
#' 
#' Exponential smoothing smoothed values using robust estimates. The time index 
#' for the computed robust smoothed values will be the same as the index of
#' \code{cleaned_series} and \code{smoothed_values}.
#' 
#' @param Lambda
#' @param cleaned_series xts of cleaned series
#' @param smoothed_values xts of smoothed values
#' @param starting_values starting values to use for initial condition when computing robust smoothed values
#' @return xts object of robust smoothed values
#' @author Ross Bennett
#' @export
smoothedValuesRobust <- function(Lambda, cleaned_series, smoothed_values, starting_values){
  
  if(!all.equal(index(cleaned_series), index(smoothed_values))){
    stop("index of cleaned_series must be equal to index of smoothed_values")
  }
  
  if(!all.equal(dim(cleaned_series), dim(smoothed_values))){
    stop("dimension of cleaned_series must be dimension to index of smoothed_values")
  }
  
  y_m <- as.numeric(last(starting_values))
  
  smoothed_robust <- matrix(0, nrow=nrow(cleaned_series), ncol=ncol(cleaned_series))
  colnames(smoothed_robust) <- colnames(cleaned_series)
  for(i in 1:nrow(cleaned_series)){
    if(i ==1){
      y_lag1 <- y_m
    } else {
      y_lag1 <- as.numeric(smoothed_values[(i-1),])
    }
    tmp_cleaned <- as.numeric(cleaned_series[i,])
    smoothed_robust[i,] <- .smoothed_values(Lambda, tmp_cleaned, y_lag1)
  }
  smoothed_robust <- xts(smoothed_robust, index(cleaned_series))
  return(smoothed_robust)
}

