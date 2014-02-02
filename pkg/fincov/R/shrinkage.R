# I should have a function called covShrink that is a wrapper for lwShrink, bsShrink, etc.
# covShrink

#' Ledoit-Wolf shrinkage covariance estimate
#' 
#' Compute the covariance matrix estimate using the Ledoit-Wolf shrinkage estimate
#' 
#' @param x xts or matrix of asset returns
#' @param shrink shrinkage constant
#' @return covariance matrix estimate
#' @references TODO
#' @note Ported to R from matlab code given at http://www.econ.uzh.ch/faculty/wolf/publications.html#9
#' @author Ross Bennett
#' @export
lwShrink <- function(x, shrink=NULL){
  # port of matlab code from http://www.econ.uzh.ch/faculty/wolf/publications.html#9
  # Ledoit, O. and Wolf, M. (2004).
  # Honey, I shrunk the sample covariance matrix.
  # Journal of Portfolio Management 30, Volume 4, 110-119.
  
  # De-mean returns
  n <- nrow(x)
  p <- ncol(x)
  meanx <- colMeans(x)
  x <- x - matrix(rep(meanx, n), ncol=p, byrow=TRUE)
  
  # Compute sample covariance matrix using the de-meaned returns
  sample <- (1 / n) * (t(x) %*% x)
  
  # Compute prior
  var <- matrix(diag(sample), ncol=1)
  sqrtvar <- sqrt(var)
  tmpMat <- matrix(rep(sqrtvar, p), nrow=p)
  rBar <- (sum(sum(sample / (tmpMat * t(tmpMat)))) - p) / (p * (p - 1))
  prior <- rBar * tmpMat * t(tmpMat)
  diag(prior) <- var
  
  if(is.null(shrink)){
    # What is called pi-hat
    y <- x^2
    phiMat <- t(y) %*% y / n - 2 * (t(x) %*% x) * sample / n + sample^2
    phi <- sum(phiMat)
    
    # What is called rho-hat
    term1 <- (t(x^3) %*% x) / n
    help <- t(x) %*% x / n
    helpDiag <- matrix(diag(help), ncol=1)
    term2 <- matrix(rep(helpDiag, p), ncol=p, byrow=FALSE) * sample
    term3 <- help * matrix(rep(var, p), ncol=p, byrow=FALSE)
    term4 <- matrix(rep(var, p), ncol=p, byrow=FALSE) * sample
    thetaMat <- term1 - term2 - term3 + term4
    diag(thetaMat) <- 0
    rho <- sum(diag(phiMat)) + rBar * sum(sum(((1 / sqrtvar) %*% t(sqrtvar)) * thetaMat))
    
    # What is called gamma-hat
    gamma <- norm(sample - prior, "F")^2
    
    # Compute shrinkage constant
    kappa <- (phi - rho) / gamma
    shrinkage <- max(0, min(1, kappa / n))
  } else {
    shrinkage <- shrink
  }
  # Compute the estimator
  sigma <- shrinkage * prior + (1 - shrinkage) * sample
  out <- list(cov=sigma, prior=prior, shrinkage=shrinkage)
  return(out)
}

# Implements a James-Stein type shrinkage estimate of covariance matrix

#' James-Stein type shrinkage estimate of covariance matrix
#' 
#' Compute the covariance matrix estimate using a James-Stein type shrinkage
#' estimate. This function is implemented using the \code{cov.shrink} function
#' from the corpcor package.
#' 
#' @param x xts or matrix of asset returns
#' @param lambda correlation shrinkage intensity
#' @param lambda.var variance shrinkage intensity
#' @param w optional: weights for each data point
#' @param verbose output status messages while computing
#' @return covariance matrix estimate
#' @references TODO
#' @author Ross Bennett
#' @export
jsShrink <- function(x, lambda, lambda.var, w, verbose=FALSE){
  stopifnot("package:corpcor" %in% search() || require("corpcor", quietly = TRUE))
  
  tmp_out <- corpcor::cov.shrink(x=x, lambda=lambda, lambda.var=lambda.var, w=w, verbose=verbose)
  tmp_attr <- attributes(tmp_out)
  cov <- unclass(tmp_out)
  lambda <- tmp_attr$lambda
  lambda.estimated <- tmp_attr$lambda.estimated
  lambda.var <- tmp_attr$lambda.var
  lambda.var.estimated <- tmp_attr$lambda.var.estimated
  
  # Structure and return
  out <- structure(list(cov=cov, 
                        lambda=lambda,
                        lambda.estimated=lambda.estimated,
                        lambda.var=lambda.var,
                        lambda.var.estimated=lambda.var.estimated))
  return(out)
}


