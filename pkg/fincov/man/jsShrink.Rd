\name{jsShrink}
\alias{jsShrink}
\title{James-Stein type shrinkage estimate of covariance matrix}
\usage{
  jsShrink(x, lambda, lambda.var, w, verbose = FALSE)
}
\arguments{
  \item{x}{xts or matrix of asset returns}

  \item{lambda}{correlation shrinkage intensity}

  \item{lambda.var}{variance shrinkage intensity}

  \item{w}{optional: weights for each data point}

  \item{verbose}{output status messages while computing}
}
\value{
  covariance matrix estimate
}
\description{
  Compute the covariance matrix estimate using a
  James-Stein type shrinkage estimate. This function is
  implemented using the \code{cov.shrink} function from the
  corpcor package.
}
\author{
  Ross Bennett
}
\references{
  TODO
}

