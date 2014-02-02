\name{objLambdaRobust}
\alias{objLambdaRobust}
\title{Objective function to find optimal smoothing matrix using robust methods}
\usage{
  objLambdaRobust(params, R, starting_values)
}
\arguments{
  \item{params}{vector of parameters for smoothing matrix}

  \item{R}{xts object of asset returns. This should be the
  training data to estimate the smoothing matrix.}

  \item{starting_values}{starting values for computing the
  smoothed values.}
}
\value{
  determinant of covariance of error matrix
}
\description{
  Objective function used in \code{\link{solveLambda}} for
  computing the optimal smoothing matrix via optimization.
  The objective value to be minimized is the determinant of
  the covariance matrix of the errors.
}
\details{
  The time index of R, smoothed values, and error matrix
  are all equal.
}
\author{
  Ross Bennett
}

