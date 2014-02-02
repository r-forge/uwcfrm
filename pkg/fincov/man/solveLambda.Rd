\name{solveLambda}
\alias{solveLambda}
\title{Solve for the optimal smoothing matrix}
\usage{
  solveLambda(R, starting_values, robust = TRUE,
    method = c("DEoptim"))
}
\arguments{
  \item{R}{xts object of asset returns. This should be the
  training data used to solve for the optimal smoothing
  matrix.}

  \item{starting_values}{starting values used as initial
  condition to estimate smoothed values.}

  \item{robust}{TRUE/FALSE (default=TRUE). If TRUE, use
  robust methods to solve for Lambda, the optimal smoothing
  matrix.}

  \item{method}{optimization method}
}
\value{
  optimal smoothing matrix
}
\description{
  Solve for the optimal smoothing matrix by minimizing the
  determinant of the covariance of the error matrix
}

