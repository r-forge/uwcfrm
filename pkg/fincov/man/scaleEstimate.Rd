\name{scaleEstimate}
\alias{scaleEstimate}
\title{Robust local estimate of scale}
\usage{
  scaleEstimate(error_matrix, starting_sigma,
    lambda_0 = 0.2)
}
\arguments{
  \item{error_matrix}{matrix of forecast or smoothed values
  erros}

  \item{starting_sigma}{initial condition covariance matrix
  of errors}

  \item{lambda_0}{priori chosen smoothing constant}
}
\value{
  list of local estimates of scale for each time period
  coinciding to \code{error_matrix}
}
\description{
  Robust local estimate of scale
}
\author{
  Ross Bennett
}

