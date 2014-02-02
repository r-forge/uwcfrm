\name{cleanedSeries}
\alias{cleanedSeries}
\title{Exponential Smoothing Robust Cleaned Series}
\usage{
  cleanedSeries(error_matrix, sigma_t, smoothed_values)
}
\arguments{
  \item{error_matrix}{xts object of errors}

  \item{sigma_t}{list of covariance matrix of errors at
  time t}

  \item{smoothed_values}{xts object of smoothed values}
}
\description{
  Compute the robust cleaned series of the exponential
  smoothing model
}
\details{
  \code{sigma_t} must be a list with length equal to the
  number of observations in \code{error_matrix} and
  \code{smoothed_values}.
}

