\name{smoothedSeries}
\alias{smoothedSeries}
\title{Exponential Smoothing Smoothed Values}
\usage{
  smoothedSeries(R, Lambda, starting_values)
}
\arguments{
  \item{R}{xts object of asset returns}

  \item{Lambda}{smoothing matrix}

  \item{starting_values}{xts object of}
}
\value{
  xts object of smoothed series in the same time index of R
}
\description{
  Compute the exponential smoothing smoothed values. The
  smoothed values are computed for the same timeframe as R.
  The starting values should be immediately before R.
}
\author{
  Ross Bennett
}

