\name{smoothedValuesRobust}
\alias{smoothedValuesRobust}
\title{Smoothed Values}
\usage{
  smoothedValuesRobust(Lambda, cleaned_series,
    smoothed_values, starting_values)
}
\arguments{
  \item{Lambda}{}

  \item{cleaned_series}{xts of cleaned series}

  \item{smoothed_values}{xts of smoothed values}

  \item{starting_values}{starting values to use for initial
  condition when computing robust smoothed values}
}
\value{
  xts object of robust smoothed values
}
\description{
  Exponential smoothing smoothed values using robust
  estimates. The time index for the computed robust
  smoothed values will be the same as the index of
  \code{cleaned_series} and \code{smoothed_values}.
}
\author{
  Ross Bennett
}

