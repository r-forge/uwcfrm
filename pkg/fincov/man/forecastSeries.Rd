\name{forecastSeries}
\alias{forecastSeries}
\title{Exponential Smoothing Forecast}
\usage{
  forecastSeries(R, Lambda, startup_period = 10)
}
\arguments{
  \item{R}{xts object of asset returns}

  \item{Lambda}{smoothing matrix}

  \item{startup_period}{periods to use for starting values}
}
\value{
  forecasted values
}
\description{
  Compute the exponential smoothing one-step-ahead forecast
}

