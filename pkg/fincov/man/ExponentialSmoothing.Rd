\name{ExponentialSmoothing}
\alias{ExponentialSmoothing}
\title{Exponential Smoothing Constructor}
\usage{
  ExponentialSmoothing(R, smoothing_matrix = NULL,
    startup_period = 10, training_period = 36,
    type = c("classic", "robust"))
}
\arguments{
  \item{R}{xts object of asset returns}

  \item{smoothing_matrix}{smoothing matrix}

  \item{startup_period}{number of periods to use for
  startup data to estimate starting values}

  \item{training_period}{number of periods to use for
  training data to estimating smoothing matrix}

  \item{type}{}
}
\description{
  If the smoothing matrix is NULL, it will be solved for
  through optimization.
}
\author{
  Ross Bennett
}

