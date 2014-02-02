\name{covEstimate}
\alias{covEstimate}
\title{Estimate covariance matrix}
\usage{
  covEstimate(R,
    method = c("classic", "mcd", "m_estimate", "mm_estimate", "mve", "ogk", "sde", "s_estimate", "lw_shrinkage", "js_shrinkage"),
    control = list())
}
\arguments{
  \item{R}{xts or matrix of asset returns}

  \item{method}{The method used to compute the covariance
  estimate}

  \item{control}{named list of arguments}
}
\description{
  Estimate covariance matrix
}
\details{
  Add detail for each method
}
\author{
  Ross Bennett
}

