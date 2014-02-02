\name{lwShrink}
\alias{lwShrink}
\title{Ledoit-Wolf shrinkage covariance estimate}
\usage{
  lwShrink(x, shrink = NULL)
}
\arguments{
  \item{x}{xts or matrix of asset returns}

  \item{shrink}{shrinkage constant}
}
\value{
  covariance matrix estimate
}
\description{
  Compute the covariance matrix estimate using the
  Ledoit-Wolf shrinkage estimate
}
\note{
  Ported to R from matlab code given at
  http://www.econ.uzh.ch/faculty/wolf/publications.html#9
}
\author{
  Ross Bennett
}
\references{
  TODO
}

