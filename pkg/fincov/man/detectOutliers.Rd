\name{detectOutliers}
\alias{detectOutliers}
\title{Outlier Detection}
\usage{
  detectOutliers(R,
    method = c("mcd", "m_estimate", "mm_estimate", "mve", "ogk", "sde", "s_estimate"),
    quantile = 0.975, threshold = NULL)
}
\arguments{
  \item{R}{}

  \item{method}{method to compute robust estimate of
  location and scalde used in computing the robust
  mahalanobis distance.}

  \item{quantile}{}
}
\description{
  Detect outliers using classical and robust estimates of
  mahalanobis distance.
}

