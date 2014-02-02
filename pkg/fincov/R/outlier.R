#' Outlier Detection
#' 
#' Detect outliers using classical and robust estimates of mahalanobis distance.
#' @param R
#' @param method method to compute robust estimate of location and scalde used
#' in computing the robust mahalanobis distance.
#' @param quantile 
#' @export
detectOutliers <- function(R, 
                           method=c("mcd", "m_estimate", "mm_estimate", 
                                    "mve", "ogk", "sde", "s_estimate"), 
                           quantile=0.975,
                           threshold=NULL){
  method = match.arg(method)
  
  if(is.null(threshold)){
    threshold <- sqrt(qchisq(quantile, ncol(R)))
  }
  
  # compute the classical mahalanobis distance estimates
  classical_md <- mahalanobis(R, colMeans(R), cov(R))
  classical_outliers <- R[which(sqrt(classical_md) > threshold),]
  
  # compute robust estimates of location and scale
  cov_rob <- covEstimate(R, method)
  robust_md <- mahalanobis(R, cov_rob$object@center, cov_rob$object@cov)
  robust_outliers <- R[which(sqrt(robust_md) > threshold),]
  
  structure(list(classical_mahalanobis=classical_md, 
                 robust_mahalanobis=robust_md,
                 classical_outliers=classical_outliers,
                 robust_outliers=robust_outliers,
                 threshold=threshold),
            class="outlier_mahalanobis")
}

# need some plot methods for outliers
