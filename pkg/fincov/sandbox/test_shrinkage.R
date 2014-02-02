library(covariance)
library(PerformanceAnalytics)

data(edhec)
R <- edhec[,1:10]

cov_js1 <- covEstimate(R=R, method="js_shrinkage")
cov_js2 <- corpcor::cov.shrink(x=R)
cov_js3 <- jsShrink(x=R)
all.equal(cov_js1$cov, unclass(cov_js2), check.attributes=FALSE)
all.equal(cov_js1$object, cov_js3, check.attributes=FALSE)

cov_lw1 <- covEstimate(R=R, method="lw_shrinkage")
cov_lw2 <- lwShrink(x=R)
all.equal(cov_lw1$object, cov_lw2)
