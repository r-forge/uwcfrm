
library(rrcov)
library(covariance)

# data(hbk)
# hbk.x <- data.matrix(hbk[, 1:3])

file <- "/Users/rossbennett/Documents/UW Comp Fin/Portfolio-Optimization-Asset-Management/Material/chap6/normal.vs.hectic.ts.csv"

dat <- read.csv(file, header=TRUE)
dat <- xts(dat[,3:6], order.by=as.Date(dat[,1]))
dat <- dat[-(1:60),]
head(dat)

##### CovClassic #####
cov_classic <- covEstimate(dat, "classic")
cov_classic_rr <- CovClassic(dat)
all.equal(cov_classic$object, cov_classic_rr, check.attributes=FALSE)


##### CovMcd #####
cov_mcd <- covEstimate(dat, "mcd", control=list(nsamp=300))
cov_mcd_rr <- CovMcd(x=dat, nsamp=300)
all.equal(cov_mcd$object, cov_mcd_rr, check.attributes=FALSE)


##### CovMest #####
cov_mest <- covEstimate(dat, "m_estimate")
cov_mest_rr <- CovMest(x=dat)
all.equal(cov_mest$object, cov_mest_rr, check.attributes=FALSE)


##### CovMMest #####
cov_mmest <- covEstimate(dat, "mm_estimate")
cov_mmest_rr <- CovMMest(x=dat)
all.equal(cov_mmest$object, cov_mmest_rr, check.attributes=FALSE)


##### CovMve #####
cov_mve <- covEstimate(dat, "mve")
cov_mve_rr <- CovMve(x=dat)
all.equal(cov_mve$object, cov_mve_rr, check.attributes=FALSE)


##### CovOgk #####
cov_ogk <- covEstimate(dat, "ogk")
cov_ogk_rr <- CovOgk(x=dat)
all.equal(cov_ogk$object, cov_ogk_rr, check.attributes=FALSE)


##### CovSde #####
cov_sde <- covEstimate(dat, "sde")
cov_sde_rr <- CovSde(x=dat)
all.equal(cov_sde$object, cov_sde_rr, check.attributes=FALSE)


##### CovSest #####
cov_sest <- covEstimate(dat, "s_estimate")
cov_sest_rr <- CovSest(x=dat)
all.equal(cov_sest$object, cov_sest_rr, check.attributes=FALSE)


