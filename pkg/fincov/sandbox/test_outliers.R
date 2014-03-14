# test the outlier detection
# Ch 6.7 in Scherer and Martin

library(lattice)

file <- "/Users/rossbennett/Documents/UW Comp Fin/Portfolio-Optimization-Asset-Management/Material/chap6/normal.vs.hectic.ts.csv"

dat <- read.csv(file, header=TRUE)
dat <- xts(dat[,3:6], order.by=as.Date(dat[,1]))
dat <- dat[-(1:60),]
head(dat)

xyplot(dat)

outliers <- detectOutliers(dat)
outliers
x <- outliers$threshold

m_dat <- data.frame(robust=outliers$robust_mahalanobis, classical=outliers$classical_mahalanobis)
idx <- 1:nrow(dat)

xyplot(m_dat, layout=c(2,1))

# par(mfrow=c(1,2))
# par(mfrow=c(1,1))

# robust:::distancePlot.covfm

# square root of classic mahalanobis distances
cl_md <- sqrt(outliers$classical_mahalanobis)

# square root of robust mahalanobis distances
rob_md <- sqrt(outliers$robust_mahalanobis)


# par(mfrow=c(1,2))
layout(matrix(c(1,2, 1,2), 2, 2, byrow = TRUE))
layout.show()
plot(rob_md, ylab="Square Root of Mahalanobis Distance")
text(x=which(rob_md > x), y=rob_md[which(rob_md > x)], labels=which(rob_md > x), pos=4, cex=0.8)
abline(h=x, lty=2)

plot(cl_md, ylab="Square Root of Mahalanobis Distance")
text(x=which(cl_md > x), y=cl_md[which(cl_md > x)], labels=which(cl_md > x), pos=4, cex=0.8)
abline(h=x, lty=2)
# par(mfrow=c(1,1))

# xyplot with vertical lines at outlier dates
xyplot(dat, panel=function(...){
  panel.abline(v=as.Date(cl_outlier_dates))
  panel.xyplot(...)
  })

# actual data of the outlier dates
dat[as.Date(cl_outlier_dates)]

# xyplot with vertical lines at outlier dates
xyplot(dat, panel=function(...){
  panel.abline(v=as.Date(rob_outlier_dates))
  panel.xyplot(...)
})

# actual data of the outlier dates
dat[as.Date(rob_outlier_dates)]


# outliers$classical_mahalanobis
# length(outliers$classical_mahalanobis)
# dim(dat)

