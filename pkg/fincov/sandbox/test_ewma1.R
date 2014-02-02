library(PerformanceAnalytics)

# Should create my own data set for the covariance package
data(edhec)
R <- edhec[, 1:2]

n <- nrow(R)
p <- ncol(R)

startup_period <- 10
training_period <- 36

##### Object Oriented Approach #####

set.seed(123)
robust_model <- ExponentialSmoothing(R, NULL, startup_period, training_period, "robust")
names(robust_model)

# extractor functions
# getLambda
getError
getSmoothedSeries
getTrainingData
getTestData

# robust exponential smoothing
getRobustScaleEstimate
getRobustCleanedSeries
getRobustSmoothedSeries


set.seed(123)
classic_model <- ExponentialSmoothing(R, NULL, startup_period, training_period, "classic")
names(classic_model)



