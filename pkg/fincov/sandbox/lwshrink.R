
library(PerformanceAnalytics)
data(edhec)

lwShrink <- function(x, shrink=NULL){
  # port of matlab code to R from http://www.econ.uzh.ch/faculty/wolf/publications.html#9
  # Ledoit, O. and Wolf, M. (2004).
  # Honey, I shrunk the sample covariance matrix.
  # Journal of Portfolio Management 30, Volume 4, 110-119.
  
  # % de-mean returns
  # [t,n]=size(x);
  # meanx=mean(x);
  # x=x-meanx(ones(t,1),:);
  
  # de-mean returns
  t <- nrow(x)
  n <- ncol(x)
  meanx <- colMeans(x)
  x <- x - matrix(rep(meanx, t), ncol=n, byrow=T)
  
  # % compute sample covariance matrix
  # sample=(1/t).*(x'*x);
  
  # compute sample covariance matrix
  sample <- (1/t) * (t(x) %*% x)
  
  # % compute prior
  # var=diag(sample);
  # sqrtvar=sqrt(var);
  # rBar=(sum(sum(sample./(sqrtvar(:,ones(n,1)).*sqrtvar(:,ones(n,1))')))-n)/(n*(n-1));
  # prior=rBar*sqrtvar(:,ones(n,1)).*sqrtvar(:,ones(n,1))';
  # prior(logical(eye(n)))=var;
  
  # compute prior
  var <- matrix(diag(sample), ncol=1)
  sqrtvar <- sqrt(var)
  tmpMat <- matrix(rep(sqrtvar, n), nrow=n)
  rBar <- (sum(sum(sample / (tmpMat * t(tmpMat))))-n) / (n * (n - 1))
  prior <- rBar * tmpMat * t(tmpMat)
  diag(prior) <- var
  
  if(is.null(shrink)){
    # % what we call pi-hat
    # y=x.^2;
    # phiMat=y'*y/t - 2*(x'*x).*sample/t + sample.^2;
    # phi=sum(sum(phiMat));
    
    # what we call pi-hat
    y <- x^2
    phiMat <- t(y) %*% y/t - 2 * (t(x) %*% x) * sample / t + sample^2
    phi <- sum(phiMat)
    
    # % what we call rho-hat
    # term1=((x.^3)'*x)/t;
    # help = x'*x/t;
    # helpDiag=diag(help);
    # term2=helpDiag(:,ones(n,1)).*sample;
    # term3=help.*var(:,ones(n,1));
    # term4=var(:,ones(n,1)).*sample;
    # thetaMat=term1-term2-term3+term4;
    # thetaMat(logical(eye(n)))=zeros(n,1);
    # rho=sum(diag(phiMat))+rBar*sum(sum(((1./sqrtvar)*sqrtvar').*thetaMat))
    
    # what we call rho-hat
    term1 <- (t(x^3) %*% x) / t
    help <- t(x) %*% x / t
    helpDiag <- matrix(diag(help), ncol=1)
    term2 <- matrix(rep(helpDiag, n), ncol=n, byrow=FALSE) * sample
    term3 <- help * matrix(rep(var, n), ncol=n, byrow=FALSE)
    term4 <- matrix(rep(var, n), ncol=n, byrow=FALSE) * sample
    thetaMat <- term1 - term2 - term3 + term4
    diag(thetaMat) <- 0
    rho <- sum(diag(phiMat)) + rBar * sum(sum(((1 / sqrtvar) %*% t(sqrtvar)) * thetaMat))
    
    # % what we call gamma-hat
    # gamma=norm(sample-prior,'fro')^2
    
    # what we call gamma-hat
    gamma <- norm(sample - prior, "F")^2
    
    # % compute shrinkage constant
    # kappa=(phi-rho)/gamma;
    # shrinkage=max(0,min(1,kappa/t))
    
    # compute shrinkage constant
    kappa <- (phi - rho) / gamma
    shrinkage <- max(0, min(1, kappa / t))
    
    # % compute the estimator
    # sigma=shrinkage*prior+(1-shrinkage)*sample
  } else {
    shrinkage <- shrink
  }
  # compute the estimator
  sigma <- shrinkage * prior + (1 - shrinkage) * sample
  return(sigma)
}

# returns data
R <- as.matrix(head(edhec[, 1:4]))

lwShrink(R)

