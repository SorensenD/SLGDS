# CODE1106
# COMPUTE HAT MATRIX FOR THE LOCAL LINEAR REGRESSION
rm(list=ls()) # CLEAR WORKSPACE
set.seed(195021)

x<-seq(from=0, to=3*pi,by=0.05)

f0<-function(x){ 100+sin(2*x)+cos(x/2) }
noise <- rnorm(n=length(x),sd=sqrt(var(f0(x))*0.5))
y<- f0(x) + noise

one <- rep(1,length(x))

X <- matrix(data=NA, nrow=length(y), ncol=(2))
W <- matrix(data=NA, nrow=length(y), ncol=length(y))
HatMat <- matrix(data=NA, nrow=length(x), ncol=length(x))

h <- seq(from = 0.1, to = 1.5, by = 0.025)
GCV <- rep(0,length(h))
LOOCV <- rep(0,length(h))
form40 <- rep(0,length(h))
dsq <- rep(0,length(h))

tr <- rep(0,length(h))
Trmse <- rep(0,length(h))

j <- seq(1:length(x))

X <- cbind(one,x)
Xt <- t(X)

Hat <- function(j,h,k,x,X){
  d <- x - x[j] 
  ele <- exp(-(d^2)/(2*h[k]^2))
  sm <- sum(ele)
  w <- diag(ele)/sm
  X[j,]%*%solve(Xt%*%w%*%X)%*%Xt%*%w
}

for (k in 1:length(h)) {
  HatMat<- t(sapply(1:length(x),Hat,h,k,x,X))
  predy <- HatMat %*% y
  tr[k] <- sum(diag(HatMat))
  Trmse[k] <- mean((predy - y) ^ 2)
  dsq[k] <- (1 - mean(diag(HatMat))) ^ 2
  GCV[k] <- Trmse[k] / dsq[k] # generalised LOOCV
  omd <- 1 - diag(HatMat)
  LOOCV[k] <- mean(((predy - y) / omd) ^ 2) # LOOCV
}

h[which.min(LOOCV)]
h[which.min(GCV)]
min(LOOCV)
min(GCV)
