# CODE1105
rm(list=ls()) # CLEAR WORKSPACE
set.seed(195021)
# LOCAL LINEAR REGRESSION SETS p = 1
p <- 1
x<-seq(from=0, to=2*pi,by=0.2)
f0<-function(x){ 100+sin(2*x)+cos(x/2) }
R2<-2/3
y<-f0(x)+rnorm(n=length(x),sd=sqrt(var(f0(x))*(1-R2)/R2))

z <- 5.12
bz <- matrix(c(1,z),nrow=2)
bzt <- t(bz)
one <- rep(1,length(y))

X <- matrix(data=NA, nrow=length(y), ncol=(p+1))
Wz <- matrix(data=NA, nrow=length(y), ncol=length(y))

X <- cbind(one,x)
Xt <- t(X)
# CHOOSE h
h <- 0.2
d <- x-z
ele <- exp(-(d^2)/(2*h^2))
sm <- sum(ele)
Wz <- diag(ele/sm)
pred_z <- bzt%*%solve(Xt%*%Wz%*%X,Xt%*%Wz%*%y)
cat("x =",z,"\n")
cat("Prediction =", pred_z,"\n")
