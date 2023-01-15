# CODE1107
rm(list=ls()) # CLEAR WORKSPACE
set.seed(195021)
library(MASS)
N <- 100
# GENERATE DATA 
x<-seq(from=1, to=2.2*pi,length=N)
signal <- cos(1.5*x)+ exp(-0.4*x)
noise <- rnorm(N,0,0.25)
y <- signal + noise
#setwd("C:/Users/au223137/Dropbox/Rsessions/MarkDown")
#pdf("C:/Users/au223137/Dropbox/Rsessions/MarkDown/Figures/
#    KregEx11.pdf")
#plot(x,y)

lambda <- 0.01
X <- cbind(1,x)
RHS <- crossprod(X,y)
LHS <- crossprod(X)
LHS[-1,-1] <- LHS[-1,-1]+diag(c(rep(1,1)))*lambda # assumes 
#          one covariate x
####  Classical LS solution 
bh1 <- solve(LHS,RHS)
bh1
yhatclassic <- X%*%bh1
yhatclassic[1:6]
yhatclassic[7:10]
#########################################################
####  kernelised (dual) solution with kernel matrix XX'
K <- x%*%t(x) # LINEAR KERNEL
X <- cbind(1,K)
RHS <- crossprod(X,y)
LHS <- crossprod(X)
LHS[-1,-1] <- LHS[-1,-1]+K*lambda
diag(LHS) <- diag(LHS) + c(0,rep(1e-8,N))
sol <- as.matrix(solve(LHS,RHS))
yhatkernellin <- sol[1]+K%*%sol[-1]
yhatkernellin[1:6]
yhatkernellin[7:10]
plot(x,y)
lines(x,yhatkernellin,col="blue")
alfa <- sol[-1]
bhkernellin <- sum(alfa*x)
bhkernellin
muhatkernellin <- sol[1]
muhatkernellin
########################################################
#######  GAUSSIAN KERNEL SOLUTION ######################
#     CONSTRUCT GAUSSIAN KERNEL
#w <- matrix(data=NA,nrow=length(x),ncol=length(x))
K <- matrix(data=NA,nrow=length(x),ncol=length(x))

d <- as.matrix(dist(x))^2
# CHOOSE h and lambda
h <- 0.7
lambda <- 0.5

K <- exp(-(1/(2*h^2))*d)

X <- cbind(1,K)
RHS <- crossprod(X,y)
LHS <- crossprod(X)
LHS[-1,-1] <- LHS[-1,-1]+K*lambda

diag(LHS) <- diag(LHS)+c(0,rep(1e-8,N))

sol <- solve(LHS,RHS)
fgaus <- sol[1]+K%*%sol[-1]
alfa <- sol[-1]
#setwd("C:/Users/au223137/Dropbox/Rsessions/MarkDown")
#pdf("C:/Users/au223137/Dropbox/Rsessions/MarkDown/Figures/
#     KregEx12.pdf")

plot(x,y)
lines(x,fgaus,col="red")

