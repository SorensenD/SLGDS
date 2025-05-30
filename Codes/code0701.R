# CODE0701
rm(list=ls()) # CLEAR WORKSPACE
set.seed(3037)
#install.packages("glmnet", .libPaths()[1])
library(glmnet)
n <- 100
p <- 100
rho <- 0.6
# GENERATE A MATRIX OF CORRELATED COVARIATES X
X <- matrix(data=0,nrow=n,ncol=p)
for (i in 1:n){
  X[i,1] <- rnorm(1)
  for (j in 1:(p-1)){
    X[i,j+1] <- X[i,j]*rho + rnorm(1)
  }
}
X <- scale(X)
beta <- rep(0,p)
betac <- rep(0,p)
beta <- sample(0:1,p,replace=TRUE,prob=c(2,1))
length(which(beta!=0))
y <- X%*%beta + rnorm(n,sd=0.4)
y <- y - mean(y)
for(i in 1:p){ betac[i]=coef(lm(y~X[,i]))[2]}
lambda=max(abs(betac))*.1
lambda
niter <- 100
bL=matrix(nrow=niter,ncol=p)
bL[1,]=0 # initial lasso estimates set to zero
r1 <- y-mean(y)
for (i in 2:niter) {
  for(j in 1:p){
    r0 <- r1+X[,j]*bL[i-1,j]
    bLS <- crossprod(X[,j],r0)/n # LEAST SQUARES ESTIMATE
    if (abs(bLS) >= lambda){bL[i,j]<-sign(bLS)*(abs(bLS)-lambda)
    } else{
      bL[i,j] <- 0
    }
    r1 <- r0-X[,j]*bL[i,j]
  }
}
fm=glmnet(y=y,x=X,alpha=1,lambda=lambda) 
# alpha=1: LASSO; alpha=0: RIDGE; 0<alpha<1: ELASTIC NET
# Number covariates included with GLMNET:
length(which(fm$beta!=0)) 
# Number included with present code
length(which(bL[niter,]!=0)) 
# PRINT A FEW ESTIMATES OBTAINED WITH GLMNET 
round(fm$beta[which(fm$beta!=0)][1:7],3)
# AND THE SAME WITH PRESENT CODE
round(bL[niter,which(bL[niter,]!=0)][1:7],3)

