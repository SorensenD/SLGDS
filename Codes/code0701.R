# CODE0701
# AN EXAMPLE USING SIMULATED CORRELATED X
rm(list=ls()) # CLEAR WORKSPACE
set.seed(3711)
#install.packages("glmnet", .libPaths()[1])
library(glmnet)
n <- 100
p <- 100
X <- matrix(rnorm(p*n),ncol=p)
X <- X*0.8 + X[,1]*0.3 # GENERATE CORRELATED COVARIATES
X <- scale(X)*sqrt(n)/sqrt(n-1)
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

