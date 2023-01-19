# CODE1302
# LikIIQ2 EM algorithm with probit model
rm(list=ls()) # Clear the workspace
set.seed(12371)
# CREATE BINARY DATA
mu <- -2
beta <- 0.7
cov <- rnorm(30,2,3) # GENERATE THE COVARIATE
xb <- cov*beta
p1 <- pnorm(mu+xb) # PROBABILITIES ACCORDING TO PROBIT MODEL
#p1 <- rbeta(30,2,2)
dat1 <- cbind(rbinom(30,1,p1),round(cov,digits=0)) # CREATE DATA
colnames(dat1) <- c("Y", "X")
d <- data.frame(dat1)
attach(d)
m <- cbind(1,X)
mt <- t(m)
mtm <- mt%*%m
miu <- -1
beta <- 1
# CONSTRUCT FUNCTION expu THAT COMPUTES CONDITIONAL EXPECTATIONS
# HERE WE ASSUME THAT Y=1 IF u>0
expu <- function(data,par,i)
{
  miu <- par[1]
  beta <- par[2]
  if(Y[i]==1)
  {
    with(data,miu+beta*X[i]+
           (dnorm(miu+beta*X[i])/pnorm(miu+beta*X[i])))
  } else
  {
    with(data,miu+beta*X[i]-
           (dnorm(miu+beta*X[i])/(1-pnorm(miu+beta*X[i]))))
  }
}
## FUNCTION loglik TO COMPUTE THE LOG-LIKELIHOOD
loglik <- function(data,par)
{
  miu <-par[1]
  beta <- par[2]
  with(data,sum(Y*log(pnorm(miu+beta*X))+
                  (1-Y)*log(1-pnorm(miu+beta*X))))
}
iter <- 100
euvec <- matrix(data=NA, nrow=30,ncol=1)
sol <- matrix(data=NA, nrow=iter, ncol=2)
llik <- matrix(data=NA,nrow=iter,ncol=1)
result <- matrix(data=NA,nrow=iter,ncol=3)
# PLACE 30 CONDITIONAL EXPECTATIONS IN euvec AND ITERATE
for (i in 1:iter)
{
  for (j in 1:length(X))
  {
    euvec[j] <- expu(d,c(miu,beta),j)
  }
  sol[i,] <- t(solve(mtm)%*%mt%*%euvec)
  miu <- sol[i,1]
  beta <- sol[i,2]
  llik[i] <-loglik(d,c(miu,beta))
  result[i,] <- c(miu,beta,llik[i])
}
# FINAL ITERATES OF THE EM ALGORITHM
tail(result)
###########################################################
### USE THE R-FUNCTION OPTIM TO COMPARE WITH THIS PROGRAMME
logl <- function(data,par)
{
  miu <- par[1]
  beta <- par[2]
  with(data,-sum(Y*log(pnorm(miu+beta*X))+
                   (1-Y)*log(1-pnorm(miu+beta*X))))
}
result1<-optim(par=c(-1,1),logl,data=d,hessian=TRUE,method="BFGS")
# ML ESTIMATES USING OPTIM
result1$par
