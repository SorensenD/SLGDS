# CODE1304
## EM PROBIT MODEL; LIKELIHOOD PROBLEMS II, QUESTION 3ii
rm(list=ls()) # Clear the workspace
set.seed(12371)
### BINOMIAL DATA (Y successes out of n=5 trials)
dat1 <- 
  matrix(c(0,-0.86,5,1,-0.3,5,3,-0.05,5,5,0.73,5),
         nrow=4,ncol=3,byrow=T)
colnames(dat1) <-c("n","X","N")
dat1
nit <- 1000
miu <- 0.1
beta <- 2.9
dat <- data.frame(dat1)
col1 <- matrix(rep(1:1,length(dat1[,1])),
               nrow=length(dat1[,1]),ncol=1)
m <- matrix(c(col1,dat1[,2]),nrow=length(dat1[,1]),ncol=2)
mt <- t(m)
dN <- diag(dat$N)
dn <- diag(dat$n)
dd <- dN-dn
## FUNCTION loglik TO COMPUTE THE LOG-LIKELIHOOD
loglik <- function(data,par)
{
  miu <-par[1]
  beta <- par[2]
  with(data, sum(n*log(pnorm(miu+beta*X))+
                   (N-n)*log(1-pnorm(miu+beta*X))))
}
# CONSTRUCT FUNCTIONS expu0 and expu1 THAT COMPUTE 
#   CONDITIONAL EXPECTATIONS
# ASSUME: Y=1 IF u>0; THEREFORE Pr[Y=1]=F (CDF of std. normal)
expu1 <- function(data,par,i)
{
  miu <- par[1]
  beta <- par[2]
  with(data,(miu+beta*X[i]+
               dnorm(miu+beta*X[i])/pnorm(miu+beta*X[i])))
}
expu0 <- function(data,par,i)
{
  miu <- par[1]
  beta <- par[2]
  with(data,(miu+beta*X[i]-
               dnorm(miu+beta*X[i])/(1-pnorm(miu+beta*X[i]))))
}
sol <- matrix(data=NA, nrow=nit, ncol=2)
e1 <- matrix(data=NA,nrow=length(dat1[,1]),ncol=1)
e0 <- matrix(data=NA,nrow=length(dat1[,1]),ncol=1)
llik <- matrix(data=NA,nrow=nit,ncol=1)
result <- matrix(data=NA,nrow=nit,ncol=3)
for (i in 1:nit)
{
  for (j in 1:length(dat1[,1]))
  {
    e1[j] <- expu1(data=dat,c(miu,beta),j)
    e0[j] <- expu0(data=dat,c(miu,beta),j)
  }
  sol[i,] <- t(solve(mt%*%dN%*%m)%*%(mt%*%dn%*%e1+mt%*%dd%*%e0))
  miu <- sol[i,1]
  beta <- sol[i,2]
  llik[i] <-loglik(dat,c(miu,beta))
  result[i,] <- c(miu,beta,llik[i])
}
tail(result)
### USE THE R-FUNCTION OPTIM TO COMPARE WITH THIS PROGRAMME
logl <- function(data,par)
{
  miu <- par[1]
  beta <- par[2]
  with(data,-sum(n*log(pnorm(miu+beta*X))+
                   (N-n)*log(1-pnorm(miu+beta*X))))
}
result1 <- optim(par=c(0.1,2.9),logl,data=dat,
                 hessian=TRUE,method= "BFGS")
result2 <- optim(par=c(0.1,2.9),logl,data=dat,
                 hessian=TRUE)
result1$par
solve(result1$hessian)
plot(result[,1],xlab="Iterate number",ylab="ML estimate")
