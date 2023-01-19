# CODE1303
## NEWTON-RAPHSON IMPLEMENTATION OF LIKELIHOOD PROBLEM II, 
rm(list=ls()) # Clear the workspace
set.seed(771133)
## BINOMIAL DATA SET (Y "successes" out of n=5 "trials")
dat1 <- 
  matrix(c(0,-0.86,5,1,-0.3,5,3,-0.05,5,5,0.73,5),
         nrow=4,ncol=3,byrow=T)
colnames(dat1) <-c("Y","X","n")
dat1
nit <- 10
miu <- matrix(data=NA, nrow=nit+1,ncol=1)
beta <- matrix(data=NA, nrow=nit+1,ncol=1)
resultnr <- matrix(data=NA, nrow=nit,ncol=2)

miu [1]<- 0.1
beta[1] <- 2.9
for (i in 1:nit)
{
  vc11 <- - sum(5*exp(miu[i]+beta[i]*dat1[,2])/
                  ((1+exp(miu[i]+ beta[i]*dat1[,2]))^2))
  vc22 <- - sum(5*dat1[,2]^2*exp(miu[i]+ beta[i]*dat1[,2])/ 
                  ((1+exp(miu[i]+ beta[i]*dat1[,2]))^2))
  vc12 <- - sum(5*dat1[,2]*exp(miu[i]+ beta[i]*dat1[,2])/ 
                  ((1+exp(miu[i]+ beta[i]*dat1[,2]))^2))
  vcmat <- matrix(c(vc11,vc12,vc12,vc22),nrow=2,ncol=2)
  vcmatinv <- solve(vcmat)
  fd1 <- sum((dat1[,1]-(5*exp(miu[i]+ beta[i]*dat1[,2]))/
                (1+exp(miu[i]+ beta[i]*dat1[,2]))))
  fd2 <- sum(((dat1[,1]*dat1[,2])-(5*dat1[,2]*exp(miu[i]+ 
                                                    beta[i]*dat1[,2]))/(1+exp(miu[i]+ beta[i]*dat1[,2]))))
  
  fd <- matrix(c(fd1,fd2),nrow=2,ncol=1)
  sol0 <- matrix(c(miu[i], beta[i]),nrow=2,ncol=1)
  sol1 <- sol0+(-vcmatinv%*%fd)
  miu[i+1] <-sol1[1]
  beta[i+1] <- sol1[2]
  resultnr[i,] <- c(miu[i+1],beta[i+1])
}
resultnr
-vcmatinv
### USE THE R-FUNCTION OPTIM TO COMPARE WITH THIS PROGRAMME
dat <- data.frame(dat1)
logl <- function(data,par)
{
  with(data,-sum(Y*(par[1]+par[2]*X)-
                   5*log(1+exp(par[1]+par[2]*X))))
}
result <- optim(par=c(0.1,2.9),logl,data=dat,
                hessian=TRUE,method="BFGS")
# IF METHOD NOT INCLUDED, OPTIM USES NELDER-MEAD ALGORITHM
# THIS CALL INCLUDES control=list(trace=1,REPORT=1): PROVIDES THE 
#  PROGRESS OF THE ITERATION
result$par
solve(result$hessian)
