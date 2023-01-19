# CODE1301
rm(list=ls()) # Clear the workspace
set.seed(12371)
# CREATE BINARY DATA
mu <- -2
beta <- 0.7
cov <- rnorm(30,2,3) # GENERATE THE COVARIATE
xb <- cov*beta
p1 <- pnorm(mu+xb) # PROBABILITIES ACCORDING TO PROBIT MODEL
dat1 <- cbind(rbinom(30,1,p1),round(cov,digits=0)) # CREATE DATA
colnames(dat1) <- c("Y", "X")
### END OF GENERATION OF DATA
nit <- 10 # NUMBER OF N-R ITERATIONS
miu <- matrix(data=NA, nrow=nit+1,ncol=1)
beta <- matrix(data=NA, nrow=nit+1,ncol=1)
resultnr <- matrix(data=NA,nrow=nit,ncol=3)
# START VALUES FOR MIU AND BETA
miu [1]<- 0.17
beta[1] <- 0.13
for (i in 1:nit)
{
  vc11 <- - sum(exp(miu[i]+beta[i]*dat1[,2])/((1+exp(miu[i]+ 
            beta[i]*dat1[,2]))^2))
  vc22 <- - sum(dat1[,2]^2*exp(miu[i]+ beta[i]*dat1[,2])/ 
                  ((1+exp(miu[i]+ beta[i]*dat1[,2]))^2))
  
  vc12 <- - sum(dat1[,2]*exp(miu[i]+ beta[i]*dat1[,2])/ 
                  ((1+exp(miu[i]+beta[i]*dat1[,2]))^2)) 
  
  vcmat <- matrix(c(vc11,vc12,vc12,vc22),nrow=2,ncol=2)
  vcmatinv <- solve(vcmat)
  fd1 <- sum((dat1[,1]-(exp(miu[i]+ beta[i]*dat1[,2]))/
                (1+exp(miu[i]+beta[i]*dat1[,2])))) 
  
  fd2 <- sum(((dat1[,1]*dat1[,2])-(dat1[,2]*exp(miu[i]+ 
         beta[i]*dat1[,2]))/(1+exp(miu[i]+beta[i]*dat1[,2]))))
  
  fd <- matrix(c(fd1,fd2),nrow=2,ncol=1)
  sol0 <- matrix(c(miu[i], beta[i]),nrow=2,ncol=1)
  sol1 <- sol0+(-vcmatinv%*%fd)
  miu[i+1] <-sol1[1]
  beta[i+1] <- sol1[2]
  resultnr[i,] <- c(i,sol1[1],sol1[2])
}
resultnr
# ASYMPTOTIC COVARIANCE MATRIX
-vcmatinv

## COMPUTE PROBABILITIES THAT Y=1, GIVEN X = -3, 1, 9
p1<-exp(miu[i+1]+beta[i+1]*(-3))/(1+exp(miu[i+1]+beta[i+1]*(-3)))
p1
p2<-exp(miu[i+1]+beta[i+1]*1)/(1+ exp(miu[i+1]+beta[i+1]*1))
p2
p3<-exp(miu[i+1]+beta[i+1]*9)/(1+ exp(miu[i+1]+beta[i+1]*9))
p3
##############################################
### USE THE R-FUNCTION OPTIM TO COMPARE WITH THIS PROGRAMME
# CODE1301 (cont)
dat <- data.frame(dat1)
logl <- function(data,par)
{
  with(data,-sum(Y*(par[1]+par[2]*X)-log(1+exp(par[1]+par[2]*X))))
}
result<-optim(par=c(-3.5,0.05),logl,data=dat,
              hessian=TRUE,method="BFGS")
# IF METHOD IS NOT INCLUDED IN THE CALL, OPTIM USES 
#  THE NELDER-MEAD ALGORITHM
# result <- optim(par=c(-3.5,0.05),logl,data=dat,hessian=TRUE,
#  method="BFGS", control=list(trace=1,REPORT=1))
# THIS CALL INCLUDES  control=list(trace=1,REPORT=1) WHICH 
#  PROVIDES THE PROGRESS OF THE ITERATION
result$par
solve(result$hessian)
##############################################
########### USE GLM in R
d <- data.frame(dat1)
logreg <- glm(d$Y~d$X, data=d, family=binomial(link="logit"))
summary(logreg)
pred <- predict(logreg,d,type="response")
logscore <- sum(d$Y*log(pred)+(1-d$Y)*log(1-pred)) # LogLik of 
#          full model evaluated at ML estimates
-2*logscore # DEVIANCE FULL MODEL
# NULL MODEL
logregnull <- glm(d$Y~1,data=d,family=binomial(link="logit")) # 
summary(logregnull)
prednull <- predict(logregnull,d,type="response")
logscorenull <- sum(d$Y*log(prednull)+(1-d$Y)*log(1-prednull)) 
# LogLik of null model evaluated at ML estimate
-2*logscorenull # DEVIANCE NULL MODEL
