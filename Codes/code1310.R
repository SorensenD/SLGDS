# CODE1310
#BINARY DATA - METROPOLIS-HASTING JOINT UPDATING - PROBIT MODEL
rm(list=ls()) # Clear the workspace
set.seed(123)
#require(graphics)
# CODE USES PACKAGE MVTNORM; IT IS INSTALLED BELOW
#install.packages("mvtnorm", .libPaths()[1])
library(mvtnorm)
#CHOOSE LENGTH OF CHAIN rep
rep<-10000
result<-matrix(data=NA,nrow=rep,ncol=4)
################################################
# CREATE BINARY DATA BASED ON PROBIT MODEL
mu <- -2
beta <- 0.7
nrec <- 30
cov <- rnorm(nrec,2,3) # GENERATE THE COVARIATE
xb <- cov*beta
p1 <- pnorm(mu+xb) # PROBABILITIES ACCORDING TO PROBIT MODEL
#dat1 <- cbind(rbinom(nrec,1,p1),round(cov,digits=0)) # CREATE DATA
d <- data.frame(Y=rbinom(nrec,1,p1),X=round(cov,digits=0)) # CREATE DATA

#colnames(dat1) <- c("Y", "X")
#d <- data.frame(dat1)
#attach(d)
#mean(Y)
mean(d$Y)

#################################################################
# CHOOSE TUNING PARAMETER LAMBDA AND COVARIANCE MATRIX C
lambda<-0.25
c<-matrix(c(1,0.0,0.0,0.1),nrow=2,ncol=2,byrow=T)

# INITIALISE THE MEAN OF THE BIVARIATE DISTRIBUTION
theta<-c(-2,1)

accept<-0
## FUNCTION TO COMPUTE THE LOG-POSTERIOR
logpost <- function(data,par)
{
  theta[1] <-par[1]
  theta[2] <- par[2]
  with(data=d,sum(Y*log(pnorm(theta[1]+theta[2]*X))+(1-Y)*log(1.000001-
                              pnorm(theta[1]+theta[2]*X))))
}
#START MH LOOP
ptm <- proc.time()
for (i in 1:rep)
{
  #SAMPLE PROPOSAL FOR THETA (Ytheta) FROM N(theta,lamdaC)
  Ytheta<- rmvnorm(1,mean=theta,sigma=lambda*c)
  lp1<- logpost(d,c(Ytheta[1],Ytheta[2]))
  lp2<- logpost(d,c(theta[1],theta[2]))
  logalfa<-lp1-lp2
  unif<-runif(1)
  if (unif<exp(logalfa))
  {
    theta[1]<-Ytheta[1]
    theta[2]<-Ytheta[2]
    result[i, ]<-c(i,theta[1],theta[2],pnorm(theta[1]+3*theta[2]))
    accept<-accept+1
  }
  else
  {
    result[i, ]<-c(i,theta[1],theta[2],pnorm(theta[1]+3*theta[2]))
  }
}
proc.time()-ptm
acceptratio<-accept/rep
meanmu<-mean(result[,2])
meanmu
varmu<-var(result[,2])
varmu
cimu<-quantile(result[,2],c(0.025,0.975))
cimu
meanbeta<-mean(result[,3])
meanbeta
varbeta<- var(result[,3])
varbeta
cibeta<-quantile(result[,3],c(0.025,0.975))
cibeta
meanprob <- mean(result[,4])
meanprob
ciprob <- quantile(result[,4],c(0.025,0.975))
ciprob
