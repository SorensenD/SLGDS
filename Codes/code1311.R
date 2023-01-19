# CODE1311
# BINARY DATA; GIBBS WITH DATA AUGMENTATION PROBIT MODEL
# LIABILITIES SAMPLED IN ONE GO AFTER JA  pg 241
rm(list=ls()) # Clear the workspace
set.seed(123)
require(graphics)
# THE CODE WILL USE THE PACKAGE MVTNORM; IT IS INSTALLED BELOW
#install.packages("mvtnorm", .libPaths()[1])
library(mvtnorm)
#CHOOSE LENGTH OF CHAIN rep
rep<-10000
result<-matrix(data=NA,nrow=rep,ncol=3)
################################################
# CREATE BINARY DATA
mu <- -2
beta <- 0.7
nrec <- 30
cov <- rnorm(nrec,2,3) # GENERATE THE COVARIATE
xb <- cov*beta
p1 <- pnorm(mu+xb) # COMPUTE PROBABILITIES FOR PROBIT MODEL
#p1 <- rbeta(30,2,2)
# CREATE DATA:
dat1 <- cbind(rbinom(nrec,1,p1),round(cov,digits=0))
colnames(dat1) <- c("Y", "X")
d <- data.frame(dat1)
#head(d)
attach(d)
mean(Y)
##################################################################
One<-rep(1,nrec)
X<-matrix(c(One,d[,2]),nrow=nrec,ncol=2)
Y<-matrix(d[,1],nrow=length(d[,1]),ncol=1)
#INITIALISE THETA (THE VECTOR WITH MU AND BETA)
theta<-matrix(data=0,nrow=2,ncol=1)
# INITIALISE VECTOR u
u<-rep(0,nrec)
#COMPUTE XTX INVERSE
xtxinv<-solve(t(X)%*%X)
#START GIBBS LOOP
ptm <- proc.time()
for (i in 1:rep)
{
  thetahat<-xtxinv%*%t(X)%*%u
  #SAMPLE THETA FROM MVN(thetahat,xtxinv)
  theta<- t(rmvnorm(1,mean=thetahat,sigma= xtxinv))
  # SAMPLE LIABILITIES u_j FROM A TN(beta+Uf,1)
  
  av<-X%*%theta # MEAN OF UNTRUNCATED NORMAL(mu+U[j,]alfa,1)
  cutoff<-pnorm(-av)
  interm<-(cutoff*(1-Y)+(1-cutoff)*Y)*runif(length(Y))+cutoff*Y
  interm[interm==1]<-0.999
  u<-qnorm(interm)+av
  # END SAMPLING LIABILITIES u
  result[i,]<-c(i,t(theta))
}
proc.time() - ptm
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
