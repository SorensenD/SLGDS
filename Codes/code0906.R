# CODE0906
# GAUSSIAN FULL SIB-FAMILY MODEL; SINGLE-SITE GIBBS SAMPLING
# AS AN APPROXIMATION TO THE TRUE PROBIT MODEL
rm(list=ls()) # Clear the workspace
set.seed(12345)
require(graphics)
# GENERATE CORRELATED (FULL-SIBS DATA 
#install.packages("MCMCpack", .libPaths()[1])
#install.packages("mvtnorm", .libPaths()[1])
library(MCMCpack)
# INITIALISE PARAMETERS
p0 <- 0.15
mu <- qnorm(p0)
iccfs<-0.15 #INTRACLASS CORRELATION FS
# VARIANCE BETWEEN FAMILIES: iccfs /(1- iccfs)
vfs <- iccfs/(1-iccfs)
nfs<-1000 # NUMBER OF FULL-SIB FAMILIES

fs<-3 #FULL-SIB FAMILY SIZE
N<-nfs*fs

c<-0
##########################################################
####  GENERATE BINARY RECORDS Y
f<-rnorm(nfs,mean=0,sd=sqrt(vfs))
p <- pnorm(mu+f)
y <- rbinom(N,1,rep(p,each=fs))
w <- rep(1:nfs,each=fs)
d<-data.frame(w,y)
family <- as.factor(w)
Z<-model.matrix(~0+family)
##########################################################
ztz<-t(Z)%*%Z
rep<-2000
resultap<-matrix(data=NA,nrow=rep,ncol=5)
transf<-matrix(data=NA,nrow=rep,ncol=4)

#INITIALISE THE VECTOR OF FAMILIY EFFECTS fe
# (Not to be confused with the TRUE family effects f)
fe<-rep(0,nfs)
# INITIALISE BETWEEN FAMILY VARIANCE COMPONENT vf
vf<-1
# INITILISE RESIDUAL VARIANCE
ve<-1
# INITIALISE k
k<-ve/vf
# INITIALISE THE MEAN
mu<-0
sumpyinvt<-0
#START GIBBS LOOP NORMAL MODEL
ptm <- proc.time()
for (i in 1:rep)
{
  print(i)
  # SAMPLE mu
  meanmu<-sum(y-Z%*%fe)/(nfs*fs)
  mu<-rnorm(1,mean=meanmu,sd=sqrt(ve/(nfs*fs)))
  # SAMPLE FAMILY EFFECTS f
  varf<-(k+fs)^(-1)
  fmean<- varf*(t(Z)%*%(y-mu))
  fe<-rnorm(nfs,mean=fmean, sd=sqrt(varf*ve))
  #SAMPLE vf
  #COMPUTE SCALE
  ftf<-sum(fe*fe)
  vfx<-ftf/rchisq(1,nfs-2)
  vf<-as.numeric(vfx)
  # SAMPLE ve
  # COMPUTE SCALE
  e<-(y-mu-Z%*%fe)
  ete<-t(e)%*%e
  vex<-ete/rchisq(1,N-2)
  ve<-as.numeric(vex)
  k<-ve/vf
  her <- (2*vf)/(vf+ve)
  resultap[i,]<-c(i,mu,vf,ve,her)
  # TRANSFORM TO PARAMETERS IN UNDERLYING SCALE
  mut <- qnorm(mu)
  vft <- vf/(dnorm(mut)**2)
  hert <- (2*vft)/(vft+1)
  transf[i,] <- c(i,mut,vft,hert)
}
proc.time()-ptm
apply(resultap[,2:5],2,mean)
apply(transf[,2:4],2,mean)

