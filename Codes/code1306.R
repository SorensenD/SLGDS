# CODE1306
#METROPOLIS-HASTING SINGLE-SITE UPDATING - USES ACF
rm(list=ls()) # Clear the workspace
set.seed(123456)
require(graphics)
y<-c(45.83,50.37,50.06,51.59,48.43,52.75,42.92,48.57,46.18,50.20)
# SET LENGTH OF CHAIN rep
rep<-10000
result<-matrix(data=NA,nrow=rep,ncol=3)

mu<- 15
v<- 10
#CHOOSE TUNING PARAMETERS kmu AND kv
#kmu<-0.19
#kv<-0.1
kmu<-19
kv<-9
acceptv<-0
acceptmu<-0
ptm <- proc.time()
for (i in 1:rep)
{
  #UPDATING THE VARIANCE 
  logYv<-rnorm(1,mean=log(v),sd=sqrt(kv))
  logalfav <-sum((y-mu)^2)/(2*v)-sum((y-mu)^2)/
    (2*exp(logYv))+(length(y)/2)*(log(v)-logYv)
  unif<-runif(1)
  if (unif<exp(logalfav))
  {
    v<-exp(logYv)
    acceptv<-acceptv+1
  }
  #UPDATING THE MEAN
  Ymu<-rnorm(1,mean=mu,sd=sqrt(kmu))
  logalfamu<- sum((y-mu)^2)/(2*v)-sum((y-Ymu)^2)/(2*v)
  unif<-runif(1)
  if (unif<exp(logalfamu))
  {
    mu<-Ymu
    acceptmu<-acceptmu+1
    result[i, ]<-c(i,mu,v)
  }
  else   {
    result[i, ]<-c(i,mu,v)
  }
}
proc.time()-ptm
acceptratiomu<-acceptmu/rep
acceptratiov<-acceptv/rep
# DISCARD THE FIRST 1500 DRAWS
mu<-result[1501:rep,2]
v<-result[1501:rep,3]
meanmus<-mean(mu)
meanmus
varmus<-var(mu)
varmus
cimus<-quantile(mu,c(0.025,0.975))
cimus
meanvs<-mean(v)
meanvs
varvs<- var(v)
varvs
civs<-quantile(v,c(0.025,0.975))
civs
