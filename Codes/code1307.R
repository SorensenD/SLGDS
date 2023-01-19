# CODE1307
#METROPOLIS-HASTING JOINT UPDATING
rm(list=ls()) # Clear the workspace
set.seed(123456)
require(graphics)
y<-c(45.83,50.37,50.06,51.59,48.43,52.75,42.92,48.57,46.18,50.20)
# SET LENGTH OF CHAIN rep
rep<-10000
result<-matrix(data=NA,nrow=rep,ncol=3)
#INITIALISE mu AND v
mu<-15
v<-10
#CHOOSE TUNING PARAMETERS kmu AND kv
kmu<-0.15
kv<-0.08
accept<-0
ptm <- proc.time()
for (i in 1:rep)
{
  Ymu<-rnorm(1,mean=mu,sd=sqrt(kmu))
  logYv<-rnorm(1,mean=log(v),sd=sqrt(kv))
  logalfa <-sum((y-mu)^2)/(2*v)-sum((y-Ymu)^2)/
    (2*exp(logYv))+(length(y)/2)*(log(v)-logYv)
  unif<-runif(1)
  if (unif<exp(logalfa))
  {
    mu<-Ymu
    v<-exp(logYv)
    result[i, ]<-c(i,mu,v)
    accept<-accept+1
  }
  else
  {
    result[i, ]<-c(i,mu,v)
  }
}
proc.time()-ptm
acceptratio<-accept/rep
mu<-result[1501:rep,2]
v<-result[1501:rep,3]
meanmuj<-mean(mu)
meanmuj
varmuj<-var(mu)
varmuj
cimuj<-quantile(mu,c(0.025,0.975))
cimuj
meanvj<-mean(v)
meanvj
varvj<- var(v)
varvj
civj<-quantile(v,c(0.025,0.975))
civj
