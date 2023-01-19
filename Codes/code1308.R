# CODE1308
#GIBBS SAMPLING ALGORITHM
rm(list=ls()) # Clear the workspace
set.seed(12345)
require(graphics)
dat<-c(45.83,50.37,50.06,51.59,48.43,52.75,42.92,48.57,46.18,50.20)
# SET LENGTH OF CHAIN rep
rep<-10000
result<-matrix(data=NA,nrow=rep,ncol=3)
#INITIALISE mu AND v
#mu<-mean(dat)
#v<-var(dat)
mu<-1
v<-2
n<-length(dat)
# START GIBBS CHAIN
ptm <- proc.time()
for (i in 1:rep)
{
  # GENERATE MIU
  mu<-rnorm(1,mean=mean(dat),sd=sqrt(v/n))
  # COMPUTE SCALE 
  s<-((n-1)*var(dat)+n*(mean(dat)-mu)^2)/n
  # GENERATE V
  v<-n*s/rchisq(1,n)
  result[i,]<-c(i,mu,v)
}
proc.time()-ptm
# END OF GIBBS CHAIN
#mu<-result[1000:rep,2]
mu<-result[,2]

#v<-result[1000:rep,3]
v<-result[,3]

meanmu<-mean(mu)
meanmu
varmu<-var(mu)
varmu
cimu<-quantile(mu,c(0.025,0.975))
cimu
meanv<-mean(v)
meanv
varv<- var(v)
varv
civ<-quantile(v,c(0.025,0.975))
civ
#################################################################
# CODE1308 (cont)
#POST-McMC ANALYSIS
# CODE FOR THE MC VARIANCE BASED ON GEYER

# CHOOSE MU OR V AND PLACE IN VECTOR Y
y<-result[,2] # READS IN mu
#y <- result[,3] # READS IN v
ns <- length(y)
svar<-var(y)*(ns-1)/ns
tau<-1
tausum<-0
for(i in 0:ns)
{
  gamaj<-0.0
  gamak<-0.0
  j<-2*i
  k<-(2*i)+1
  lag1<-j
  lag2<-k
  #USE THE R-FUNCTION ACF TO COMPUTE AUTOCORRELATIONS
  cov1<-acf(y,type="covariance",lag.max=lag1,plot=FALSE)
  cov2<-acf(y,type="covariance",lag.max=lag2,plot=FALSE)
  gamaj<-cov1$acf[lag1+1]
  gamak<-cov2$acf[lag2+1]
  tau<-gamaj+gamak
  if(tau<0)
  {
    break
  }
  tausum<-tausum+tau
}

varch<- -svar+2*tausum
mcvar<-varch/ns
mcvar
efchsize<-svar/mcvar
efchsize
integrautoc<-varch/svar
integrautoc
