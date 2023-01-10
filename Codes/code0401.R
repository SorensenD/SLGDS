# CODE0401
# CODE FOR THE MC VARIANCE TESTED ON 
#A FIRST-ORDER AUTOREGRESSIVE PROCESS WITH CORRELATION RHO
rm(list=ls()) # Clear the workspace
set.seed(1237)
#GENERATE DATA: AUTOREGRESSIVE PROCESS WITH CORRELATION RHO
ns <- 10000
y <- rep(NA,ns)
rho<-0.8
sum <- 0
y[1] <- rnorm(1,0,1)
for(i in 2:ns)
{
  y[i] <- rho*y[i-1] + rnorm(1,0,1)
  sum <- sum + y[i]*y[i-1]
}
cov <- sum/ns
rhohat <- cov/var(y)
muhat <- mean(y)
gama0 <- var(y)
rhohat
# CODE FOR THE MC VARIANCE BASED ON GEYER
svar<-var(y)*(ns-1)/ns
tau<-1
tausum<-0
ptm <- proc.time()
for(i in 0:ns)
{
  gamaj<-0.0
  gamak<-0.0
  j<-2*i
  k<-(2*i)+1
  # FASTER CODE: JUMP THE LOOP AND USE FUINCTION ACF
  #  for (ii in 1:(ns-j))
  #  {
  #    cov<-(y[ii]*y[ii+j]-mean(y)*(y[ii]+y[ii+j])+mean(y)*mean(y))
  #    gamaj<-gamaj+cov
  #  }
  #  for(ii in 1:(ns-k))
  #  {
  #    gamak<-gamak+(y[ii]*y[ii+k]-mean(y)*(y[ii]+y[ii+k])+
  #        mean(y)*mean(y))
  #  }
  #  gamaj<-gamaj/ns
  #  gamak<-gamak/ns
  
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
proc.time()-ptm
muhat
gama0
varch<- -svar+2*tausum
varch
mcvar<-varch/ns
mcvar
efchsize<-svar/mcvar
efchsize
integrautoc<-varch/svar
integrautoc
