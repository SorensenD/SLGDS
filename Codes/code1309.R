# CODE1309
#BINARY DATA - METROPOLIS-HASTING JOINT UPDATING - LOGISTIC MODEL
rm(list=ls()) # Clear the workspace
set.seed(123)
require(graphics)
# THE CODE WILL USE THE PACKAGE MVTNORM; IT IS INSTALLED BELOW
#install.packages("mvtnorm", .libPaths()[1])
library(mvtnorm)
#CHOOSE LENGTH OF CHAIN rep
rep<-10000
result<-matrix(data=NA,nrow=rep,ncol=4)
################################################
# CREATE BINARY DATA
mu <- -2
beta <- 0.7
nrec <- 30
cov <- rnorm(nrec,2,3) # GENERATE THE COVARIATE
xb <- cov*beta
p1 <- pnorm(mu+xb) # PROBABILITIES ACCORDING TO PROBIT MODEL
dat1 <- cbind(rbinom(nrec,1,p1),round(cov,digits=0)) # CREATE DATA
colnames(dat1) <- c("Y", "X")
d <- data.frame(dat1)
attach(d)
mean(Y)
# CHOOSE TUNING PARAMETER LAMBDA AND COVARIANCE MATRIX C
lambda<-1
c<-matrix(c(1,0.0,0.0,0.1),nrow=2,ncol=2,byrow=T)

# INITIALISE THE MEAN OF THE BIVARIATE DISTRIBUTION
theta<-c(-2,1)
accept<-0
## FUNCTION TO COMPUTE THE LOG-POSTERIOR
logpost <- function(data,par)
{
  theta[1] <-par[1]
  theta[2] <- par[2]
  with(data=d,sum(Y*( theta[1] + theta[2]*X)-log(1+exp(theta[1] + 
                                                         theta[2]*X))))
}
#START MH LOOP
ptm <- proc.time()
for (i in 1:rep)
{
  #SAMPLE PROPOSAL FOR THETA (Ytheta) FROM N(theta,lamdaC)
  Ytheta<- rmvnorm(1,mean=theta,sigma=lambda*c)
  logalfa<-logpost(d,c(Ytheta[1],Ytheta[2])) - 
    logpost(d,c(theta[1],theta[2]))
  unif<-runif(1)
  if (unif<exp(logalfa))
  {
    theta[1]<-Ytheta[1]
    theta[2]<-Ytheta[2]
    result[i, ]<-c(i,theta[1],theta[2],exp(theta[1]+3*theta[2])/
                     (1+exp(theta[1]+3*theta[2])))
    accept<-accept+1
  }
  else
  {
    result[i, ]<-c(i,theta[1],theta[2],exp(theta[1]+3*theta[2])/
                     (1+exp(theta[1]+3*theta[2])))
  }
}
proc.time()-ptm
acceptratio<-accept/rep
acceptratio
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
#################################################################
# CODE1309 (cont)
### USE THE R-FUNCTION OPTIM TO COMPARE WITH THIS PROGRAMME
## FUNCTION TO COMPUTE THE NEGATIVE OF THE LOG-POSTERIOR
logpostoptim <- function(data,par)
{
  theta[1] <-par[1]
  theta[2] <- par[2]
  with(data=d,-sum(Y*( theta[1] + theta[2]*X)-log(1+exp(theta[1] + 
                                                          theta[2]*X))))
}
result1 <- optim(par=c(-3,1),logpostoptim,
                 data=d,hessian=TRUE,method= "BFGS")
result1$par
solve(result1$hessian)
#################################################################
# CODE1309 (cont)
# computing MODE with package modeest and with density estimation
#install.packages("modeest", .libPaths()[1])
y<- result[,2]
x<-result[,3]
library(modeest)
myDensity<-density(y)
modeEstmu <- mlv(y,method = "kernel", kernel = "gaussian")
modeEstmu
modeEstbeta <- mlv(x,method = "kernel", kernel = "gaussian")
modeEstbeta
################################################################
# CODE1309 (cont)
#POST-McMC ANALYSIS
# CODE FOR THE MC VARIANCE BASED ON GEYER
ns<-rep
# CHOOSE MU OR V AND PLACE IN VECTOR Y
y<- result[,2]
#y<-result[,3]
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
