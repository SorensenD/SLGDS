# CODE0502
rm(list=ls()) # CLEAR WORKSPACE
########### DRAWING ALL U's IN ONE GO
### REQUIRES TO GENERATE BINARY RECORDS y
##  USING PARAMETERISATION A OR B BELOW
# A
##### y = 0 -> TN(mean,1)(0,Infinity)
##### y = 1 -> TN(mean,1)(-Infinity,0)
########## OR #################
# B
##### y = 1 -> TN(mean,1)(0,Infinity)
##### y = 0 -> TN(mean,1)(-Infinity,0)
#####################
set.seed(237777)
nrow <- 10
ncol <- 5
mu <- 0
# GENERATE X MATRIX
X<-matrix(nrow= nrow,ncol= ncol,rbinom(nrow*ncol,size=2,p=.5))
# GENERATE VECTOR b
b <- rnorm(ncol,0.5,1)
xb<-X%*%b
# LOGIT MODEL
#p1<-exp(mu+xb)/(1+exp(mu+xb)) # IF B
#p1 <- 1 - exp(mu+xb)/(1+exp(mu+xb)) # IF A
# PROBIT MODEL
#p1 <- pnorm(mu+xb) # IF B
p1 <- 1 - pnorm(mu+xb) # IF A
y <- rbinom(nrow,1,p1)
mean <- mu+xb
sd <- 1

interm<-(1-y)*pnorm(0,mean=mean,sd=sd)+runif(length(y))*
  (pnorm(0,mean=mean,sd=sd)*(y)+(1-pnorm(0,mean=mean,sd=sd))*(1-y))

u <- qnorm(interm,mean=mean,sd=sd)
p1[1:5]
p1[6:nrow]
y[1:5]
y[6:nrow]
u[1:5]
u[6:nrow]

