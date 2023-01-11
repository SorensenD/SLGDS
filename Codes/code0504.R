# CODE0504
#CPO EXAMPLE
rm(list=ls()) # Clear the workspace
set.seed(123771)
ptm<-proc.time()
require(graphics)
# GENERATE CORRELATED (FULL-SIBS DATA 
#install.packages("MCMCpack", .libPaths()[1])
#install.packages("mvtnorm", .libPaths()[1])
library(MCMCpack)
# INITIALISE PARAMETERS
mus<-10 # MEAN
vfs<-10 #VARIANCE BETWEEN FULL-SIBS
# RESIDUAL VARIANCE
ves<-50
nf<-400 # NUMBER OF FULL-SIB FAMILIES
n<-3 #FULL-SIB FAMILY SIZE
N<-nf*n
y<-matrix(data=0,nrow=nf*n,ncol=1)
# z IS COLUMN MATRIX WITH FAMILY ID (ID=1,.,nfs)
z<-matrix(data=0,nrow=nf*n,ncol=1)
# GENERATE nf FULL-SIB EFFECTS f
fs<-rnorm(nf,mean=0,sd=sqrt(vfs))
# GENERATE nf*n RESIDUAL EFFECTS f
es<-rnorm(nf*n,mean=0,sd=sqrt(ves))
#c<-0
#for (i in 1:nf)
#{
#  for (j in 1:n) # GENERATE FS FULL SIBS RECORDS
#  {
#    c<-c+1
#    y[c]<-mus+fs[i]+es[c]
#    z[c]<-i
#  }
#}
# GENERATE FULL SIBS (CAN CHOOSE MORE TRANSPARENT LOOP ABOVE)
z <- rep(1:nf,each=n)
y <- mus+fs[z]+es
d<-data.frame(y,z)
# GENERATE INCIDENCE MATRIX Z
family<-z
family <- as.factor(family)
Z<-model.matrix(~0+family)
# WITH INDEPENDENT FAMLIES Z'Z IS DIAGONAL
ztz<-rep(n,nf)
#END OF GENERATION OF DATA Y
#CHOOSE LENGTH OF CHAIN 
rep<-1000
resultt<-matrix(data=NA,nrow=rep,ncol=4)
resultw<-matrix(data=NA,nrow=rep,ncol=3)
#INITIALISE THE VECTOR OF FAMILIY EFFECTS f
f<-rep(0,nf)
# INITIALISE BETWEEN FAMILY VARIANCE COMPONENT vf
vf<-5
# INITILISE RESIDUAL VARIANCE
ve<-5
# INITIALISE k
k<-ve/vf
# INITIALISE THE MEAN
mu<-0
sumpyinvt<-0
#START GIBBS LOOP TRUE MODEL
for (i in 1:rep)
{
  # SAMPLE mu
  meanmu<-sum(y-Z%*%f)/(nf*n)
  mu<-rnorm(1,mean=meanmu,sd=sqrt(ve/(nf*n)))
  # SAMPLE FAMILY EFFECTS f
  varf<-(k+n)^(-1)
  fmean<- varf*(t(Z)%*%(y-mu))
  f<-rnorm(nf,mean=fmean, sd=sqrt(varf*ve))
  #SAMPLE vf
  #COMPUTE SCALE
  ftf<-sum(f*f)
  vfx<-ftf/rchisq(1,nf-2)
  vf<-as.numeric(vfx)
  # SAMPLE ve
  # COMPUTE SCALE
  e<-(y-mu-Z%*%f)
  ete<-t(e)%*%e
  vex<-ete/rchisq(1,N-2)
  ve<-as.numeric(vex)
  k<-ve/vf
  resultt[i,]<-c(i,mu,vf,ve)
  #  print(resultt[i,])
  # COMPUTE CPOs FOR TRUE MODEL
  pyinvt<-1/(dnorm(y,mean=mu+Z%*%f,sd=sqrt(ve)))
  sumpyinvt<-sumpyinvt+pyinvt
}
phatyt<-rep*(sumpyinvt)^(-1)
logcpot<-sum(log(phatyt))
proc.time()-ptm
#START GIBBS LOOP WRONG MODEL
vew<-20
sumpyinvw<-0
for (i in 1:rep)
{
  # SAMPLE muw
  meanmuw<-sum(y)/N
  muw<- rnorm(1,mean=meanmuw,sd=sqrt(vew/(N)))
  # SAMPLE vew
  vew<-sum((y-muw)*(y-muw))/rchisq(1,N-2)
  resultw[i,]<-c(i,muw,vew)
  #  print(resultw[i,])
  # COMPUTE CPOs FOR WRONG MODEL
  pyinvw<-1/(dnorm(y,mean=mu,sd=sqrt(vew)))
  sumpyinvw<-sumpyinvw+pyinvw
}
phatyw<- rep*(sumpyinvw)^(-1)
logcpow<-sum(log(phatyw))
# PRINT OUTPUT
# LOG CPO TRUE MODEL
logcpot
# 95% POSTERIOR INTERVALS FOR THE MEAN AND THE TWO VARIANCES
quantile(resultt[,2],c(0.025,0.975))
quantile(resultt[,3],c(0.025,0.975))
quantile(resultt[,4],c(0.025,0.975))

# LOG CPO FOR WRONG MODEL
logcpow
# 95% POSTERIOR INTERVALS FOR THE MEAN AND VARIANCE
quantile(resultw[,2],c(0.025,0.975))
quantile(resultw[,3],c(0.025,0.975))
