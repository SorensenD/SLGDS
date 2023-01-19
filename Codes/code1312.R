# CODE1312
#SINGLE-SITE GIBBS - CORRELATED PROBIT MODEL
# DOES NOT USE THE SVD OF ZZ'
rm(list=ls()) # Clear the workspace
set.seed(7713)

require(graphics)
# GENERATE CORRELATED (FULL-SIBS) BINARY DATA  (THRESHOLD MODEL)
#install.packages("MCMCpack", .libPaths()[1])
#install.packages("mvtnorm", .libPaths()[1])
library(mvtnorm)
#library(MCMCpack)
# INITIALISE PARAMETERS
p0 <- 0.5
mu <- qnorm(p0)
iccfs<-0.15 #INTRACLASS CORRELATION FS
# VARIANCE BETWEEN FAMILIES: iccfs /(1- iccfs)
# PHENOTYPIC VARIANCE: 1/(1-iccfs)
nfs<-500 # NUMBER OF FULL-SIB FAMILIES
fs<-2 #FULL-SIB FAMILY SIZE
# SET DATA Y= 0
y<-matrix(data=0,nrow=fs*nfs,ncol=1)
# x IS COLUMN MATRIX WITH FAMILY ID (ID=1,.,nfs)
x<-matrix(data=0,nrow=fs*nfs,ncol=1)
# GENERATE NFS FULL-SIB EFFECTS f
f<-rnorm(nfs,mean=0,sd=sqrt(iccfs/(1-iccfs)))
# SET COUNTER c EQUAL TO ZERO 
# (c = nfs*fs IS EQUAL TO THE LENGTH OF BINARY DATA VECTOR y)
c<-0
##########################################################
####  GENERATE BINARY RECORDS Y
f<-rnorm(nfs,mean=0,sd=sqrt(iccfs/(1-iccfs)))
p <- pnorm(mu+f)
y <- rbinom(nfs*fs,1,rep(p,each=fs))
w <- rep(1:nfs,each=fs)
d<-data.frame(w,y)
family <- w
family <- as.factor(family)
Z<-model.matrix(~0+family)
# WITH INDEPENDENT FAMLIES Z'Z IS DIAGONAL
ztz<-rep(fs,nfs)
#CHOOSE LENGTH OF CHAIN rep
rep<-2000
result<-matrix(data=NA,nrow=rep,ncol=3)
# INITIALISE LIABILITY VECTOR u
u<-rep(0,fs*nfs)
#INITIALISE THE VECTOR OF FAMILIY EFFECTS fam
fam<-rep(0,nfs)
# INITIALISE BETWEEN FAMILY VARIANCE COMPONENT vf
vf<-0.2
# INITIALISE THE MEAN (HERE BETA, A SCALAR)
beta<-0
#START GIBBS LOOP
ptm<-proc.time()
for (i in 1:rep)
{
  zfam <- Z%*%fam
  # SAMPLE BETA
  betahat<-sum(u-zfam)/(fs*nfs)
  beta<-rnorm(1,mean=betahat,sd=sqrt(1/(fs*nfs)))
  # SAMPLE LIABILITIES u FROM A TN(beta+Zf,1)
  av<-beta+ zfam # MEAN OF UNTRUNCATED NORMAL(mu+Zf,1)
  prob<-pnorm(-av)
  interm<-( prob *(1-y)+(1- prob)*y)*runif(fs*nfs)+ prob*y
  interm[interm==1]<-0.999
  u<-qnorm(interm)+av
  # END SAMPLING LIABILITIES u
  # SAMPLE FAMILY EFFECTS fam
  varfam<-1/(fs+(1/vf))
  fammean<-varfam*(t(Z)%*%(u-beta))
  fam<-rnorm(nfs,mean=fammean, sd=sqrt(varfam))
  #SAMPLE vf
  #COMPUTE SCALE
  ftf<-sum(fam*fam)
  vf<-ftf/rchisq(1,nfs-2)
  result[i,]<-c(i,beta,vf)
}
proc.time()-ptm
# CONSTRUCT THE DERIVED PARAMETER "HERITABILITY" herit
herit <- 2*result[,3]/(result[,3]+1)
meanbeta <- mean(result[,2])
meanbeta
meanvf <- mean(result[,3])
meanvf
meanher <- mean(herit)
meanher
cibeta <- quantile(result[,2],c(0.025,0.975))
cibeta
civf <- quantile(result[,3],c(0.025,0.975))
civf
ciher <- quantile(herit,c(0.025,0.975))
ciher
