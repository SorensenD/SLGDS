# CODE0902
rm(list=ls()) # CLEAR WORKSPACE
set.seed(3033)
rep<-10000
se<-0.85
sp<-0.95
a_se<-4
b_se<-1.5
a_sp<-10
b_sp<-1.5

# THEORETICAL MODE OF BETA PRIORS
theoretmodese<-(a_se-1)/(a_se+b_se-2)
theoretmodesp<-(a_sp-1)/(a_sp+b_sp-2)
resultg<-matrix(data=NA,nrow=rep,ncol=4)

# NUMBER OF TESTS: n
n<-10000
# NUMBER OF POSITIVE TESTS: t
t<-900
for (i in 1:rep){
  psi<-rbeta(1,(t+1),(n-t+1))
  a<-1-sp
  b<-se+sp-1
  if((psi-a) > 0)
  {
    theta <- (psi-a)/b
  }
  else {
    theta <- 0
  }
  
  # DRAW SE SP FROM TRUNCATED BETA'S USING DEVROYE'S ALGORITHM
  se<-qbeta(runif(1,pbeta(0.8,a_se,b_se),pbeta(1,a_se,b_se)),
            a_se,b_se)
  sp<-qbeta(runif(1,pbeta(0.925,a_sp,b_sp),pbeta(1,a_sp,b_sp)),
            a_sp,b_sp)
  #se <- 0.85
  #sp <- 0.95
  # OR
  # DRAW SE AND SP FROM APPROPRIATE UNIFORM DISTRIBUTIONS
  #se<-runif(1,0.8,1)
  #sp<-runif(1,0.925,1)
  resultg[i,]<-c(i,theta,se,sp)
}
# BAYESIAN POSTERIOR MEAN OF THETA AND POSTERIOR INTERVAL
postmean <- mean(resultg[,2])
postinterval <- quantile(resultg[,2],c(0.025,0.975))
postmean
postinterval
