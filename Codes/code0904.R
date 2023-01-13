# CODE0904
# PREDICTION OF A GENETIC DISEASE - NUMERICAL INTEGRALS
rm(list=ls()) # CLEAR WORKSPACE
library(mvtnorm)
# EXAMPLE: FATHER-MOTHER-CHILD
p<-0.02 # INCIDENCE IN THE POPULATION
# p <- 0.002
mean<-c(0,0,0) # MEAN OF THE THREE LIABILITIES
her <- 1/3 # heritability
var <- 1.5 # variance of liability
cov <- 0.5*her*var # covariance single parent-child
t<-qnorm((1-p),mean=0,sd=sqrt(1.5)) # THRESHOLD
# # VAR-COV MATRIX OF LIABILITY:
sigma <- matrix(c(var,cov,cov,cov,var,0,cov,0,var),3,3) 
sigma
# CASE 1: FATHER AFFECTED, MOTHER AFFECTED
den<-pmvnorm(lower=c(-Inf,t,t),upper=c(Inf,Inf,Inf),
             mean=mean,sigma=sigma)
num<-pmvnorm(lower=c(t,t,t),upper=c(Inf,Inf,Inf),
             mean=mean,sigma=sigma)
proboffsprcase1 <- num/den
proboffsprcase1
# CASE 2: FATHER AFFECTED, MOTHER UNAFFECTED
den<-pmvnorm(lower=c(-Inf,t,-Inf),upper=c(Inf,Inf,t),
             mean=mean,sigma=sigma)
num<-pmvnorm(lower=c(t,t,-Inf),upper=c(Inf,Inf,t),
             mean=mean,sigma=sigma)
proboffsprcase2 <- num/den
proboffsprcase2
# CASE 3: FATHER UNAFFECTED, MOTHER UNAFFECTED
den<-pmvnorm(lower=c(-Inf,-Inf,-Inf),upper=c(Inf,t,t),
             mean=mean,sigma=sigma)
num<-pmvnorm(lower=c(t,-Inf,-Inf),upper=c(Inf,t,t),
             mean=mean,sigma=sigma)
proboffsprcase3 <- num/den
proboffsprcase3
