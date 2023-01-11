# CODE0501
# ABO BLOOD GROUP DATA - GIBBS SAMPLING
rm(list=ls()) # Clear the workspace
set.seed(1237)
#install.packages("MCMCpack", .libPaths()[1])
# to access function rdirichlet
library(MCMCpack)
#CHOOSE LENGTH OF GIBBS CHAIN rep
rep<-3000
result<-matrix(data=NA,nrow=rep,ncol=4)
# INITIALISE PARAMETERS
p_A<-0.33
p_B<-0.33
p_0<-1-p_A-p_B
alfa<-2
# DATA
n_A<-725
n_AB<-72
n_B<-258
n_0<-1073
#START WITH THE GIBBS LOOP
for (i in 1:rep)
{
  # SAMPLE n_AA AND n_BB
  n_AA<-rbinom(1,n_A,p_A^2/(p_A^2+2*p_A*p_0))
  n_BB<-rbinom(1,n_B,p_B^2/(p_B^2+2*p_B*p_0))
  n_A0<-n_A-n_AA
  n_B0<-n_B-n_BB
  # SAMPLE p_A,p_B,p_0
  a<-2*n_AA+n_A0+n_AB+alfa
  b<-2*n_BB+n_B0+n_AB+alfa
  c<-2*n_0+n_A0+n_B0+alfa
  draws<- rdirichlet(1,c(a,b,c))
  p_A<-draws[1,1]
  p_B<-draws[1,2]
  p_0<-draws[1,3]
  result[i, ]<-c(i,p_A,p_B,p_0)
}
# END OF GIBBS LOOP
meanp_A<-mean(result[,2])
meanp_A
varp_A<-var(result[,2])
cip_A<- quantile(result[,2],c(0.025,0.975))
cip_A
meanp_B<-mean(result[,3])
meanp_B
varp_B<-var(result[,3])
cip_B<- quantile(result[,3],c(0.025,0.975))
cip_B
covp_Ap_B<-cov(result[,2],result[,3])