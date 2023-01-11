# CODE0306
#### EM algorithm ##########################
rm(list=ls())
set.seed(30371)
niter<-9
result <- matrix(data=NA,nrow=niter,ncol=3)

# DATA
n_A<-725
n_AB<-72
n_B<-258
n_0<-1073
# START VALUES FOR P_A and p_B
p_A <- 0.2
p_B <- 0.2
for (i in 1:niter){
  # E-step
  n_AA <- n_A * p_A^2/(p_A^2 + 2*p_A*(1-p_A-p_B))
  n_A0 <- n_A * (2*p_A*(1-p_A-p_B))/(p_A^2 + 2*p_A*(1-p_A-p_B))
  n_BB <- n_B * p_B^2/(p_B^2 + 2*p_B*(1-p_A-p_B))
  n_B0 <- n_B * (2*p_B*(1-p_A-p_B))/(p_B^2 + 2*p_B*(1-p_A-p_B))
  # M-step
  p_A <- (2*n_AA + n_AB + n_A0)/(2*(n_A + n_AB + n_B + n_0))
  p_B <- (2*n_BB + n_AB + n_B0)/(2*(n_A + n_AB + n_B + n_0))
  loglik <- (725*log(p_A*(2 - p_A - 2* p_B)) + 72*log(2*p_A*p_B) + 
               258*log(p_B*(2 - p_B - 2*p_A)) + 2*1073*log(1 - p_A - p_B))
  result[i,] <- c(p_A,p_B,loglik)
  
}
result
