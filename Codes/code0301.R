##################################################
# CODE0301
rm(list=ls())
set.seed(30371)
fd<-matrix(data=NA,nrow=2,ncol=1)
sd<-matrix(data=NA,nrow=2,ncol=2)
freq<-matrix(data=NA,nrow=2,ncol=1)
niter<-20
# DATA
n_A<-725
n_AB<-72
n_B<-258
n_O<-1073
# INITIALIZE GENE FREQ
freq[1,1]<-0.3
freq[2,1]<-0.2
p_A<-freq[1,1]
p_B<-freq[2,1]
# ITERATION LOOP
for (i in 1:niter){
  fd[1,1] <- n_AB/p_A + n_A*(2-2*p_A-2*p_B)/(p_A*(2-p_A-2*p_B)) -
    2*n_B/(2-2*p_A-p_B) - 2*n_O/(1-p_A-p_B)
  fd[2,1] <- n_AB/p_B + n_B*(2-2*p_A-2*p_B)/(p_B*(2-2*p_A-p_B)) -
    2*n_A/(2-p_A-2*p_B) - 2*n_O/(1-p_A-p_B)
  s11a <- -n_AB/((p_A)^2)
  s11b <- (n_A*(2-2*p_A-2*p_B))/((p_A*(2-p_A-2*p_B)^2))
  s11c <- - (2*n_A)/((p_A*(2-p_A-2*p_B)))
  s11d <-  - (n_A*(2-2*p_A-2*p_B))/((p_A^2*(2-p_A-2*p_B)))
  s11e <- - (4*n_B)/((2-2*p_A-p_B)^2)
  s11f <- - (2*n_O)/((1-p_A-p_B)^2)
  sd[1,1] <- s11a + s11b + s11c + s11d + s11e + s11f
  s22a <- - n_AB/((p_B)^2)
  s22b <-  n_B*(2-2*p_A-2*p_B)/(((2-2*p_A-p_B)^2) * p_B)
  s22c <- - (2*n_B)/(p_B*(2-2*p_A-p_B))
  s22d <- - (n_B*(2-2*p_A-2*p_B))/((p_B^2*(2-2*p_A-p_B)))
  s22e <- -(4*n_A)/((2-p_A-2*p_B)^2)
  s22f <- -(2*n_O)/((1-p_A-p_B)^2)
  sd[2,2] <- s22a + s22b + s22c + s22d + s22e + s22f
  sd[1,2] <- -2*n_A/((2-p_A-2*p_B)^2) - 2*n_O/((1-p_A-p_B)^2) - 
    2*n_B/((2-2*p_A-p_B)^2)
  sd[2,1] <- sd[1,2]
  freq<-freq-solve(sd)%*%fd
  p_A<-freq[1,1]
  p_B<-freq[2,1]
}
# ML ESTIMATES ARE
freq
# OBSERVED INFORMATION IS
-sd
# ASYMPTOTIC VAR-COVAR MATRIX IS
-solve(sd)
