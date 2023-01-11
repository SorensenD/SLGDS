# CODE0503
# EM FOR TRUNCATED DATA; ESTIMATE MEAN OF UNTRUNCATED
# GENERATE Y ~ N(MEAN,VAR)
# TRUNCATE AT T SO THAT Z = Y > T ARE OBSERVED
# Y < T ARE MISSING (KNOWN INFORMATION)
rm(list=ls()) # CLEAR WORKSPACE
set.seed(12371)
nindiv<-50000
mean <- 10
var <- 3
T <- mean + 1.5*sqrt(var) # ASSUMED KNOWN
# CREATE COMPLETE DATA
y <- rnorm(nindiv,mean,sqrt(var))
# TRUNCATE: OBSERVED DATA
z <- y[y>T]
#length(z)
m <- length(y)-length(z)
#mean(y)
#mean(z)
#var(z)
#####################  McMC ######################
nrep <- 1000
resmc <- matrix(data=NA,nrow=nrep,ncol=2)
w <- rep(0,m)
# START VALUES FOR MEAN (mu) AND VARIANCE (sigmasq)
mu <- 0
sigmasq <- 2
sigma <- sqrt(sigmasq)
ptm <- proc.time()
for (j in 1:nrep){
  #  print(j)
  T_star <- (T-mu)/sigma
  std <- sqrt(var)
  # sample m missing records in one go (left from threshold T)
  w <- mu + std*qnorm(runif(m)*pnorm(T_star))
  # sample the variance
  scale <- sum((w-mu)^2) + sum((z-mu)^2)
  sigmasq <- scale/rchisq(1,length(y)-2)
  sigma <- sqrt(sigmasq)
  # sample the mean
  xbar <- (sum(w)+sum(z))/(length(w)+length(z))
  disp <- sigmasq/(length(w)+length(z))
  mu <- rnorm(1,xbar,sqrt(disp))
  resmc[j,] <- c(mu,sigmasq)
}
proc.time()-ptm
postmean <- mean(resmc[100:nrep,1])
postvar <- mean(resmc[100:nrep,2])
postmean
postvar
# 95% POSTERIOR INTERVAL FOR THE MEAN
pimean <- quantile(resmc[100:nrep,1],c(0.025,0.975))
# 95% POSTERIOR INTERVALFOR THE VARIANCE
pivar <- quantile(resmc[100:nrep,2],c(0.025,0.975))
pimean
pivar

