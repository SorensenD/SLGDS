# CODE0307
# EM FOR TRUNCATED DATA; ESTIMATE MEAN OF UNTRUNCATED
# GENERATE Y ~ N(MEAN,VAR)
# TRUNCATE AT T SO THAT Z = Y > T ARE OBSERVED
# Y < T ARE MISSING (KNOWN INFORMATION)
rm(list=ls()) # CLEAR WORKSPACE
set.seed(12371)
nindiv<-50000
mean <- 10
var <- 3
T <- mean + 1.5*sqrt(var)
# CREATE DATA
y <- rnorm(nindiv,mean,sqrt(var))
z <- y[y>T]
length(z)
m <- length(y)-length(z)
mean(y)
mean(z)
var(z)
# START VALUES FOR MEAN (mu) AND VARIANCE (sigmasq)
mu <- 0
sigmasq <- 2
sigma <- sqrt(sigmasq)
iter <- 1000
res <- matrix(data=NA, nrow=iter,ncol=2)
###  EM LOOP  ##########################
for (i in 1:iter){
  T_star <- (T-mu)/sigma
  expymiss <- mu - (sigma * dnorm(T_star))/(pnorm((T_star)))
  mu <- (sum(z)+m*expymiss)/length(y)
  e <- z-mu
  sigmasq <-(m*sigmasq*(1-T_star*dnorm(T_star)/pnorm(T_star))+
               sum(e*e))/length(y)
  sigma <- sqrt(sigmasq)
  res[i,] <- c(mu,sigmasq)
}
tail(res)
emmu <- res[iter,1]
emmu
emsigmasq <- res[iter,2]
emsigmasq
emsel <- mu + sigma*dnorm(T_star)/(1-pnorm(T_star))
emsel
i <- dnorm(T_star)/(1-pnorm(T_star))
varsel <- sigmasq*(1-i*(i-T_star))
varsel
