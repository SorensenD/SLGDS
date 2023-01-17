# CODE1114
# GENERATING AN EPISTATIC MODEL AND FIT RANDOM FOREST, 
# LASSO, BAYESIAN MIXTURE
rm(list=ls()) # CLEAR WORKSPACE
set.seed(303371)
nindiv<-10000
nmark <- 1000
nqtl <- 10
nintqtl <- 10
mu_y<-0
Xq<-matrix(data=NA,nrow= nindiv,ncol= nqtl)
########  GENETIC MARKERS Xm  #########################
Xm<-matrix(nrow= nindiv,ncol= nmark,
           rbinom(nindiv*nmark,size=2,p=.5)-1)
# from the nmark markers, choose nqtl as QTL:
IDq<-sample(1:nmark,nqtl,replace=F)  
Xq <- Xm[,IDq] # QTL GENOTYPIC MATRIX
# INTERACTION GENOTYPIC MATRIX:
Xi<-matrix(data=NA,nrow= nindiv,ncol= nintqtl) 
b <- rep(0,nmark+nintqtl)
nr <- ncol(Xq)
i1 <- combn(nr,2)
i2 <- sample(ncol(i1),nintqtl,replace=FALSE)
i3 <- as.matrix(i1[,i2])
for (i in 1:nintqtl){
  Xi[,i] <- Xq[,i3[1,i]]*Xq[,i3[2,i]]+1
}
# GENERATE GENOTYPIC VALUES g
b[IDq] <- 0.5
# #######################
# BELOW: ZERO OUT length(idzero) ADDITIVE EFFECTS
idzero <- sample(IDq,floor(0.9*nqtl),replace=FALSE)
b[idzero] <- 0
# #########################
lb <- nmark+1
ub <- nmark+nintqtl
b[lb:ub] <- 1.0
gi <-Xi%*%b[lb:ub] 
ga <- Xq%*%b[IDq]
g <- ga + gi
va <- var(ga)
vi <- var(gi)
vg <- var(g)
y <- ga+gi+rnorm(nindiv,0,sqrt(5-va-vi))
vy <- var(y)
her_a <- va/vy
her_i <- vi/vy
V <- vy*(1-her_a-her_i) # CONDITIONAL VARIANCE
cov(ga,gi)
va
vi
vy
