# CODE0101
# EIGEN DECOMPOSITION
rm(list=ls()) # CLEAR WORKSPACE
set.seed(30371171)
# SIMULATE Z FROM AN UNSTRUCTURED POPULATION
Z <- matrix(nrow= 500,ncol= 1300,rbinom(500*1300,size=2,p=.5))
Gz <- tcrossprod(scale(Z)) # CENTRING AND SCALING
EVD <- eigen(Gz)
U <- EVD$vector
tU<-t(U)
val <- EVD$values
qr(Gz)$rank
plot(U[,1],U[,2],xlab='U1',ylab='U2')

# SIMULATE Z FROM TWO POPULATIONS WITH DIFFERENT GENE FREQUENCIES
Z1 <- matrix(nrow= 250,ncol= 1300,rbinom(250*1300,size=2,p=.5))
Z2 <- matrix(nrow= 250,ncol= 1300,rbinom(250*1300,size=2,p=.3))
Z <- rbind(Z1,Z2)
Gz <- tcrossprod(scale(Z))
EVD <- eigen(Gz)
U <- EVD$vector
tU<-t(U)
val <- EVD$values
qr(Gz)$rank
plot(U[,1],U[,2],xlab='U1',ylab='U2')