# CODE0303
# DATA BASED ON GENOMIC MODEL; OBTAIN THE SVD OF WW'(1/m)
rm(list=ls()) # CLEAR WORKSPACE
set.seed(1953)
nindiv<-10
nmark<-20
X<-matrix(nrow= nindiv,ncol= nmark,
          rbinom(n=nindiv*nmark,size=2,p=.5))
W <- matrix(data=NA,nrow= nindiv,ncol=nmark)
U <- matrix(data=NA,nrow= nindiv,ncol= nindiv)
G<-matrix(data=NA,nrow= nindiv,ncol= nindiv)
cm <- colMeans(X)
# CREATE MATRIX OF STANDARDISED MARKER GENOTYPE CODES
for (i in 1:nmark)
{
  W[,i] <-( X[,i]-cm[i]) / sd(X[,i])
}
# THIS IS MORE EFFICIENT THAN THE LOOP:
# W <- scale(X, center=TRUE, scale=TRUE) 
qr(X)$rank
qr(W)$rank

# GENOMIC RELATIONSHIP MATRIX G
G <- (1/nmark)*W%*%t(W)
# THIS IS MORE EFFICIENT THAN THE LINE ABOVE:
# G <- (1/nmark)*tcrossprod(W) 
# SVD OF G
EVD <- eigen(G)
names(EVD)
head(EVD$values[1:5])
U <- EVD$vector
val <- EVD$values
val[nindiv] <-0
D <- diag(val,nrow=nindiv)
# CHECK THAT G = UDU':
identical(G, U%*%D%*%t(U))
max(abs(G - U%*%D%*%t(U)))
