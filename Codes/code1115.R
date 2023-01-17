# CODE1115
###      GENETIC MODELS
rm(list=ls()) # CLEAR WORKSPACE
set.seed(303371)
nindiv<-1000000
nqtl <- 2
nintqtl <- 1
ba <- rep(1,nqtl)
bi <- rep(1,nintqtl)
mu_y<-0
# INTERACTION GENOTYPIC MATRIX
Xi<-matrix(data=NA,nrow= nindiv,ncol= nintqtl) 
Xq<-matrix(nrow= nindiv,ncol= nqtl,
           rbinom(nindiv*nqtl,size=2,p=.5)-1) # LINEAR TERMS 
nr <- ncol(Xq)
i1 <- combn(nr,2)
i2 <- sample(ncol(i1),nintqtl,replace=FALSE)
i3 <- as.matrix(i1[,i2])
# CONSTRUCT INTERACTION GENOTYPE
for (i in 1:nintqtl){
  Xi[,i] <- Xq[,i3[1,i]]*Xq[,i3[2,i]]+1
}
gi <-Xi%*%bi # INTERACTION GENETIC VALUES
ga <- Xq%*%ba # ADDITIVE GENETIC VALUES
g <- ga + gi # TOTAL GENETIC VALUES
va <- var(ga)
vi <- var(gi)
vg <- var(g)
y_i <- gi+rnorm(nindiv,0,sqrt(5-vi))
y_ai <- ga + gi+rnorm(nindiv,0,sqrt(5-va-vi))
vy <- var(y_i)
cor(Xq)
vy
vi
va
vg
# FIT LINEAR REGRESSION OF INTERACTION GENETIC VALUES ON QTL LOCI
fa <- lm(gi ~ Xq)
# CONFIRM THAT THE MODEL DOES NOT CAPTURE ANY 
# (ADDITIVE GENETIC) VARIATION
summary(fa)

