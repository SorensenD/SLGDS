# CODE0103
rm(list=ls()) # CLEAR WORKSPACE
set.seed(30331)
#install.packages("tree")
library(sda)
library(tree)
# library(glmnet)
data(singh2002)
d <- data.frame(singh2002$x)
d$y <- singh2002$y
nrep <- 1 # NUMBER OF REPLICATES
res <- matrix(data=NA,nrow=nrep,ncol=3)
ptm<-proc.time()
for ( i in 1:nrep ) {
  cat(i,"\n",sep="")
  train <- c(sample(1:50,25),sample(51:102,26))
  # FIT THE TREE TO THE TRAINING DATA    
  trees <- tree(y ~ . ,  data=d[train,])
  # FIT FUNCTION PREDICT TO THE TRAINING AND VALIDATING DATA    
  predtreev <- predict(trees,d[-train,],type="class")
  predtreet <- predict(trees,d[train,],type="class")
  # CALCULATE 1-CLASSIFICATION ERROR IN TRAINIMG AND VALIDATING DATA    
  predv <- sum(predtreev==d$y[-train])/length(d$y[-train])
  predt <- sum(predtreet==d$y[train])/length(d$y[train])
  # RECORD TRAINING / VALIDATING CLASSIFICATION ERROR AND 
  # NUMBER OF COVARIATES IN TREE 
  res[i,]<-c((1-predt),(1-predv),length(summary(trees)$used))
}
proc.time()-ptm
res
tab <- table(predtreev,d$y[-train])
tab
# CHECK CLASSIFICATION ERROR
(tab[1,2]+tab[2,1])/(length(d$y[-train]))
#summary(res)
plot(trees)
text(trees,pretty=0)
