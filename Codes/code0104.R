# CODE0104
#install.packages("randomForest")
rm(list=ls()) # CLEAR WORKSPACE
library(sda)
library(randomForest)
data(singh2002)
d <- data.frame(X=singh2002$x)
d$y <- singh2002$y
n0 <- sum(d$y=="healthy")
n1 <- sum(d$y=="cancer")
set.seed(3037)
p <- .5
nrep <- 1
# nrep <- 200 (to construct figure 1.11 in chapter 1)
#mtry <- c(2,3,4,  5,  8, 10, 13, 17, 20, 30, 40, 50, 60, 70, 80)
mtry <- c(5,20,50,80,120)
#mtry <- c(20,80)
#mtry <- c(80)
sumd <- data.frame()
res <- rep(NA,nrep)
ptm<-proc.time()
for ( m in mtry) {
  cat("mtry ",m,"\n",sep="")
  for ( rep in 1:nrep ) {
    cat("Replicate ",rep,"\n",sep="")
    train <- c(sample( 1:n0,floor(p*n0) ),
               sample( (n0+1):(n0+n1),floor(p*n1) ))
    rf.singh =randomForest(y ~.,
                           data=d,
                           subset =train,
                           mtry=m,
                           importance =TRUE)
    predict <- predict(rf.singh,d[-train,])
    observed <- d$y[-train]
    t <- table(observed,predict)
    print(t)
    res[rep] <- (t[1,1]+t[2,2])/sum(t)
  }
  sumd <- rbind(sumd,c(m,min(res),mean(res),median(res),
                       max(res),var(res)))
}
proc.time()-ptm
names(sumd) <- c("mtry","min","mean","median","max","var")

with(sumd,plot(mtry,mean,type="l",col="red",ylim=c(min(min),1),
               ylab="1-Validating MSE",
               xlab="Number of Predictors/Split"))
with(sumd,lines(mtry,min,lty=2,col="blue"))
with(sumd,lines(mtry,max,lty=2,col="blue"))
