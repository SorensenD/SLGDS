########################################################
# code0302
rm(list=ls()) # CLEAR WORKSPACE
set.seed(37)
fr<-function(par){
  p_A <- par[1]
  p_B <- par[2]
  -(725*log(p_A*(2 - p_A - 2* p_B)) + 72*log(2*p_A*p_B) + 
      258*log(p_B*(2 - p_B - 2*p_A)) + 2*1073*log(1 - p_A - p_B))}
result <- optim(par=c(0.3,0.2),fr,hessian=TRUE)
result$par
solve(result$hessian)
-result$value