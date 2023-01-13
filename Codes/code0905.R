# CODE0905
# rm(list=ls()) # CLEAR WORKSPACE
#        Pearson-Aitken formula
# Input
mean<-0 # marginal mean of liability
var<-1.5 # marginal variance of liability
varg<-0.5 # additive genetic variance
p<-0.02 # Incidence in population
t<-qnorm((1-p),mean=0,sd=sqrt(1.5)) # THRESHOLD
varfm<-var*diag(2) # VAR-COV PARENTS
covop<-matrix(0.5*varg,nrow=1,ncol=2) # COV OFFS-PARENTS
# CASE 1: FATHER AFFECTED, MOTHER AFFECTED
# FATHER: 
alfa<-(t-mean)/sqrt(var) # LOWER TRUNCATION
p_alfa<-dnorm(alfa) # PDF AT LOWER TRUNCATION POINT
cum_alfa<-pnorm(alfa) # CUM. DIST. FUNCTION AT LOWER TRUNCATION
intsel<-p_alfa/(1-cum_alfa) # "INTENSITY OF SELECTION"
# MEAN AND VARIANCE OF (SELECTED) FATHER (MOTHER BELOW)
minfather<-mean+sqrt(var)*intsel
varfather<-var*(1-intsel*(intsel-alfa))
# MOTHER
minmother<-minfather
varmother<-varfather
# MEAN AND VARIANCE OF SELECTED FATHER AND MOTHER
expcase1 <- matrix(c(minfather,minmother),2,1)
varcase1 <- matrix(c(varfather,0,0,varmother),2,2)
# CONDITIONAL MEAN AND VARIANCE OF LIABILITY OF OFFSPRING
condmin<-mean+covop%*%solve(varfm)%*%(expcase1)
int1<-solve(varfm)-(solve(varfm)%*%varcase1%*%solve(varfm))
condvar<-var-(covop%*%(int1)%*%t(covop))
proboffPAcase1<-1-pnorm(t,mean=condmin,sd=sqrt(condvar))
proboffPAcase1
# *********************************************************
# CASE 2: FATHER AFFECTED, MOTHER UNAFFECTED
# FATHER: AS IN CASE 1
beta<-(t-mean)/sqrt(var) # upper TRUNCATION
p_beta<-dnorm(beta)
cum_beta<-pnorm(beta)
# MOTHER MEAN AND VARIANCE
minmother<-mean-sqrt(var)*(p_beta/cum_beta)
varmother<-var*(1-(p_beta/cum_beta)*((p_beta/cum_beta)+beta))
expcase2 <- matrix(c(minfather,minmother),2,1)
varcase2 <- matrix(c(varfather,0,0,varmother),2,2)
# CONDITIONAL MEAN AND VARIANCE OF LIABILITY OF OFFSPRING
condmin<-mean+covop%*%solve(varfm)%*%(expcase2)
int2<-solve(varfm)-(solve(varfm)%*%varcase2%*%solve(varfm))
condvar<-var-(covop%*%(int2)%*%t(covop))
proboffPAcase2<-1-pnorm(t,mean=condmin,sd=sqrt(condvar))
proboffPAcase2
# ********************************************************
# CASE 3: FATHER UNAFFECTED, MOTHER UNAFFECTED
# FATHER: AS MOTHER IN CASE 2
expcase3 <- matrix(c(minmother,minmother),2,1)
varcase3 <- matrix(c(varmother,0,0,varmother),2,2)

# CONDITIONAL MEAN AND VARIANCE OF LIABILITY OF OFFSPRING
condmin<-mean+covop%*%solve(varfm)%*%(expcase3)
int3<-solve(varfm)-(solve(varfm)%*%varcase3%*%solve(varfm))
condvar<-var-(covop%*%(int3)%*%t(covop))
proboffPAcase3<-1-pnorm(t,mean=condmin,sd=sqrt(condvar))
proboffPAcase3
