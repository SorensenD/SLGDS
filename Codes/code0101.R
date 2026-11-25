# CODE0101
# READING SINGH ET AL. (2002) DATA
rm(list=ls()) # CLEAR WORKSPACE
# Lasso solutions using package glmnet
#install.packages("glmnet", .libPaths()[1])
#install.packages("sda")
library("sda")
library(glmnet)
data(singh2002)

X <- singh2002$x
y <- ifelse(singh2002$y == "cancer", 1, 0)
n <- nrow(X)

# DIVIDE DATA INTO TRAINING AND VALIDATING SETS

set.seed(3037)
train <- sample(1:n, n/2)
test <- setdiff(1:n, train)

y.train <- y[train]
y.test <- y[test]

# CHOOSE THE FIRST 15 COVARIATES

p <- 15
X1 <- X[, 1:p]

# FIT LOGISTIC REGRESSION TO THE TRAINING DATA

glm.fit <- glm(y[train] ~ X1[train, ],
               family = binomial(link = "logit"))

b <- coef(glm.fit)

# COMPUTE ESTIMATED LINEAR PREDICTORS

etaT <- b[1] + X1[train, ] %*% b[-1]
etaV <- b[1] + X1[test, ] %*% b[-1]

# COMPUTE ESTIMATED PROBABILITIES

probT <- plogis(etaT)
probV <- plogis(etaV)

# COMPUTE PREDICTIONS ON THE BINARY SCALE

predT <- ifelse(probT >= 0.5, 1, 0)
predV <- ifelse(probV >= 0.5, 1, 0)

# COMPUTE MISCLASSIFICATION ERRORS

MSEt <- mean((y.train - predT)^2)
MSEv <- mean((y.test - predV)^2)

c(p = p, MSEt = MSEt, MSEv = MSEv)

