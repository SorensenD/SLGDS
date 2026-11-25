# CODE0102
# READING SINGH ET AL. (2002) DATA
rm(list = ls())
suppressPackageStartupMessages(library(sda))
suppressPackageStartupMessages(library(glmnet))
data(singh2002)
X <- singh2002$x
y <- ifelse(singh2002$y == "cancer", 1, 0)
n <- nrow(X)
# DIVIDE DATA INTO TRAINING AND VALIDATING SETS
set.seed(3037)
train <- sample(1:n, n/2)
# NUMBER OF TRAINING/VALIDATING SPLITS
repl <- 1
#repl <- 50
result <- matrix(NA, nrow = repl, ncol = 4)
colnames(result) <- c("MSEt", "MSEv",
                      "lambda", "n.selected")
set.seed(3037)
for (i in 1:repl) {
  if (i > 1) {
    train <- sample(1:n, n/2)
  }
  test <- setdiff(1:n, train)
  y.train <- y[train]
  y.test  <- y[test]
  
  # STEP 1: CROSS-VALIDATION USING THE TRAINING DATA
  # FIND THE VALUE OF LAMBDA THAT MINIMISES
  # THE CROSS-VALIDATED MISCLASSIFICATION ERROR
  # alpha = 1: LASSO
  # alpha = 0: RIDGE REGRESSION
  
  cv.out <- cv.glmnet(X[train, ], y.train,
                      alpha = 1,
                      family = "binomial",
                      type.measure = "class")
  bestlam <- cv.out$lambda.min
  
  # STEP 2: EXTRACT ESTIMATED REGRESSION COEFFICIENTS
  # ASSOCIATED WITH THE SELECTED VALUE OF LAMBDA
  b <- as.numeric(coef(cv.out, s = "lambda.min"))
  # NUMBER OF NON-ZERO REGRESSION COEFFICIENTS
  # EXCLUDING THE INTERCEPT
  n.selected <- sum(b[-1] != 0)
  # COMPUTE ESTIMATED LINEAR PREDICTORS
  # IN THE TRAINING AND VALIDATING DATA
  etaT <- b[1] + X[train, ] %*% b[-1]
  etaV <- b[1] + X[test, ] %*% b[-1]
  # COMPUTE ESTIMATED PROBABILITIES
  probT <- plogis(etaT)
  probV <- plogis(etaV)
  # COMPUTE PREDICTIONS ON THE BINARY SCALE
  predT <- ifelse(probT >= 0.5, 1, 0)
  predV <- ifelse(probV >= 0.5, 1, 0)
  # COMPUTE TRAINING AND VALIDATING MSE
  # (MISCLASSIFICATION ERRORS)
  MSEt <- mean((y.train - predT)^2)
  MSEv <- mean((y.test - predV)^2)
  # STORE RESULTS
  result[i, ] <- c(MSEt, MSEv,
                   bestlam, n.selected)
}
out <- as.data.frame(round(result, 4))
out$n.selected <- as.integer(out$n.selected)
print(out, row.names = FALSE)
# SUMMARY OF RESULTS OVER REPLICATES
if (repl > 1) {
  # VALIDATING MEAN SQUARED ERROR
  mean.MSEv <- mean(result[, "MSEv"])
  range.MSEv <- range(result[, "MSEv"])
  # NUMBER OF SELECTED COVARIATES
  median.selected <- median(result[, "n.selected"])
  range.selected <- range(result[, "n.selected"])
  print(mean.MSEv)
  print(range.MSEv)
  print(median.selected)
  print(range.selected)
}
# ALTERNATIVE IMPLEMENTATION USING predict()
# WITHOUT EXPLICITLY EXTRACTING THE ESTIMATED
# REGRESSION COEFFICIENTS
predT <- predict(cv.out,
                 newx = X[train, ],
                 s = "lambda.min",
                 type = "class")
predV <- predict(cv.out,
                 newx = X[test, ],
                 s = "lambda.min",
                 type = "class")
MSEt <- mean((y.train - as.numeric(predT))^2)
MSEv <- mean((y.test - as.numeric(predV))^2)
c(MSEt = MSEt, MSEv = MSEv)

