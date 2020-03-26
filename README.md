# c2plasso: the structural varying-coefficient regression and the pliable lasso 
The R package `c2plasso` implements the structural varying-coefficient regression (svReg) for the model selection of a varying-coefficient model with structured main predictors or modifying variables. It also provides functions for performing the pliable lasso (plasso) by Tibshirani and Friedman (2019).

The reference for the pliable lasso can be found at:
* [A Pliable Lasso](https://doi.org/10.1080/10618600.2019.1648271) by Tibshirani and Friedman (2019).

## Installation

```
devtools::install_github("rakheon/c2plasso", force = TRUE)
```

## Example

```
# data generation
x=matrix(rnorm(100*5, 0, 1),100,5)
z1=matrix(rnorm(100*3, 0, 1),100,3)
z2=matrix(as.factor(sample(0:3, 100*2, prob=c(1/4,1/4,1/4,1/4), replace = TRUE)),100,2)
z2=as.data.frame(model.matrix(~., data=as.data.frame(z2))[,-1])
z=cbind(z1, z2)
z=as.matrix(z)
y=2*x[,1] - (2+2*z[,1])*x[,2] + (2+3*z[,4]+2*z[,5]-2*z[,6])*x[,3] + rnorm(100, 0, 1)

# fitting the plasso and the svReg for a sequence of tuning parameters
plasso_res = plasso(X = x, Z = z, Y = y, lambda_seq = c(1, 0.5), alpha = 0.5)
svReg_res = svReg(X = x, Z = z, Y = y, df_X = rep(1,5), df_Z = c(1,1,1,3,3), lambda_seq = c(1, 0.5), alpha = 0.5)

# Perform k-fold cross validation
cv.plasso_res = cv.plasso(X = x, Z = z, Y = y, lambda_seq = c(1, 0.5), alpha = 0.5)
cv.svReg_res = cv.svReg(X = x, Z = z, Y = y, df_X = rep(1,5), df_Z = c(1,1,1,3,3), lambda_seq = c(1, 0.5), alpha = 0.5)

# Correlated design for grouped main predictors
x=matrix(rnorm(100*5, 0, 1),100,5)
z1=matrix(rnorm(100*3, 0, 1),100,3)
z2=matrix(as.factor(sample(0:3, 100*2, prob=c(1/4,1/4,1/4,1/4), replace = TRUE)),100,2)
z2=as.data.frame(model.matrix(~., data=as.data.frame(z2))[,-1])
z=cbind(z1, z2)
z=as.matrix(z)
x[,3] = 2/3*x[,1] + 2/3*x[,2] + 1/3*rnorm(100, 0, 1)
y = x[,1] + x[,2] + (2+3*z[,4]+2*z[,5]-2*z[,6])*x[,4] + rnorm(100, 0, 1)

# fitting the svReg for a sequence of tuning parameters and perform k-fold cross validation
svReg_res = svReg(X = x, Z = z, Y = y, df_X = c(3,1,1), df_Z = c(1,1,1,3,3), lambda_seq = c(1, 0.5), alpha = 0.5)
cv.svReg_res = cv.svReg(X = x, Z = z, Y = y, df_X = c(3,1,1), df_Z = c(1,1,1,3,3), lambda_seq = c(1, 0.5), alpha = 0.5)
```
