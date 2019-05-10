# c2plasso: continuous-categorical pliable lasso
The R package c2plasso implements continuous-categorical pliable lasso for the data with modifying variable of the form of continuous variables, categorical variable or both. It also provides functions for performing the originial pliable lasso.

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

# fitting plasso and c2plasso for a sequence of tuning parameters
plasso(X = x, Z = z, Y = y, lambda_seq = c(1, 0.5), alpha = 0.5)
c2plasso(X = x, Z = z, Y = y, df_Z = c(1,1,1,3,3), lambda_seq = c(1, 0.5), alpha = 0.5)

# Perform k-fold cross validation
cv.plasso(X = x, Z = z, Y = y, lambda_seq = c(1, 0.5), alpha = 0.5)
cv.c2plasso(X = x, Z = z, Y = y, df_Z = c(1,1,1,3,3), lambda_seq = c(1, 0.5), alpha = 0.5)
```
