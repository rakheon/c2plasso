plasso1 <- function(X, Z, Y, lambda = 0.5, alpha = 0.5, tt = 0.1, beta = NULL, theta = NULL, zlinear = TRUE, tol = 1e-7, iter = 500){

    # X: main predictor; Z: modifying variable; Y: response;
    # lambda: penalty parameter; alpha: weight between group penalty and individual penalty;
    # tt: learning rate for gradient descent
    # tol: tolerance
    # zlinear: if true, the linear terms of the modifying variables are included
    # iter: maximum number of iteration

    # Check that N is the same between Xtilde and Ytilde
    if (nrow(X) != length(Y) || nrow(X) != nrow(Z)){
        stop("Sample size of response and predictors are not equal.")
    }

    # Standardize inputs (Y: center, X & Z: center and scale)
    SXYZ <- standardizeXYZ(X, Z, Y)
    Xtilde <- SXYZ$Xtilde; Ztilde <- SXYZ$Ztilde; Ytilde <- SXYZ$Ytilde

    # p: number of main predictors; K: number of modifying variables, N: sample size
    p <- ncol(Xtilde)
    K <- ncol(Ztilde)
    N <- nrow(Xtilde)

    # W: tensor (list of matrix) of component-wise multiplication between X and Z (interaction of X and Z)
    # (same notation as in the paper "A pliable lasso")
    W <- get_empty_list(paste0("W_",1:p))
    for (i in 1:p){
        W[[i]] <- Xtilde[,i] * Ztilde
    }

    # beta: coefficient for main predictors; theta: coefficient for modifying variables
    if (is.null(beta)){
        beta <- matrix(0, ncol = 1, nrow = p)
    } else {
        beta <- beta
    }
    if (is.null(theta)){
        theta <- matrix(0, ncol = p, nrow = K)
    } else {
        theta <- theta
    }

    # Initialize the settings
    itr <- 0
    error <- 10000
    # full_res: residual of the current model before beta0 and theta0 calculation
    # full_res2: actual residual (used for full residual approach in the coordinate descent)
    full_res <- Ytilde - (Xtilde %*% beta)
    for (jj in 1:p){
        full_res <- full_res - (as.matrix(W[[jj]]) %*% theta[,jj])
    }
    Ytilde0 <- Ytilde
    lmmodel <- lm(full_res~Ztilde)
    if (zlinear == FALSE){
        beta0 <- mean(full_res)
        theta0 <- rep(0, K)
    } else {
        beta0 <- lmmodel$coefficients[1]
        theta0 <- lmmodel$coefficients[-1]
    }
    Ytilde <- Ytilde0 - beta0 - Ztilde %*% theta0
    full_res2 <- Ytilde

    while (error>tol && itr < iter){

        itr <- itr + 1
        beta_old <- beta
        theta_old <- theta
        beta0_old <- beta0
        theta0_old <- theta0

        for (j in 1:p){

            # check (beta,theta) = (0,0)
            b_tmp <- beta[j] + t(Xtilde[,j]) %*% (full_res2 + as.matrix(W[[j]]) %*% theta[,j])/N
            t_tmp <- crossprod(W[[j]]) %*% theta[,j]/N + t(as.matrix(W[[j]])) %*% (full_res2 + Xtilde[,j]*matrix(beta[j], N) )/N
            screen_cond_1 <- (abs(b_tmp) <= (1-alpha)*lambda)
            screen_cond_2 <- (sqrt(sum(soft_thresh(t_tmp, alpha*lambda)^2)) <= 2*(1-alpha)*lambda)
            if (screen_cond_1 == TRUE & screen_cond_2 == TRUE){
                # If (beta,theta) = (0,0), skip to the next predictor
                beta[j] <- 0
                theta[,j] <- 0
            } else {
                # If (beta,theta) != (0,0), compute beta_hat and check theta=0
                beta_check <- beta
                beta_check[j] <- N/sum(Xtilde[,j]^2) * soft_thresh(b_tmp, (1-alpha)*lambda)
                screen_cond_3 <- sqrt(sum(soft_thresh( t_tmp - (t(as.matrix(W[[j]])) %*% (Xtilde[,j]*matrix(beta_check[j], N)) )/N, alpha*lambda)^2)) <= (1-alpha)*lambda
                if (screen_cond_3 == TRUE){
                    # With beta_hat, if theta=0, set beta = beta_hat, theta=0 and skip to the next predictor
                    beta[j] <- beta_check[j]
                    theta[,j] <- 0
                } else{
                    # If beta != 0 and theta != 0, use gradient descent to update beta and theta
                    t <- tt
                    c <- t*(1-alpha)*lambda
                    res <- Ytilde - (Xtilde %*% beta)
                    for (jj in 1:p){
                        res <- res - (as.matrix(W[[jj]]) %*% theta[,jj])
                    }
                    grad_beta <- -(1/N) * ( t(Xtilde[,j]) %*% res )
                    g1 <- abs(beta[j] - t*grad_beta)
                    grad_theta <- -(1/N) * ( t(as.matrix(W[[j]])) %*% res )
                    g2 <- sqrt(sum(soft_thresh(theta[,j] - t*grad_theta, t*alpha*lambda)^2))
                    r1 <- -c+sqrt(c^2-2*c*g2+g1^2+g2^2)
                    r2 <- -c-sqrt(c^2-2*c*g2+g1^2+g2^2)

                    # a: norm of beta, b: norm of theta
                    # Hence, we choose the largest value of a and b to take positive value
                    a <- max(g1*r1/(c+r1), g1*r2/(c+r2), g1*r1/(c+r2), g1*r2/(c+r1))
                    b <- max((g2-c)*r2/(c+r2), (g2-c)*r1/(c+r1), (g2-c)*r1/(c+r2), (g2-c)*r2/(c+r1))
                    c1 <- 1+t*(1-alpha)*lambda/sqrt(a^2+b^2)
                    c2 <- 1+t*(1-alpha)*lambda*(1/b+1/sqrt(a^2+b^2))

                    beta[j] <- (beta[j] - t*grad_beta)/c1
                    theta[,j] <- soft_thresh(theta[,j] - t*grad_theta, t*alpha*lambda)/rep(c2, length(theta[,j]))
                }
            }
            full_res <- full_res + Xtilde[,j] * matrix((beta_old[j] - beta[j]), N) + W[[j]] %*% (theta_old[,j] - theta[,j])

            # beta0 and theta0 update
            lmmodel <- lm(full_res~Ztilde)
            if (zlinear == FALSE){
                beta0 <- mean(full_res)
                theta0 <- rep(0, K)
            } else {
                beta0 <- lmmodel$coefficients[1]
                theta0 <- lmmodel$coefficients[-1]
            }
            Ytilde <- Ytilde0 - beta0 - Ztilde %*% theta0

            # actual full residual calculation
            full_res2 <- full_res - beta0 - Ztilde %*% theta0
        }

        fmin=object_plasso(beta, theta, beta0, theta0, Xtilde, Ytilde0, Ztilde, W, alpha, lambda)
        error=abs(object_plasso(beta_old, theta_old, beta0_old, theta0_old, Xtilde, Ytilde0, Ztilde, W, alpha, lambda)-object_plasso(beta, theta, beta0, theta0, Xtilde, Ytilde0, Ztilde, W, alpha, lambda))
    }
    #print(c("iteration number: ",itr))

    beta_raw <- beta*(1/SXYZ$Xweights)
    theta_raw <- theta*matrix((1/SXYZ$Xweights), K, p, byrow = TRUE)*matrix((1/SXYZ$Zweights), K, p, byrow = FALSE)

    WW <- get_empty_list(paste0("WW_",1:p))
    for (i in 1:p){
        WW[[i]] <- X[,i] * Z
    }
    full_res_raw <- Y - (as.matrix(X) %*% beta_raw)
    for (jj in 1:p){
        full_res_raw <- full_res_raw - (as.matrix(WW[[jj]]) %*% theta_raw[,jj])
    }
    lmmodel_raw <- lm(full_res_raw~as.matrix(Z))
    if (zlinear == FALSE){
        beta0_raw <- mean(full_res_raw)
        theta0_raw <- rep(0, K)
    } else {
        beta0_raw <- lmmodel_raw$coefficients[1]
        theta0_raw <- lmmodel_raw$coefficients[-1]
    }

    return(list("average_coef"=c(beta, rowSums(theta)), "actual_coef"=list("main_coef"=beta, "modifying_coef"=theta), "raw_coef"=list("main_coef"=beta_raw, "modifying_coef"=theta_raw), "intercept"=list("beta0"=beta0, "theta0"=theta0), "intercept_raw"=list("beta0_raw"=beta0_raw, "theta0_raw"=theta0_raw), "fmin"=fmin))

}


#' Fit the pliable lasso in the linear regression setting
#'
#' @param X N by p matrix of main predictors
#' @param Z N by K matrix of modifying variables. Modifying variables can take the form of continuous variables or categorical variables or both. Categorical variable should be coded by dummy variables (0-1).
#' @param Y vector of response variable
#' @param lambda_seq sequence of the tuning parameter, lambda. Can take the form of a sequence or a scalar.
#' @param alpha weight parameter between group penalty and individual penalty. Default value is 0.5.
#' @param tt learning rate for the gradient descent procedure. Default value is 0.1.
#' @param zlinear if true, the linear terms of the modifying variables are included. These terms are not regularized. Default value is TRUE.
#' @param tol tolerance for convergence. Convergence is determined by the value of the objective function: abs(objective_old - objective_new) is compared with the tolerance value. Default value is 1e-7.
#' @param iter maximum number of iteration for one lambda. Default value is 500.
#'
#' @return lambda_seq: lambda sequence used in the algorithm
#' @return beta_mat: p by (length of lambda_seq) matrix of estimated beta for scaled and centered main predictors. Each column represents the vector of fitted beta for each lambda value. The order of lambda is the order of lambda_seq. For a scalar value of lambda_seq, the output is a p-dim vector of fitted beta.
#' @return theta_mat: p by K by (length of lambda_seq) array of estimated theta for scaled and centered main predictors and modifying variables.
#' @return beta0_vec: intercept term
#' @return theta0_vec: coefficient for the linear terms of the modifying variables. If zlinear = FALSE, the output is the vector of zeros.
#' @return beta_raw_mat: estimated beta for raw main predictors (non-standardized)
#' @return theta_raw_mat: estimated theta for raw modifying variables (non-standardized)
#' @return beta0_raw_vec: intercept term (non-standardized)
#' @return theta0_raw_vec: coefficient for the linear terms of the modifying variables (non-standardized)
#' @return fmin_vec: vector of objective function values for the lambda_seq values.
#' @export
#'
#' @examples
#' x=matrix(rnorm(100*5, 0, 1),100,5)
#' z=matrix(rnorm(100*3, 0, 1),100,3)
#' y=2*x[,1] - 2*x[,2] + (2+2*z[,1]-2*z[,2])*x[,3] + rnorm(100, 0, 1)
#' plasso(X = x, Z = z, Y = y, lambda_seq = c(1, 0.5), alpha = 0.5)
#' plasso(X = x, Z = z, Y = y, lambda_seq = 0.5, alpha = 0.5)
#' plasso(X = x, Z = z, Y = y, lambda_seq = 0.5, alpha = 0.5, zlinear = FALSE)
plasso <- function(X, Z, Y, lambda_seq = NULL, alpha = 0.5, tt = 0.1, zlinear = TRUE, tol = 1e-7, iter = 500){

    # p: number of main predictors; K: number of modifying variables, N: sample size
    p <- ncol(X)
    K <- ncol(Z)
    N <- nrow(X)

    if (!is.null(lambda_seq)){
        lambda_seq <- sort(lambda_seq[lambda_seq >= 0], decreasing = TRUE)
        if (length(lambda_seq) == 0){
            stop("All lambda values are negative")
        }
    } else {
        stop("lambda_seq must be specified")
    }

    para_array=array(NA, c(p, K+1,length(lambda_seq))); para_array_raw=array(NA, c(p, K+1,length(lambda_seq)))
    fmin_vec=rep(NA,length(lambda_seq))
    beta0_vec=rep(NA,length(lambda_seq)); beta0_raw_vec=rep(NA,length(lambda_seq))
    theta0_vec=matrix(NA, K, length(lambda_seq)); theta0_raw_vec=matrix(NA, K, length(lambda_seq))

    # Starting beta and theta of Warm Start is zero vector and matrix
    print(c("lambda: ",lambda_seq[1]))
    fit <- plasso1(X, Z, Y, lambda_seq[1], alpha = alpha, tt = tt, beta = NULL, theta = NULL, zlinear = zlinear, tol = tol, iter = iter)
    para_array[,1,1] <- fit$actual_coef$main_coef; para_array[,-1,1] <- t(fit$actual_coef$modifying_coef)
    para_array_raw[,1,1] <- fit$raw_coef$main_coef; para_array_raw[,-1,1] <- t(fit$raw_coef$modifying_coef)
    beta0_vec[1] <- fit$intercept$beta0; theta0_vec[,1] <- fit$intercept$theta0
    beta0_raw_vec[1] <- fit$intercept_raw$beta0_raw; theta0_raw_vec[,1] <- fit$intercept_raw$theta0_raw
    fmin_vec[1] <- fit$fmin
    # Carry over previous beta for Warm Start
    if (length(lambda_seq) > 1){
        for (i in 2:length(lambda_seq)){
            print(c("lambda: ",lambda_seq[i]))
            fit <- plasso1(X, Z, Y, lambda_seq[i], alpha = alpha, tt = tt, beta = para_array[,1,i-1], theta = t(para_array[,-1,i-1]), zlinear = zlinear, tol = tol, iter = iter)
            para_array[,1,i] <- fit$actual_coef$main_coef; para_array[,-1,i] <- t(fit$actual_coef$modifying_coef)
            para_array_raw[,1,i] <- fit$raw_coef$main_coef; para_array_raw[,-1,i] <- t(fit$raw_coef$modifying_coef)
            beta0_vec[i] <- fit$intercept$beta0; theta0_vec[,i] <- fit$intercept$theta0
            beta0_raw_vec[i] <- fit$intercept_raw$beta0_raw; theta0_raw_vec[,i] <- fit$intercept_raw$theta0_raw
            fmin_vec[i] <- fit$fmin
        }
    }

    return(list(lambda_seq = lambda_seq, beta_mat = para_array[,1,], theta_mat = para_array[,-1,], beta0_vec = beta0_vec, theta0_vec = theta0_vec, beta_raw_mat = para_array_raw[,1,], theta_raw_mat = para_array_raw[,-1,], beta0_raw_vec = beta0_raw_vec, theta0_raw_vec = theta0_raw_vec, fmin_vec = fmin_vec))

}


#' Perform k-fold cross validation for the pliable lasso model over a sequence of regularization parameter
#'
#' @param X N by p matrix of main predictors
#' @param Z N by K matrix of modifying variables. Modifying variables can take the form of continuous variables or categorical variables or both. Categorical variable should be coded by dummy variables (0-1).
#' @param Y vector of response variable
#' @param kfold the number of folds (=k) for the k-fold cross-validation. Default value is 10.
#' @param lambda_seq sequence of the tuning parameter, lambda. Can take the form of a sequence or a scalar.
#' @param alpha weight parameter between group penalty and individual penalty. Default value is 0.5.
#' @param tt learning rate for the gradient descent procedure. Default value is 0.1.
#' @param zlinear if true, the linear terms of the modifying variables are included. These terms are not regularized. Default value is TRUE.
#' @param tol tolerance for convergence. Convergence is determined by the value of the objective function: abs(objective_old - objective_new) is compared with the tolerance value. Default value is 1e-7.
#' @param iter maximum number of iteration for one lambda. Default value is 500.
#' @param cvseed if specified, seed number for random sampling in the cross-validation procedure is fixed so the result is reproducible. If unspecified, the result is not reproducible. Default value is NULL.
#'
#' @return lambda_seq: lambda sequence used in the algorithm
#' @return beta_mat: p by (length of lambda_seq) matrix of estimated beta for scaled and centered main predictors. Each column represents the vector of fitted beta for each lambda value. The order of lambda is the order of lambda_seq. For a scalar value of lambda_seq, the output is a p-dim vector of fitted beta.
#' @return theta_mat: p by K by (length of lambda_seq) array of estimated theta for scaled and centered main predictors and modifying variables.
#' @return beta0_vec: intercept term
#' @return theta0_vec: coefficient for the linear terms of the modifying variables. If zlinear = FALSE, the output is the vector of zeros.
#' @return beta_raw_mat: estimated beta for raw main predictors (non-standardized)
#' @return theta_raw_mat: estimated theta for raw modifying variables (non-standardized)
#' @return beta0_raw_vec: intercept term (non-standardized)
#' @return theta0_raw_vec: coefficient for the linear terms of the modifying variables (non-standardized)
#' @return lambda_min: the lambda value which minimizes the pliable lasso objective function among the values in the lambda_seq.
#' @return lambda_1se: the largest lambda value such that the difference with minimum objective function value is within 1 standard error of the minimum
#' @return cvm: the sequence of objective function values for lambda_seq
#' @return cvse: the sequence of the standard error of objective function values for lambda_seq
#' @return cvfold: (kfold) by (length of lambda_seq) matrix of the mean squared error of the test set for each fold
#' @return sqerror: N by (length of lambda_seq) matrix of the squared error
#' @export
#'
#' @examples
#' x=matrix(rnorm(100*5, 0, 1),100,5)
#' z=matrix(rnorm(100*3, 0, 1),100,3)
#' y=2*x[,1] - 2*x[,2] + (2+2*z[,1]-2*z[,2])*x[,3] + rnorm(100, 0, 1)
#' cv.plasso(X = x, Z = z, Y = y, lambda_seq = c(1, 0.5), alpha = 0.5)
#' cv.plasso(X = x, Z = z, Y = y, lambda_seq = c(1, 0.5), alpha = 0.5, cvseed = 1234)
#' cv.plasso(X = x, Z = z, Y = y, lambda_seq = 0.5, alpha = 0.5)
#' cv.plasso(X = x, Z = z, Y = y, lambda_seq = 0.5, alpha = 0.5, zlinear = FALSE)
cv.plasso <- function(X, Z, Y, kfold = 10, lambda_seq = NULL, alpha = 0.5, tt = 0.1, zlinear = TRUE, tol = 1e-7, iter = 500, cvseed = NULL){

    # p: number of main predictors; K: number of modifying variables, N: sample size
    p <- ncol(X)
    K <- ncol(Z)
    N <- nrow(X)

    # Fit Pliable Lasso on original data using plasso
    plassofit <- plasso(X, Z, Y, lambda_seq = lambda_seq, alpha = alpha, tt = tt, zlinear = zlinear, tol = tol, iter = iter)

    # Split the data into K folds
    if (!is.null(cvseed)) set.seed(cvseed)
    idfold <- sample(1:N) %% kfold + 1

    # Calculate Pliable Lasso for each fold removed
    n_lambda <- length(plassofit$lambda_seq)
    sqerror <- matrix(NA, N, n_lambda)
    cvfold <- matrix(NA, kfold, n_lambda)
    for (fold in 1:kfold){

        print(c("fold: ", fold))

        # Training data
        xtrain = X[idfold != fold, ]
        ztrain = Z[idfold != fold, ]
        ytrain = Y[idfold != fold]

        # Test data
        xtest = X[idfold == fold, ]
        ztest = Z[idfold == fold, ]
        ytest = Y[idfold == fold]

        SXYZtest <- standardizeXYZ(xtest, ztest, ytest)
        xtest <- SXYZtest$Xtilde; ztest <- SXYZtest$Ztilde; ytest <- SXYZtest$Ytilde

        # Calculate LASSO on that fold using fitLASSO
        cvfit <- plasso(xtrain, ztrain, ytrain, lambda_seq, alpha = alpha, tt = tt, zlinear = zlinear, tol = tol, iter = iter)

        # Any additional calculations that are needed for calculating CV and SE_CV(lambda)
        testfitted <- matrix(rep(cvfit$beta0_vec, length(ytest)), nrow = length(ytest), byrow = TRUE) + ztest %*% cvfit$theta0_vec + xtest %*% cvfit$beta_mat
        for (j in 1:p){
            if (n_lambda == 1){
                testfitted <- testfitted + (xtest[,j]*ztest) %*% cvfit$theta_mat[j,]
            } else{
                testfitted <- testfitted + (xtest[,j]*ztest) %*% cvfit$theta_mat[j,,]
            }
        }
        cvfold[fold, ] <- colMeans((ytest - testfitted)^2)
        sqerror[idfold == fold, ] <- (ytest - testfitted)^2
    }

    # Calculate CV(lambda) and SE_CV(lambda) for each value of lambda
    cvm <- colMeans(sqerror)
    cvse <- apply(cvfold, 2, sd)/sqrt(kfold)

    # Find lambda_min
    lambda_min <- plassofit$lambda_seq[which.min(cvm)]

    # Find lambda_1SE
    lambda_1se <- plassofit$lambda_seq[cvm <= (min(cvm) + cvse[which.min(cvm)])][1]


    lambda_seq <- plassofit$lambda_seq
    beta_mat <- plassofit$beta_mat
    beta0_vec <- plassofit$beta0_vec
    theta_mat <- plassofit$theta_mat
    theta0_vec <- plassofit$theta0_vec
    beta_raw_mat <- plassofit$beta_raw_mat
    beta0_raw_vec <- plassofit$beta0_raw_vec
    theta_raw_mat <- plassofit$theta_raw_mat
    theta0_raw_vec <- plassofit$theta0_raw_vec

    return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, theta_mat = theta_mat, beta0_vec = beta0_vec, theta0_vec = theta0_vec, beta_raw_mat = beta_raw_mat, beta0_raw_vec = beta0_raw_vec, theta_raw_mat = theta_raw_mat, theta0_raw_vec = theta0_raw_vec, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse, cvfold = cvfold, sqerror = sqerror))
}
