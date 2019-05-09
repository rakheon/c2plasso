c2plasso1 <- function(X, Z, Y, df_Z, lambda = 0.5, alpha = 0.5, tt = 0.1, beta = NULL, theta = NULL, zlinear = TRUE, tol = 1e-7, iter = 500){

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

    ## check degrees of freedom of modifying variables
    df_Z <- round(df_Z,0)
    error_positive_value(df_Z)
    error_modifying_variable_df(sum(df_Z), ncol(Z))

    df_Z_cum <- cumsum(df_Z)
    group_partition <- get_empty_list(paste0("g_",1:length(df_Z)))
    group_partition[[1]] <- 1:df_Z_cum[1]
    if (length(df_Z) > 1){
        for (i in 2:length(df_Z)){
            group_partition[[i]] <- (df_Z_cum[i-1]+1):df_Z_cum[i]
        }
    }

    # Standardize inputs (Y: center, X & Z: center and scale)
    # If Z = NULL, plasso is equivalent to plain lasso
    SXYZ <- standardizeXYZ(X, Z, Y)
    Xtilde <- SXYZ$Xtilde; Ztilde <- SXYZ$Ztilde; Ytilde <- SXYZ$Ytilde

    # p: number of main predictors; K: number of modifying variables, N: sample size
    p <- ncol(Xtilde)
    K <- ncol(Ztilde)
    N <- nrow(Xtilde)
    G <- length(df_Z)

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
            tg_tmp <- get_empty_list(paste0("g_",1:G))
            for (kk in 1:G){
                tg_tmp[[kk]] <- crossprod(W[[j]][,group_partition[[kk]]]) %*% theta[,j][group_partition[[kk]]]/N + t(as.matrix(W[[j]][,group_partition[[kk]]])) %*% (full_res2 + Xtilde[,j]*matrix(beta[j], N) )/N
            }
            screen_cond_1 <- (abs(b_tmp) <= (1-alpha)*lambda)
            screen_cond_2 <- (sqrt(sum(soft_thresh(t_tmp, alpha*lambda)^2)) <= (1+sum(sqrt(df_Z))/sqrt(1+K))*(1-alpha)*lambda)
            if (screen_cond_1 == TRUE & screen_cond_2 == TRUE){
                # If (beta,theta) = (0,0), skip to the next predictor
                beta[j] <- 0
                theta[,j] <- 0
            } else {
                # If (beta,theta) != (0,0), compute beta_hat and check theta=0
                beta_check <- beta
                beta_check[j] <- N/sum(Xtilde[,j]^2) * soft_thresh(b_tmp, (1-alpha)*lambda)

                screen_cond_3G <- logical(G)
                screen_cond_3G_new <- !logical(G)
                while (sum(ifelse(screen_cond_3G==screen_cond_3G_new,1,0)) < G){

                    screen_cond_3G <- screen_cond_3G_new

                    for (kk in 1:G){
                        screen_cond_3G[kk] <- (sqrt(sum(soft_thresh( tg_tmp[[kk]] - ( t(as.matrix(W[[j]][,group_partition[[kk]]])) %*% (Xtilde[,j]*matrix(beta_check[j], N)) )/N, alpha*lambda)^2)) <= (1-alpha)*lambda*sqrt(df_Z[kk])/sqrt(1+K))
                    }

                    zG <- sum(screen_cond_3G)
                    nzG <- G - zG
                    nzdf_Z <- df_Z[!screen_cond_3G]

                    if (nzG==0){
                        # With beta_hat, if theta=0 (i.e. no non-zero group of Z), set beta = beta_hat, theta=0 and skip to the next predictor
                        beta[j] <- beta_check[j]
                        theta[,j] <- 0
                    } else{

                        # If beta != 0 and theta != 0, use gradient descent to update beta and theta for the non-zero group of Z
                        if (zG!=0){
                            z_group_partition <- group_partition[screen_cond_3G]
                            for (kk in 1:zG){
                                theta[,j][z_group_partition[[kk]]] <- 0
                            }
                        }

                        nz_group_partition <- group_partition[!screen_cond_3G]
                        t <- tt
                        c <- t*(1-alpha)*lambda
                        res <- Ytilde - (Xtilde %*% beta)
                        for (jj in 1:p){
                            res <- res - (as.matrix(W[[jj]]) %*% theta[,jj])
                        }
                        grad_beta <- -(1/N) * ( t(Xtilde[,j]) %*% res )
                        g1 <- abs(beta[j] - t*grad_beta)
                        grad_theta <- get_empty_list(paste0("grad_theta_",1:nzG))
                        g2 <- numeric(nzG)
                        for (i in 1:nzG){
                            grad_theta[[i]] <- -(1/N) * ( t(as.matrix(W[[j]][,nz_group_partition[[i]]])) %*% res )
                            g2[i] <- sqrt(sum(soft_thresh(theta[,j][nz_group_partition[[i]]] - t*grad_theta[[i]], t*alpha*lambda)^2))
                        }
                        r1 <- -c+sqrt(c^2-2*c*sum(g2*sqrt(nzdf_Z)/sqrt(1+K))+g1^2+sum(g2^2)-(1-sum(nzdf_Z)/(1+K))*(c^2))
                        r2 <- -c-sqrt(c^2-2*c*sum(g2*sqrt(nzdf_Z)/sqrt(1+K))+g1^2+sum(g2^2)-(1-sum(nzdf_Z)/(1+K))*(c^2))
                        # a: norm of beta, b: norm of theta
                        # Hence, we choose the largest value of a and b to take positive value
                        #a <- max(g1*r1/(c+r1), g1*r2/(c+r2), g1*r1/(c+r2), g1*r2/(c+r1))
                        #b <- max((g2-c)*r2/(c+r2), (g2-c)*r1/(c+r1), (g2-c)*r1/(c+r2), (g2-c)*r2/(c+r1))
                        a <- max(g1*r1/(c+r1), g1*r2/(c+r2), g1*r1/(c+r2), g1*r2/(c+r1))
                        b <- numeric(nzG)
                        for (i in 1:nzG){
                            b[i] <- max((g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r2/(c+r2), (g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r1/(c+r1), (g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r1/(c+r2), (g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r2/(c+r1))
                        }
                        c1 <- 1+t*(1-alpha)*lambda/sqrt(a^2+sum(b^2))
                        c2 <- numeric(nzG)
                        for (i in 1:nzG){
                            c2[i] <- 1+t*(1-alpha)*lambda*(sqrt(nzdf_Z[i])/sqrt(1+K)*1/b[i]+1/sqrt(a^2+sum(b^2)))
                        }
                        beta[j] <- (beta[j] - t*grad_beta)/c1
                        for (i in 1:nzG){
                            theta[,j][nz_group_partition[[i]]] <- soft_thresh(theta[,j][nz_group_partition[[i]]] - t*grad_theta[[i]], t*alpha*lambda)/rep(c2[i], length(theta[,j][nz_group_partition[[i]]]))
                        }

                        beta_check[j] <- beta[j]
                        for (kk in 1:G){
                            screen_cond_3G_new[kk] <- (sqrt(sum(soft_thresh( tg_tmp[[kk]] - ( t(as.matrix(W[[j]][,group_partition[[kk]]])) %*% Xtilde[,j]*beta_check[j] )/N, alpha*lambda)^2)) <= (1-alpha)*lambda*sqrt(df_Z[kk])/sqrt(1+K))
                        }
                    }
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

        fmin=object_c2plasso(beta, theta, beta0, theta0, Xtilde, Ytilde0, Ztilde, W, group_partition, G, alpha, lambda)
        error=abs(object_c2plasso(beta_old, theta_old, beta0_old, theta0_old, Xtilde, Ytilde0, Ztilde, W, group_partition, G, alpha, lambda)-object_c2plasso(beta, theta, beta0, theta0, Xtilde, Ytilde0, Ztilde, W, group_partition, G, alpha, lambda))
    }
    print(c("iteration number: ",itr))

    beta_raw <- beta*(1/SXYZ$Xweights)
    theta_raw <- theta*matrix((1/SXYZ$Xweights), K, p, byrow = TRUE)*matrix((1/SXYZ$Zweights), K, p, byrow = FALSE)

    return(list("average_coef"=c(beta, rowSums(theta)), "actual_coef"=list("main_coef"=beta, "modifying_coef"=theta), "raw_coef"=list("main_coef"=beta_raw, "modifying_coef"=theta_raw), "intercept"=list("beta0"=beta0, "theta0"=theta0), "fmin"=fmin))

}


c2plasso <- function(X, Z, Y, df_Z, lambda_seq = NULL, alpha = 0.5, tt = 0.1, zlinear = TRUE, tol = 1e-7, iter = 500){

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

    para_array=array(NA, c(p, K+1,length(lambda_seq)))
    fmin_vec=rep(NA,length(lambda_seq))
    beta0_vec=rep(NA,length(lambda_seq))
    theta0_vec=matrix(NA, K, length(lambda_seq))

    # Starting beta and theta of Warm Start is zero vector and matrix
    fit <- c2plasso1(X, Z, Y, df_Z, lambda_seq[1], alpha = alpha, tt = tt, beta = NULL, theta = NULL, zlinear = zlinear, tol = tol, iter = iter)
    para_array[,1,1] <- fit$actual_coef$main_coef; para_array[,-1,1] <- t(fit$actual_coef$modifying_coef)
    beta0_vec[1] <- fit$intercept$beta0; theta0_vec[,1] <- fit$intercept$theta0
    fmin_vec[1] <- fit$fmin
    # Carry over previous beta for Warm Start
    if (length(lambda_seq) > 1){
        for (i in 2:length(lambda_seq)){
            fit <- c2plasso1(X, Z, Y, df_Z, lambda_seq[i], alpha = alpha, tt = tt, beta = para_array[,1,i-1], theta = t(para_array[,-1,i-1]), zlinear = zlinear, tol = tol, iter = iter)
            para_array[,1,i] <- fit$actual_coef$main_coef; para_array[,-1,i] <- t(fit$actual_coef$modifying_coef)
            beta0_vec[i] <- fit$intercept$beta0; theta0_vec[,i] <- fit$intercept$theta0
            fmin_vec[i] <- fit$fmin
        }
    }

    return(list(lambda_seq = lambda_seq, beta_mat = para_array[,1,], theta_mat = para_array[,-1,], beta0_vec = beta0_vec, theta0_vec = theta0_vec, fmin_vec = fmin_vec))

}


cv.c2plasso <- function(X, Z, Y, df_Z, kfold = 10, lambda_seq = NULL, alpha = 0.5, tt = 0.1, zlinear = TRUE, tol = 1e-7, iter = 500, cvseed = NULL){

    # p: number of main predictors; K: number of modifying variables, N: sample size
    p <- ncol(X)
    K <- ncol(Z)
    N <- nrow(X)

    # Fit Pliable Lasso on original data using plasso
    c2plassofit <- c2plasso(X, Z, Y, df_Z, lambda_seq = lambda_seq, alpha = alpha, tt = tt, zlinear = zlinear, tol = tol, iter = iter)

    # Split the data into K folds
    if (!is.null(cvseed)) set.seed(cvseed)
    idfold <- sample(1:N) %% kfold + 1

    # Calculate Pliable Lasso for each fold removed
    n_lambda <- length(c2plassofit$lambda_seq)
    sqerror <- matrix(NA, N, n_lambda)
    cvfold <- matrix(NA, kfold, n_lambda)
    for (fold in 1:kfold){

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
        cvfit <- c2plasso(xtrain, ztrain, ytrain, df_Z, lambda_seq, alpha = alpha, tt = tt, zlinear = zlinear, tol = tol, iter = iter)

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
    lambda_min <- c2plassofit$lambda_seq[which.min(cvm)]

    # Find lambda_1SE
    lambda_1se <- c2plassofit$lambda_seq[cvm <= (min(cvm) + cvse[which.min(cvm)])][1]


    lambda_seq <- c2plassofit$lambda_seq
    beta_mat <- c2plassofit$beta_mat
    beta0_vec <- c2plassofit$beta0_vec
    theta_mat <- c2plassofit$theta_mat
    theta0_vec <- c2plassofit$theta0_vec

    return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, theta_mat = theta_mat, beta0_vec = beta0_vec, theta0_vec = theta0_vec, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse, cvfold = cvfold, sqerror = sqerror))
}
