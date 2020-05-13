object_onedim <- function(bt, bt0, Xtilde, Ytilde, lambda){

    p <- ncol(Xtilde)
    N <- nrow(Xtilde)

    diff <- Ytilde - Xtilde %*% bt - bt0
    obj <- (1/(2*N))*sum((diff)**2) + lambda*(sqrt(sum(bt^2)))
    return(obj)
}
object_onedim1 <- function(bt, ot, Xtilde_b, Xtilde_o, Ytilde, lambda){

    N <- length(Ytilde)
    #print(Xtilde_o)
    #print(ot)
    diff <- Ytilde - Xtilde_b * bt - as.matrix(Xtilde_o) %*% ot
    obj <- (1/(2*N))*sum((diff)**2) + lambda*(sqrt(bt^2+sum(ot^2)))
    return(obj)
}

svReg1 <- function(X, Z = NULL, Y, df_X, df_Z, lambda = 0.5, alpha = 0.5, tt = 0.1, beta = NULL, theta = NULL, zlinear = TRUE, tol = 1e-7, iter = 500){

    # X: main predictor; Z: modifying variable; Y: response;
    # lambda: penalty parameter; alpha: weight between group penalty and individual penalty;
    # tol: tolerance

    ## check sample_size >=1
    df_X <- round(df_X,0)
    df_Z <- round(df_Z,0)
    #error_positive_value(df_Z,get_variable_name(df_Z))

    #error_modifying_variable_df(sum(df_Z), ncol(Z))

    df_X_cum <- cumsum(df_X)
    df_Z_cum <- cumsum(df_Z)
    main_partition <- get_empty_list(paste0("l_",1:length(df_X)))
    main_partition[[1]] <- 1:df_X_cum[1]
    for (i in 2:length(df_X)){
        main_partition[[i]] <- (df_X_cum[i-1]+1):df_X_cum[i]
    }
    group_partition <- get_empty_list(paste0("g_",1:length(df_Z)))
    group_partition[[1]] <- 1:df_Z_cum[1]
    for (i in 2:length(df_Z)){
        group_partition[[i]] <- (df_Z_cum[i-1]+1):df_Z_cum[i]
    }

    # Standardize inputs (Y: center, X & Z: center and scale)
    # If Z = NULL, plasso is equivalent to plain lasso
    if (is.null(Z)){
        SXYZ <- standardizeXYZ(X, Z, Y)
        Xtilde <- SXYZ$Xtilde; Ytilde <- SXYZ$Ytilde
        Ztilde <- matrix(rep(0,nrow(X)), ncol=1)
    } else {
        SXYZ <- standardizeXYZ(X, Z, Y)
        Xtilde <- SXYZ$Xtilde; Ztilde <- SXYZ$Ztilde; Ytilde <- SXYZ$Ytilde
    }

    # p: number of main predictors; K: number of modifying variables, N: sample size, L: number of main groups, G: number of modifying groups
    p <- ncol(Xtilde)
    K <- ncol(Ztilde)
    N <- nrow(Xtilde)
    L <- length(df_X)
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
        #beta <- matrix(0, ncol = 1, nrow = p)
    }
    if (is.null(theta)){
        theta <- matrix(0, ncol = p, nrow = K)
    } else {
        theta <- theta
        #theta <- matrix(0, ncol = p, nrow = K)
    }

    itr <- 0
    error <- 10000
    #fmin <- numeric(iter)
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
    Ytilde <- Ytilde0 - beta0 - Ztilde %*% theta0 - Xtilde %*% beta
    for (jj in 1:p){
        Ytilde <- Ytilde - (as.matrix(W[[jj]]) %*% theta[,jj])
    }
    full_res2 <- Ytilde
    #xsqmean <- colMeans(Xtilde^2)
    while (error>tol && itr < iter){
        #while ((sum(beta, theta)==0 || sum((beta-beta_old)^2)>tol || sum((theta-theta_old)^2)>tol) && itr < iter){

        itr <- itr + 1
        beta_old <- beta
        theta_old <- theta
        beta0_old <- beta0
        theta0_old <- theta0

        for (l in 1:L){

            main_group <- main_partition[[l]]
            j <- main_group
            # check (beta,theta) = (0,0)
            b_update = full_res2*0
            for (m in 1:length(j)){
                b_update = b_update + as.matrix(W[j][[m]]) %*% as.matrix(as.matrix(theta[,j])[,m])
            }
            #b_tmp <- beta[j] + solve( crossprod(Xtilde[,j])/N ) %*% t(Xtilde[,j]) %*% (full_res2 + b_update)/N
            #b_tmp <- beta[j] + t(Xtilde[,j]) %*% (full_res2 + b_update)/N
            t_update = full_res2*0
            for (m in 1:length(j)){
                t_update = t_update + as.matrix(as.matrix(Xtilde[,j])[,m])*matrix(beta[j][m], N)
            }
            #b_tmp_f <- beta[j] + t(Xtilde[,j]) %*% (full_res2)/N
            #print(c("full res:", b_tmp))
            b_tmp <- t(Xtilde[,j]) %*% (full_res2 + b_update + t_update)/N
            #print(c("partial res:", b_tmp_p))
            #b_tmp <- beta[j] + t(Xtilde[,j]) %*% (full_res2 + b_update)/N
            #b_tmp <- beta[j] + t(Xtilde[,j]) %*% (full_res2 + b_update + t_update)/N
            tg_tmp <- get_empty_list(paste0("g_",1:G))
            #tg_tmp_f <- get_empty_list(paste0("g_",1:G))
            screen_cond_2 <- logical(G)
            for (kk in 1:G){
                tg_tmp[[kk]] <- matrix(0, nrow = length(group_partition[[kk]]), ncol = length(j))
                #tg_tmp_f[[kk]] <- matrix(0, nrow = length(group_partition[[kk]]), ncol = length(j))
                screen_cond_2_val = tg_tmp[[kk]]*0
                for (m in 1:length(j)){
                    #tg_tmp_f[[kk]][,m] <- crossprod(W[j][[m]][,group_partition[[kk]]]) %*% as.matrix(as.matrix(theta[,j])[,m])[group_partition[[kk]]]/N + t(as.matrix(W[j][[m]][,group_partition[[kk]]])) %*% (full_res2)/N
                    #print(c("full res tg:", tg_tmp[[kk]][,m]))
                    tg_update = b_update - (as.matrix(W[j][[m]][,group_partition[[kk]]]) %*% as.matrix(as.matrix(theta[,j])[,m])[group_partition[[kk]]])
                    tg_tmp[[kk]][,m] <- t(as.matrix(W[j][[m]][,group_partition[[kk]]])) %*% (full_res2 + b_update + t_update - tg_update)/N
                    #print(c("partial res tg:", tg_tmp_p[[kk]][,m]))
                    #tg_tmp[[kk]][,m] <- crossprod(W[j][[m]][,group_partition[[kk]]]) %*% as.matrix(as.matrix(theta[,j])[,m])[group_partition[[kk]]]/N + t(as.matrix(W[j][[m]][,group_partition[[kk]]])) %*% (full_res2 + t_update)/N
                    #tg_tmp[[kk]][,m] <- crossprod(W[j][[m]][,group_partition[[kk]]]) %*% as.matrix(as.matrix(theta[,j])[,m])[group_partition[[kk]]]/N + t(as.matrix(W[j][[m]][,group_partition[[kk]]])) %*% (full_res2 + b_update + t_update)/N
                    #tg_tmp[[kk]][,m] <- crossprod(W[j][[m]][,group_partition[[kk]]]) %*% as.matrix(as.matrix(theta[,j])[,m])[group_partition[[kk]]]/N + t(as.matrix(W[j][[m]][,group_partition[[kk]]])) %*% (full_res2 + as.matrix(as.matrix(Xtilde[,j])[,m])*matrix(beta[j][m], N) )/N
                    screen_cond_2_val[,m] <- soft_thresh(tg_tmp[[kk]][,m], alpha*lambda)^2
                }
                #tg_tmp[[kk]] <- crossprod(W[[j]][,group_partition[[kk]]]) %*% theta[,j][group_partition[[kk]]]/N + t(as.matrix(W[[j]][,group_partition[[kk]]])) %*% (full_res2 + Xtilde[,j]*matrix(beta[j], N) )/N
                screen_cond_2[kk] <- (sqrt(sum(screen_cond_2_val)) <= (sqrt(df_X[l])+sqrt(df_X[l]*df_Z[kk])/sqrt(1+K))*(1-alpha)*lambda)
                #screen_cond_2[kk] <- (sqrt(sum(screen_cond_2_val)) <= (sqrt(df_X[l])+sqrt(df_Z[kk])/sqrt(1+K))*(1-alpha)*lambda)
                #screen_cond_2[kk] <- (sqrt(sum(screen_cond_2_val)) <= (1+sqrt(df_Z[kk])/sqrt(1+K))*(1-alpha)*lambda)
            }
            #print(b_tmp)
            screen_cond_1 <- (sqrt(sum(b_tmp^2)) <= sqrt(df_X[l])*(1-alpha)*lambda)
            #screen_cond_1 <- (sqrt(sum(b_tmp^2)) <= (1-alpha)*lambda)
            #screen_cond_1 <- (abs(b_tmp) <= (1-alpha)*lambda)
            #print(screen_cond_1)
            #print(screen_cond_2)
            if (screen_cond_1 == TRUE & prod(screen_cond_2) == TRUE){
                #print("start step1")
                # If (beta,theta) = (0,0), skip to the next predictor
                beta[j] <- beta[j]*0
                theta[,j] <- theta[,j]*0
            } else {
                #print("start step2")
                #b_tmp <- beta[j] + solve( crossprod(Xtilde[,j])/N ) %*% t(Xtilde[,j]) %*% (full_res2 + b_update + t_update)/N
                #b_tmp <- beta[j] + t(Xtilde[,j]) %*% (full_res2 + b_update + t_update)/N
                # If (beta,theta) != (0,0), compute beta_hat and check theta=0
                beta_check <- beta
                if (TRUE){
                    if (df_X[l]>1){
                        x=as.matrix(Xtilde[,j])
                        y=full_res2 + b_update + t_update
                        betain = solve(t(x)%*%x)%*%t(x)%*%y #matrix(0, nrow = p, ncol = 1)
                        nn = dim(x)[1]; pp = dim(x)[2]
                        lambdas = lambda*(1-alpha)*sqrt(df_X[l])
                        obj = 0; obj_old=1
                        while (abs(obj-obj_old) > 1e-5){
                            obj_old = obj
                            for (i in 1:pp){
                                betain[i,] = stats::optimize(object_onedim1, c(-100,100), ot=betain[-i,], Xtilde_b=x[,i], Xtilde_o=x[,-i], Ytilde=y, lambda=lambdas)$minimum
                            }
                            obj = object_onedim(bt=betain, bt0=0, Xtilde=x, Ytilde=y, lambda=lambdas)
                            #print(abs(obj-obj_old))
                        }
                        beta_check[j] <- betain
                    } else {
                        #b_tmp <- t(Xtilde[,j]) %*% (full_res2 + b_update + t_update)/N
                        beta_check[j] <- N/sum(Xtilde[,j]^2) * soft_thresh(b_tmp, sqrt(df_X[l])*(1-alpha)*lambda)
                    }
                }
                if (FALSE){
                if (df_X[l]>1){
                    ### Use group lasso
                    # Assuming theta = 0, partial residual is fullres + beta_l
                    grp_data = list(x=as.matrix(Xtilde[,j]), y=full_res2 + b_update + t_update, n=N)
                    #grp_data = list(x=as.matrix(Xtilde[,j]), y=full_res2 + b_update, n=N)
                    #grp_data = list(x=as.matrix(Xtilde[,j]), y=full_res2 + t_update, n=N)
                    grp_index = rep(1, df_X[l])
                    grp_pen = (1-alpha)*lambda
                    #sgl <- SGL(grp_data, grp_index, type = "linear", nlam = 1, alpha = 0, lambdas = grp_pen, standardize = FALSE)
                    sgl <- oneDim(data = grp_data, index = grp_index, alpha = 0, lambdas = grp_pen)
                    #print(sgl$beta)
                    beta_check[j] <- sgl$beta
                } else {
                    #b_tmp <- t(Xtilde[,j]) %*% (full_res2 + b_update + t_update)/N
                    beta_check[j] <- N/sum(Xtilde[,j]^2) * soft_thresh(b_tmp, sqrt(df_X[l])*(1-alpha)*lambda)
                }
                }
                #print(c("b_tmp", b_tmp))
                #b_tmp_scale <- solve( crossprod(Xtilde[,j])/N ) %*% b_tmp
                #beta_check[j] <- max(1-sqrt(df_X[l])*(1-alpha)*lambda/sqrt(sum(b_tmp_scale^2)), 0)*b_tmp_scale
                #beta_check[j] <- solve( crossprod(Xtilde[,j])/N ) %*% (max(1-sqrt(df_X[l])*(1-alpha)*lambda/sqrt(sum(b_tmp^2)), 0)*b_tmp)
                #beta_check[j] <- N/sum(Xtilde[,j]^2) * max(1-sqrt(df_X[l])*(1-alpha)*lambda/sqrt(sum(b_tmp^2)), 0)*b_tmp
                # approximate for non-orthogonal case
                #beta_check[j] <- max(1-sqrt(df_X[l])*(1-alpha)*lambda/sqrt(sum(b_tmp^2)), 0)*(solve( crossprod(Xtilde[,j])/N ) %*% b_tmp)
                # orthogonal case
                #beta_check[j] <- max(1-sqrt(df_X[l])*(1-alpha)*lambda/sqrt(sum(b_tmp^2)), 0)*b_tmp
                #print(beta_check[j])
                #beta_check[j] <- N/sum(Xtilde[,j]^2) * soft_thresh(b_tmp, sqrt(df_X[l])*(1-alpha)*lambda)
                #beta_check[j] <- N/sum(Xtilde[,j]^2) * soft_thresh(b_tmp, (1-alpha)*lambda)

                screen_cond_3G <- logical(G)
                screen_cond_3G_new <- !logical(G)
                while (sum(ifelse(screen_cond_3G==screen_cond_3G_new,1,0)) < G){
                    #print(screen_cond_3G)
                    #print(screen_cond_3G_new)
                    screen_cond_3G <- screen_cond_3G_new
                    #print(screen_cond_3G)
                    for (kk in 1:G){
                        screen_cond_3G_val = tg_tmp[[kk]]*0
                        for (m in 1:length(j)){
                            screen_cond_3G_val[,m] <- soft_thresh( tg_tmp[[kk]][,m] - ( t(as.matrix(W[j][[m]][,group_partition[[kk]]])) %*% (as.matrix(as.matrix(Xtilde[,j])[,m])*matrix(beta_check[j][m], N)) )/N, alpha*lambda)^2
                        }
                        screen_cond_3G[kk] <- (sqrt(sum(screen_cond_3G_val)) <= (1-alpha)*lambda*sqrt(df_X[l]*df_Z[kk])/sqrt(1+K))
                        #screen_cond_3G[kk] <- (sqrt(sum(screen_cond_3G_val)) <= (1-alpha)*lambda*sqrt(df_Z[kk])/sqrt(1+K))
                        #screen_cond_3G[kk] <- (sqrt(sum(soft_thresh( tg_tmp[[kk]] - ( t(as.matrix(W[[j]][,group_partition[[kk]]])) %*% (Xtilde[,j]*matrix(beta_check[j], N)) )/N, alpha*lambda)^2)) <= (1-alpha)*lambda*sqrt(df_Z[kk])/sqrt(1+K))
                    }
                    #print(screen_cond_3G)
                    zG <- sum(screen_cond_3G)
                    nzG <- G - zG
                    nzdf_Z <- df_Z[!screen_cond_3G]

                    if (nzG==0){
                        beta[j] <- beta_check[j]
                        theta[,j] <- theta[,j]*0
                    } else{
                        #print("start grad")
                        if (zG!=0){
                            z_group_partition <- group_partition[screen_cond_3G]
                            #print(z_group_partition)
                            for (kk in 1:zG){
                                theta[z_group_partition[[kk]],j] <- theta[z_group_partition[[kk]],j]*0
                            }
                        }

                        nz_group_partition <- group_partition[!screen_cond_3G]
                        t <- tt#/xsqmean[j]
                        c <- t*(1-alpha)*lambda*sqrt(df_X[l])
                        #c <- t*(1-alpha)*lambda
                        res <- Ytilde - (Xtilde %*% beta)
                        for (jj in 1:p){
                            res <- res - (as.matrix(W[[jj]]) %*% theta[,jj])
                        }
                        #grad_beta <- -(1/N) * ( solve(crossprod(Xtilde[,j])/N) %*% t(Xtilde[,j]) %*% res )
                        grad_beta <- -(1/N) * ( t(Xtilde[,j]) %*% res )
                        g1 <- sqrt(sum((beta[j] - t*grad_beta)^2))
                        #g1 <- abs(beta[j] - t*grad_beta)
                        grad_theta <- get_empty_list(paste0("grad_theta_",1:nzG))
                        g2 <- numeric(nzG)
                        for (i in 1:nzG){
                            grad_theta[[i]] <- matrix(0, nrow = length(nz_group_partition[[i]]), ncol = length(j))
                            g2_val <- numeric(length(j))
                            for (m in 1:length(j)){
                                grad_theta[[i]][,m] <- -(1/N) * ( t(as.matrix(W[j][[m]][,nz_group_partition[[i]]])) %*% res )
                                #print(soft_thresh(as.matrix(as.matrix(theta[,j])[,m])[nz_group_partition[[i]]] - t*grad_theta[[i]][,m], t*alpha*lambda)^2)
                                g2_val[m] <- sum(soft_thresh(as.matrix(as.matrix(theta[,j])[,m])[nz_group_partition[[i]]] - t*grad_theta[[i]][,m], t*alpha*lambda)^2)
                            }
                            #print(c("g2_val", g2_val))
                            g2[i] <- sqrt(sum(g2_val))
                            #grad_theta[[i]] <- -(1/N) * ( t(as.matrix(W[[j]][,nz_group_partition[[i]]])) %*% res )
                            #g2[i] <- sqrt(sum(soft_thresh(theta[,j][nz_group_partition[[i]]] - t*grad_theta[[i]], t*alpha*lambda)^2))
                        }
                        #r1 <- -c*sqrt(df_X[l])+sqrt((c^2)*df_X[l]-2*c*sum(g2*sqrt(nzdf_Z)/sqrt(1+K))+g1^2+sum(g2^2)-(df_X[l]-sum(nzdf_Z)/(1+K))*(c^2))
                        #r2 <- -c*sqrt(df_X[l])-sqrt((c^2)*df_X[l]-2*c*sum(g2*sqrt(nzdf_Z)/sqrt(1+K))+g1^2+sum(g2^2)-(df_X[l]-sum(nzdf_Z)/(1+K))*(c^2))
                        r1 <- -c+sqrt((c^2)-2*c*sum(g2*sqrt(nzdf_Z)/sqrt(1+K))+g1^2+sum(g2^2)-(1-sum(nzdf_Z)/(1+K))*(c^2))
                        r2 <- -c-sqrt((c^2)-2*c*sum(g2*sqrt(nzdf_Z)/sqrt(1+K))+g1^2+sum(g2^2)-(1-sum(nzdf_Z)/(1+K))*(c^2))
                        # a: norm of beta, b: norm of theta
                        # Hence, we choose the largest value of a and b to take positive value
                        a <- max(g1*r1/(c+r1), g1*r2/(c+r2), g1*r1/(c+r2), g1*r2/(c+r1))
                        #a <- max(g1*r1/(c*sqrt(df_X[l])+r1), g1*r2/(c*sqrt(df_X[l])+r2), g1*r1/(c*sqrt(df_X[l])+r2), g1*r2/(c*sqrt(df_X[l])+r1))
                        b <- numeric(nzG)
                        for (i in 1:nzG){
                            #b[i] <- max((g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r2/(c*sqrt(df_X[l])+r2), (g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r1/(c*sqrt(df_X[l])+r1), (g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r1/(c*sqrt(df_X[l])+r2), (g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r2/(c*sqrt(df_X[l])+r1))
                            b[i] <- max((g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r2/(c+r2), (g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r1/(c+r1), (g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r1/(c+r2), (g2[i]-c*sqrt(nzdf_Z[i])/sqrt(1+K))*r2/(c+r1))
                        }
                        c1 <- 1+t*sqrt(df_X[l])*(1-alpha)*lambda/sqrt(a^2+sum(b^2))
                        #c1 <- 1+t*(1-alpha)*lambda/sqrt(a^2+sum(b^2))
                        c2 <- numeric(nzG)
                        for (i in 1:nzG){
                            c2[i] <- 1+t*sqrt(df_X[l])*(1-alpha)*lambda*(sqrt(nzdf_Z[i])/sqrt(1+K)*1/b[i]+1/sqrt(a^2+sum(b^2)))
                            #c2[i] <- 1+t*(1-alpha)*lambda*(sqrt(nzdf_Z[i])/sqrt(1+K)*1/b[i]+1/sqrt(a^2+sum(b^2)))
                        }
                        beta[j] <- (beta[j] - t*grad_beta)/c1
                        for (i in 1:nzG){
                            theta[nz_group_partition[[i]],j] <- soft_thresh(theta[nz_group_partition[[i]],j] - t*grad_theta[[i]], t*alpha*lambda)/c2[i]
                            #theta[,j][nz_group_partition[[i]]] <- soft_thresh(theta[,j][nz_group_partition[[i]]] - t*grad_theta[[i]], t*alpha*lambda)/rep(c2[i], length(theta[,j][nz_group_partition[[i]]]))
                        }
                        #print(c("updated beta", beta[j]))
                        beta_check[j] <- beta[j]
                        for (kk in 1:G){
                            screen_cond_3G_val = tg_tmp[[kk]]*0
                            for (m in 1:length(j)){
                                screen_cond_3G_val[,m] <- soft_thresh( tg_tmp[[kk]][,m] - ( t(as.matrix(W[j][[m]][,group_partition[[kk]]])) %*% (as.matrix(as.matrix(Xtilde[,j])[,m])*matrix(beta_check[j][m], N)) )/N, alpha*lambda)^2
                            }
                            screen_cond_3G_new[kk] <- (sqrt(sum(screen_cond_3G_val)) <= (1-alpha)*lambda*sqrt(df_X[l]*df_Z[kk])/sqrt(1+K))
                            #screen_cond_3G_new[kk] <- (sqrt(sum(soft_thresh( tg_tmp[[kk]] - ( t(as.matrix(W[[j]][,group_partition[[kk]]])) %*% Xtilde[,j]*beta_check[j] )/N, alpha*lambda)^2)) <= (1-alpha)*lambda*sqrt(df_Z[kk])/sqrt(1+K))
                        }
                    }
                    #print(screen_cond_3G)
                    #print(screen_cond_3G_new)
                    #print(c("again?", (sum(ifelse(screen_cond_3G==screen_cond_3G_new,1,0)) < G), (sum(ifelse(screen_cond_3G==screen_cond_3G_new,1,0)))))
                }
            }
            #print(Xtilde[,j])
            #print(matrix(c(beta_old[j] - beta[j]), nrow = N, ncol = length(beta[j]), byrow = TRUE))
            full_res <- full_res + Xtilde[,j] %*% as.matrix(beta_old[j] - beta[j]) #+ W[[j]] %*% (theta_old[,j] - theta[,j])
            for (m in 1:length(j)){
                full_res <- full_res + W[j][[m]] %*% (as.matrix(theta_old[,j])[,m] - as.matrix(theta[,j])[,m])
            }
            #full_res <- Ytilde0 - (Xtilde %*% beta) # - beta0 - Ztilde %*% theta0
            #for (jj in 1:p){
            #    full_res <- full_res - (as.matrix(W[[jj]]) %*% theta[,jj])
            #}

            if (!is.null(Z)){
                lmmodel <- lm(full_res~Ztilde)
                if (zlinear == FALSE){
                    beta0 <- mean(full_res)
                    theta0 <- rep(0, K)
                } else {
                    beta0 <- lmmodel$coefficients[1]
                    theta0 <- lmmodel$coefficients[-1]
                }
                Ytilde <- Ytilde0 - beta0 - Ztilde %*% theta0
            }

            full_res2 <- full_res - beta0 - Ztilde %*% theta0
        }

        fmin=object_svReg(beta, theta, beta0, theta0, Xtilde, Ytilde0, Ztilde, W, main_partition, L, group_partition, G, alpha, lambda)
        error=abs(object_svReg(beta_old, theta_old, beta0_old, theta0_old, Xtilde, Ytilde0, Ztilde, W, main_partition, L, group_partition, G, alpha, lambda)-object_svReg(beta, theta, beta0, theta0, Xtilde, Ytilde0, Ztilde, W, main_partition, L, group_partition, G, alpha, lambda))
        #print(error)
        #print(beta[1:8])
        #print(theta)
        #if (sum(beta, theta)==0){itr <- iter}
    }
    print(c("iteration number: ",itr, "error: ", error))
    #print(beta[1:8])
    beta_raw <- beta*(1/SXYZ$Xweights)
    theta_raw <- theta*matrix((1/SXYZ$Xweights), K, p, byrow = TRUE)*matrix((1/SXYZ$Zweights), K, p, byrow = FALSE)

    WW <- get_empty_list(paste0("WW_",1:p))
    for (i in 1:p){
        WW[[i]] <- X[,i] * Z
    }
    #print(beta_raw)
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

    if (is.null(Z)){
        return(list("average_coef"=c(beta), "actual_coef"=list("main_coef"=beta, "modifying_coef"=theta), "raw_coef"=list("main_coef"=beta_raw, "modifying_coef"=theta_raw)))
    } else {
        return(list("average_coef"=c(beta, rowSums(theta)), "actual_coef"=list("main_coef"=beta, "modifying_coef"=theta), "raw_coef"=list("main_coef"=beta_raw, "modifying_coef"=theta_raw), "intercept"=list("beta0"=beta0, "theta0"=theta0), "intercept_raw"=list("beta0_raw"=beta0_raw, "theta0_raw"=theta0_raw), "fmin"=fmin))
    }
}


#' Fit the structural varying-coefficient regression (svReg) in the linear regression setting
#'
#' @param X N by p matrix of main predictors
#' @param Z N by K matrix of modifying variables. Modifying variables can take the form of continuous variables or categorical variables or both. Categorical variable should be coded by dummy variables (0-1).
#' @param Y vector of response variable
#' @param df_X vector of degrees of freedom for each group of main predictors.
#' @param df_Z vector of degrees of freedom for each group of modifying variables.
#' @param lambda_seq sequence of the tuning parameter, lambda. Can take the form of a sequence or a scalar.
#' @param alpha weight parameter between group penalty and individual penalty. Default value is 0.5.
#' @param tt learning rate for the gradient descent procedure. Default value is 0.1.
#' @param zlinear if true, the linear terms of the modifying variables are included. These terms are not regularized. Default value is TRUE.
#' @param tol tolerance for convergence. Convergence is determined by the value of the objective function: abs(objective_old - objective_new) is compared with the tolerance value. Default value is 1e-7.
#' @param iter maximum number of iteration for one lambda.  Default value is 500.
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
#' x = matrix(rnorm(100*5, 0, 1),100,5)
#' z1 = matrix(rnorm(100*3, 0, 1),100,3)
#' z2 = matrix(as.factor(sample(0:3, 100*2, prob=c(1/4,1/4,1/4,1/4), replace = TRUE)),100,2)
#' z2 = as.data.frame(model.matrix(~., data=as.data.frame(z2))[,-1])
#' z = cbind(z1, z2)
#' z = as.matrix(z)
#' y = 2*x[,1] - (2+2*z[,1])*x[,2] + (2+3*z[,4]+2*z[,5]-2*z[,6])*x[,3] + rnorm(100, 0, 1)
#' sv1=svReg(X=x,Z=z,Y=y,df_X=rep(1,5),df_Z=c(1,1,1,3,3),lambda_seq=c(1,0.5))
#' sv2=svReg(X=x,Z=z,Y=y,df_X=rep(1,5),df_Z=c(1,1,1,3,3),lambda_seq=0.5)
#' sv3=svReg(X=x,Z=z,Y=y,df_X=rep(1,5),df_Z=c(1,1,1,3,3),lambda_seq=0.5,zlinear=FALSE)
#' x[,3] = 2/3*x[,1] + 2/3*x[,2] + 1/3*rnorm(100, 0, 1)
#' y = x[,1] + x[,2] + (2+3*z[,4]+2*z[,5]-2*z[,6])*x[,4] + rnorm(100, 0, 1)
#' sv4=svReg(X=x,Z=z,Y=y,df_X=c(3,1,1),df_Z=c(1,1,1,3,3),lambda_seq=c(1,0.5))
svReg <- function(X, Z = NULL, Y, df_X, df_Z, lambda_seq = NULL, alpha = 0.5, tt = 0.1, zlinear = TRUE, tol = 1e-7, iter = 500){

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
    fit <- svReg1(X, Z, Y, df_X, df_Z, lambda_seq[1], alpha = alpha, tt = tt, beta = NULL, theta = NULL, zlinear = zlinear, tol = tol, iter = iter)
    para_array[,1,1] <- fit$actual_coef$main_coef; para_array[,-1,1] <- t(fit$actual_coef$modifying_coef)
    para_array_raw[,1,1] <- fit$raw_coef$main_coef; para_array_raw[,-1,1] <- t(fit$raw_coef$modifying_coef)
    beta0_vec[1] <- fit$intercept$beta0; theta0_vec[,1] <- fit$intercept$theta0
    beta0_raw_vec[1] <- fit$intercept_raw$beta0_raw; theta0_raw_vec[,1] <- fit$intercept_raw$theta0_raw
    fmin_vec[1] <- fit$fmin
    # Carry over previous beta for Warm Start
    if (length(lambda_seq) > 1){
        for (i in 2:length(lambda_seq)){
            print(c("lambda: ",lambda_seq[i]))
            fit <- svReg1(X, Z, Y, df_X, df_Z, lambda_seq[i], alpha = alpha, tt = tt, beta = para_array[,1,i-1], theta = t(para_array[,-1,i-1]), zlinear = zlinear, tol = tol, iter = iter)
            para_array[,1,i] <- fit$actual_coef$main_coef; para_array[,-1,i] <- t(fit$actual_coef$modifying_coef)
            para_array_raw[,1,i] <- fit$raw_coef$main_coef; para_array_raw[,-1,i] <- t(fit$raw_coef$modifying_coef)
            beta0_vec[i] <- fit$intercept$beta0; theta0_vec[,i] <- fit$intercept$theta0
            beta0_raw_vec[i] <- fit$intercept_raw$beta0_raw; theta0_raw_vec[,i] <- fit$intercept_raw$theta0_raw
            fmin_vec[i] <- fit$fmin
        }
    }

    return(list(lambda_seq = lambda_seq, beta_mat = para_array[,1,], theta_mat = para_array[,-1,], beta0_vec = beta0_vec, theta0_vec = theta0_vec, beta_raw_mat = para_array_raw[,1,], theta_raw_mat = para_array_raw[,-1,], beta0_raw_vec = beta0_raw_vec, theta0_raw_vec = theta0_raw_vec, fmin_vec = fmin_vec))

}


#' Perform k-fold cross validation for the structural varying-coefficient regression (svReg) model over a sequence of regularization parameter
#'
#' @param X N by p matrix of main predictors
#' @param Z N by K matrix of modifying variables. Modifying variables can take the form of continuous variables or categorical variables or both. Categorical variable should be coded by dummy variables (0-1).
#' @param Y vector of response variable
#' @param df_X vector of degrees of freedom for each group of main predictors.
#' @param df_Z vector of degrees of freedom for each group of modifying variables.
#' @param kfold the number of folds (=k) for the k-fold cross-validation. Default value is 10.
#' @param lambda_seq sequence of the tuning parameter, lambda. The length of the sequence should be greater than 1.
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
#' @return lambda_min: the lambda value which minimizes the continuous-categorical pliable lasso objective function among the values in the lambda_seq.
#' @return lambda_1se: the largest lambda value such that the difference with minimum objective function value is within 1 standard error of the minimum
#' @return cvm: the sequence of objective function values for lambda_seq
#' @return cvse: the sequence of the standard error of objective function values for lambda_seq
#' @return cvfold: (kfold) by (length of lambda_seq) matrix of the mean squared error of the test set for each fold
#' @return sqerror: N by (length of lambda_seq) matrix of the squared error
#' @export
#'
#' @examples
#' x=matrix(rnorm(100*5, 0, 1),100,5)
#' z1=matrix(rnorm(100*3, 0, 1),100,3)
#' z2=matrix(as.factor(sample(0:3, 100*2, prob=c(1/4,1/4,1/4,1/4), replace = TRUE)),100,2)
#' z2=as.data.frame(model.matrix(~., data=as.data.frame(z2))[,-1])
#' z=cbind(z1, z2)
#' z=as.matrix(z)
#' y=2*x[,1] - (2+2*z[,1])*x[,2] + (2+3*z[,4]+2*z[,5]-2*z[,6])*x[,3] + rnorm(100, 0, 1)
#' cv1=cv.svReg(X=x,Z=z,Y=y,df_X=rep(1,5),df_Z=c(1,1,1,3,3),lambda_seq=c(1,0.5))
#' cv2=cv.svReg(X=x,Z=z,Y=y,df_X=rep(1,5),df_Z=c(1,1,1,3,3),lambda_seq=c(1,0.5),cvseed=1)
#' cv3=cv.svReg(X=x,Z=z,Y=y,df_X=rep(1,5),df_Z=c(1,1,1,3,3),lambda_seq=0.5,zlinear=FALSE)
#' x[,3] = 2/3*x[,1] + 2/3*x[,2] + 1/3*rnorm(100, 0, 1)
#' y = x[,1] + x[,2] + (2+3*z[,4]+2*z[,5]-2*z[,6])*x[,4] + rnorm(100, 0, 1)
#' cv4=cv.svReg(X=x,Z=z,Y=y,df_X=c(3,1,1),df_Z=c(1,1,1,3,3),lambda_seq=c(1,0.5))
cv.svReg <- function(X, Z = NULL, Y, df_X, df_Z, kfold = 10, lambda_seq = NULL, alpha = 0.5, tt = 0.1, zlinear = TRUE, tol = 1e-7, iter = 500, cvseed = NULL){

    # p: number of main predictors; K: number of modifying variables, N: sample size
    p <- ncol(X)
    K <- ncol(Z)
    N <- nrow(X)

    # Fit Pliable Lasso on original data using plasso
    c2plassofit <- svReg(X, Z, Y, df_X, df_Z, lambda_seq = lambda_seq, alpha = alpha, tt = tt, zlinear = zlinear, tol = tol, iter = iter)

    # Split the data into K folds
    if (!is.null(cvseed)) set.seed(cvseed)
    idfold <- sample(1:N) %% kfold + 1

    # Calculate Pliable Lasso for each fold removed
    n_lambda <- length(c2plassofit$lambda_seq)
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
        cvfit <- svReg(xtrain, ztrain, ytrain, df_X, df_Z, lambda_seq, alpha = alpha, tt = tt, zlinear = zlinear, tol = tol, iter = iter)

        # Any additional calculations that are needed for calculating CV and SE_CV(lambda)
        testfitted <- matrix(rep(cvfit$beta0_vec, length(ytest)), nrow = length(ytest), byrow = TRUE) + ztest %*% cvfit$theta0_vec + xtest %*% cvfit$beta_mat
        for (j in 1:p){
            testfitted <- testfitted + (xtest[,j]*ztest) %*% cvfit$theta_mat[j,,]
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
    beta_raw_mat <- c2plassofit$beta_raw_mat
    beta0_raw_vec <- c2plassofit$beta0_raw_vec
    theta_raw_mat <- c2plassofit$theta_raw_mat
    theta0_raw_vec <- c2plassofit$theta0_raw_vec

    return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, theta_mat = theta_mat, theta0_vec = theta0_vec, beta_raw_mat = beta_raw_mat, beta0_raw_vec = beta0_raw_vec, theta_raw_mat = theta_raw_mat, theta0_raw_vec = theta0_raw_vec, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse, cvfold = cvfold, sqerror = sqerror))
}
