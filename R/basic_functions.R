# Center the response variable, scale and center the main predictors and modifying variables
standardizeXYZ <- function(X, Z, Y){

    n = length(Y)
    p = ncol(X)
    K = ncol(Z)

    # Center Y
    Ymean <- mean(Y)
    Ytilde <- Y - Ymean

    # Center and scale X and Z
    Xmeans <- colMeans(X)
    Xcentered <- X - matrix(Xmeans, n, p, byrow = TRUE)
    Xweights <- sqrt(colSums(Xcentered^2)/n)
    Xweights[Xweights==0] <- 1
    Xtilde <- as.matrix(Xcentered) %*% diag(1/Xweights, nrow = p, ncol = p)
    #Xtilde <- scale(X)*sqrt(n/(n-1))

    # Data values for all participants may be zero
    Xtilde[Xtilde=="NaN"] <- 0

    Zmeans <- colMeans(Z)
    Zcentered <- Z - matrix(Zmeans, n, K, byrow = TRUE)
    Zweights <- sqrt(colSums(Zcentered^2)/n)
    Zweights[Zweights==0] <- 1
    if (ncol(Zcentered)>1){
        Ztilde <- as.matrix(Zcentered) %*% diag(1/Zweights)
    } else{
        Ztilde <- as.matrix(Zcentered) * 1/Zweights
    }
    #Ztilde <- scale(Z)*sqrt(n/(n-1))
    Ztilde[Ztilde=="NaN"] <- 0

    # Return the mean of Y and means of columns of X, as well as weights to be used in back-scaling (that is sqrt(X_j'X_j/n))
    return(list(Ytilde = Ytilde, Xtilde = Xtilde, Ztilde = Ztilde, Ymean = Ymean, Xmeans = Xmeans, Zmeans = Zmeans, Xweights = Xweights, Zweights = Zweights))

}

# Vector soft-thresholding operator
soft_thresh <- function(a, lambda){
    return(as.numeric(sign(a)*pmax(0,abs(a) - lambda)))
}

# Compute plasso objective function value
object_plasso <- function(bt, th, bt0, th0, Xtilde, Ytilde, Ztilde, W, alpha, lambda){

    p <- ncol(Xtilde)
    N <- nrow(Xtilde)

    diff <- Ytilde - Xtilde %*% bt - bt0 - Ztilde %*% th0
    for (i in 1:p){
        diff <- diff - W[[i]] %*% th[,i]
    }
    obj <- (1/(2*N))*sum((diff)**2)
    for (i in 1:p){
        obj <- obj + (1-alpha)*lambda*(sqrt(bt[i]^2+sum(th[,i]^2)) + sqrt(sum(th[,i]^2)))
    }
    obj <- obj + alpha*lambda*sum(abs(th))
    return(obj)
}

# Compute c2plasso objective function value
object_c2plasso <- function(bt, th, bt0, th0, Xtilde, Ytilde, Ztilde, W, group_partition, G, alpha, lambda){

    p <- ncol(Xtilde)
    N <- nrow(Xtilde)
    K <- ncol(Ztilde)

    diff <- Ytilde - Xtilde %*% bt - bt0 - Ztilde %*% th0
    for (i in 1:p){
        diff <- diff - W[[i]] %*% th[,i]
    }
    obj <- (1/(2*N))*sum((diff)**2)
    for (i in 1:p){
        obj <- obj + (1-alpha)*lambda*(sqrt(bt[i]^2+sum(th[,i]^2)))
        for (kk in 1:G){
            obj <- obj + (1-alpha)*lambda*sqrt(length(group_partition[[kk]])/(1+K))*sqrt(sum(th[,i][group_partition[[kk]]]^2))
        }
    }
    obj <- obj + alpha*lambda*sum(abs(th))
    return(obj)
}

# Compute svReg objective function value
object_svReg <- function(bt, th, bt0, th0, Xtilde, Ytilde, Ztilde, W, main_partition, L, group_partition, G, alpha, lambda){

    p <- ncol(Xtilde)
    N <- nrow(Xtilde)
    K <- ncol(Ztilde)

    diff <- Ytilde - Xtilde %*% bt - bt0 - Ztilde %*% th0
    for (i in 1:p){
        diff <- diff - W[[i]] %*% th[,i]
    }
    obj <- (1/(2*N))*sum((diff)**2)
    for (l in 1:L){
        main_group <- main_partition[[l]]
        i <- main_group
        #obj <- obj + (1-alpha)*lambda*(sqrt(bt[i]^2+sum(th[,i]^2)))
        obj <- obj + (1-alpha)*lambda*sqrt(length(main_partition[[l]]))*(sqrt(sum(bt[i]^2)+sum(th[,i]^2)))
        for (kk in 1:G){
            obj <- obj + (1-alpha)*lambda*sqrt(length(main_partition[[l]]))*sqrt(length(group_partition[[kk]])/(1+K))*sqrt(sum(th[,i][group_partition[[kk]]]^2))
        }
    }
    obj <- obj + alpha*lambda*sum(abs(th))
    return(obj)
}


# make empty lists for storage
get_empty_list <- function(names){
    out <-  vector("list",length(names))
    names(out) <- names
    return(out)
}

# Report an error if x <=0
error_positive_value <- function(x){
    return(invisible(if(any(x <=0)){
        stop("The number of elements of each modifying variable group must be >0.")
    }
    ))
}

# Check the sum of the degrees of freedom of modifying variables
error_modifying_variable_df <- function(x,y){
    return(invisible((
        if(x != y){
            stop("The sum of the degrees of freedom of Z should be equal to the numbers of columns of Z")
        })))
}
