mean3 <- function(X, na.action = na.fail)
    {
    X<-na.action(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    p <- dim(X)[2]
    if (p<2) stop("'X' must be at least bivariate")
    r.i.2 <-  mahalanobis(X, colMeans(X), cov(X))
    Y <- sweep(X,1,r.i.2,"*")
    colMeans(Y)/p
    }
