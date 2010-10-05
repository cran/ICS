### Location estimate based on 3rd moments
###
###

mean3 <- function(X, na.action = na.fail)
    {
    X<-na.action(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    p <- dim(X)[2]
    if (p<2) stop("'X' must be at least bivariate")
    n <- dim(X)[1]
    cor.fac <- (n-1)/n
    r.i.2 <- mahalanobis(X, colMeans(X), cor.fac*cov(X))
    Y <- sweep(X,1,r.i.2,"*")
    colMeans(Y)/p
    }
