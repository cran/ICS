`covAxis` <-
function(X, na.action = na.fail)
    {
    X<-na.action(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    
    n <- dim(X)[1]
    p <- dim(X)[2]                                                
    if (p<2) stop("'X' must be at least bivariate")  
    
    Xmeans <- colMeans(X)
    di<-sqrt(mahalanobis(X,Xmeans,cov(X)))
    X.centered <- sweep(X, 2, Xmeans)
    Y<-sweep(X.centered,1,di,FUN="/")

    v.tilde <- p*t(Y)%*%Y / n
    return(v.tilde)
    }
