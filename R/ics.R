`ics` <-
function(X, S1=cov, S2=cov4, S1args = list(), S2args = list(), standardize = "Z", na.action = na.fail)
    {
    X<-na.action(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    data.matrix<-as.matrix(X)
    
    if (is.numeric(data.matrix)==FALSE ) stop("'X' must be numeric")
    
    if (standardize != "B" & standardize != "Z") stop("'standardize' must be 'B' or 'Z'")
    
    p <- dim(X)[2]                                                                                                
    if (p<2) stop("'X' must be at least bivariate")  
    
    
    B1 <- solve(mat.sqrt(do.call("S1", c(list(X), S1args))))
    X1 <- data.matrix %*% B1
    B2 <- do.call("S2", c(list(X1), S2args))
    B2.eigen <- eigen(B2)
    U2 <- B2.eigen$vectors
    DiagB2 <- B2.eigen$values
    X2 <- X1 %*% U2
    B <- t(U2) %*% B1
    
    if (standardize=="B")
        {
        row.signs <- apply(B, 1, .sign.max)
        row.norms <- sqrt(rowSums((B)^2))
        B.res <- sweep(B, 1, row.norms*row.signs, "/")
        Z <- as.data.frame(data.matrix %*% t(B.res))
        }
    if (standardize=="Z")
        {
        Z1 <- data.matrix %*% t(B)
        skewness <- colMeans(Z1) - apply(Z1, 2, median)
        skew.signs <- ifelse(skewness > 0, 1, -1)
        B.res <- sweep(B, 1, skew.signs, "*")
        Z <- as.data.frame(data.matrix %*% t(B.res))
        }
        
    names(Z)<-paste(rep("IC",p),1:p,sep=".")
    
    if(is.null(colnames(X))==TRUE) 
        names.X<-paste(rep("X",p),1:p,sep=".") 
    else 
        names.X<-colnames(X)
    
    res <- new("ics", Kurt = DiagB2, UnMix = B.res,
               S1 = deparse(substitute(S1)),
               S2 = deparse(substitute(S2)),
               Scores = Z,
               DataNames = names.X,
               Standardize = standardize)
    return(res)
}
