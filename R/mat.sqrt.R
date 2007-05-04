`mat.sqrt` <-
function(A) 
    {
    eigen.A<-eigen(A)
    sqrt.A<-eigen.A$vectors%*%(diag(eigen.A$values))^0.5%*%t(eigen.A$vectors)
    return(sqrt.A)
    }

