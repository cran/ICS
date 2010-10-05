### function to compute a symmetric square root of a matrix
### internal function
###

`mat.sqrt` <-
function(A) 
    {
    eigen.A<-eigen(A)
    sqrt.A<-eigen.A$vectors%*%(diag(eigen.A$values))^0.5%*%t(eigen.A$vectors)
    return(sqrt.A)
    }
