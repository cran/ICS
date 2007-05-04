# covariance matrix based on 4th moments wrt to the mean vector

.cov4moments.mean<-function(X)
    {
    p<-dim(X)[2]
    n<-dim(X)[1]
    data.centered<-sweep(X,2,colMeans(X),"-")
    Sigma.data.sqrt<-mat.sqrt(cov(X)) 
    radius<-sqrt(rowSums((data.centered%*%solve(Sigma.data.sqrt))^2))
    V<-((n+1)/(n*(p+2)))*cov(radius*data.centered)  
    return(V) 
    }

# covariance matrix based on 4th moments wrt to origin

.cov4moments.origin<-function(X)
    {
    p<-dim(X)[2]
    n<-dim(X)[1]
    Sigma.data.sqrt<-mat.sqrt(covOrigin(X)) 
    radius<-sqrt(rowSums((X%*%solve(Sigma.data.sqrt))^2))
    V<-(1/(p+2))*covOrigin(radius*X)  
    return(V) 
    }
    
# Iteration step for Tyler's shape matrix
    
.tyler.step<-function(V.old,datas,p,n)
        {
        sqrt.V.old<-mat.sqrt(V.old)
        r<-sqrt(rowSums((datas %*% sqrt.V.old)^2))
        M.V.old<-p/n*(t(((1/r)*datas%*%sqrt.V.old))%*%((1/r)*datas%*%sqrt.V.old))
        V.new<-sum(diag(sqrt.V.old %*% solve(M.V.old)))^(-1)*(sqrt.V.old %*% solve(M.V.old) %*% sqrt.V.old)
        return(V.new)
        }

# Sign of the maximum element of a row in a matrix
      
.sign.max<-function(x)
 {
 ifelse(identical(max(x),max(abs(x))),1,-1)
 }
        
