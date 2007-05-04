`frobenius.norm` <-
function(A)
    {
    sqrt(sum(diag(t(A)%*%A)))
    }

