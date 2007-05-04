# defining the ics S4-class

setClass("ics",representation(Kurt="numeric",UnMix="matrix",S1="character",S2="character",Scores="data.frame",
DataNames="character",Standardize="character"))

setValidity("ics",function(object){
    if(!is(object@Kurt, "numeric")) return("Kurtosis values of ics objects must be numeric")
    if(!(is(object@UnMix, "matrix") )) return("slot 'UnMix' of a ics object must be a numeric matrix")
    if(!(is.numeric(object@UnMix) )) return("slot 'UnMix' of a ics object must be a numeric matrix")
    if(!is(object@S1, "character") || length(object@S1)!=1) return("slot 'S1' of a ics object must be the name of a scatter function")

    if(!is(object@S2, "character") || length(object@S2)!=1) return("slot 'S2' of a ics object must be the name of a scatter function")
    if(!(is(object@Scores, "data.frame"))) return("slot 'Scores' of a ics object must be a numeric data frame")
    if(!all(sapply(object@Scores, is.numeric))) return("slot 'Scores' of a ics object must be a numeric data frame")
    if(!is(object@DataNames, "character")) return("slot 'DataNames' of a ics object must give the column names of the data matrix")
    if(!is(object@Standardize, "character") || length(object@Standardize)!=1) return("slot 'Standardize' of a ics object must be the name of a standardization method of 'UnMix'")

    if(length(object@Kurt)!=dim(object@UnMix)[2]) return("length of 'Kurt' must correspond to the number of columns of 'UnMix'")
    if(length(object@Kurt)!=dim(object@Scores)[2]) return("length of 'Kurt' must correspond to the number of columns of 'Scores'")
    if(length(object@Kurt)!=length(object@DataNames)) return("length of 'Kurt' must be the same as length of 'DataNames'")
    if(length(object@Kurt)<2) return("at least bivariate data is needed")

    return(TRUE)
})

setMethod("show",signature(object="ics"),
function(object)
    {
    tmp <- list(Kurt=object@Kurt,
                UnMix=object@UnMix)
    print(tmp,quote=F)
    invisible(tmp)  
    }
)




setMethod("plot",signature(x="ics",y="missing"),
function(x,index=NULL,...)
    {
    p<-ncol(x@UnMix)
    if (is.null(index) & p<=6) pairs(x@Scores,...)
    if (is.null(index) & p>6) pairs(x@Scores[,c(1:3,p-2:0)],...)
    if (length(index)==1) stop("index must be NULL or at least a vector of length 2")
    if (length(index)>1) pairs(x@Scores[,index],...)
    }
)

setMethod("summary",signature(object="ics"),
function(object,digits=4)
    {   
    cat("\nICS based on two scatter matrices \n")
    cat("S1: ", object@S1) 
    cat("\nS2: ",object@S2)
    cat("\n")
    cat("\nThe Kurtosis measures of the components are:\n")
    print(format(round(object@Kurt,digits)),quote=F)
    cat("\n")
    cat("\nThe Unmixing matrix is:\n")
    print(round(object@UnMix,digits))
    invisible(object)
    }
)


setMethod("fitted",signature(object="ics"),
function(object, index=NULL)
    {
    p<-ncol(object@Scores)
    if (is.null(index)==FALSE && max(index)>p) stop("undefined columns selected")
    Mix<-solve(object@UnMix)
    if (is.null(index)) index=1:p
    fits<-as.matrix(object@Scores[,index])%*%t(Mix[,index])
    return(as.data.frame(fits))
    }
)

setMethod("coef",signature(object="ics"),
function(object)
    {
    return(object@UnMix)
    }
)
