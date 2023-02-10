### screeplot method for an ics object
### barplot or lineplot are options
###
#' @export
#' @importFrom graphics plot barplot axis
`screeplot.ics` <-
    function(x,index=NULL,type="barplot",
               main = deparse(substitute(x)),ylab="generalized kurtosis",xlab= "component", names.arg=index, labels=TRUE,...)
    {
    #if(class(x)!="ics") stop("'x' must be of class ics")
    type<-match.arg(type,c("barplot","lines"))
    if (is.null(index)) index=1:length(x@gKurt)
    if (type=="barplot")
        {
        barplot(x@gKurt[index],names.arg=names.arg,ylab=ylab,xlab=xlab,main=main,...)
        }
    else
        {
        plot(index,x@gKurt[index],type="b",ylab=ylab,xlab=xlab,axes=FALSE,main=main,...)
        axis(2)
        axis(1, at = seq_along(index), labels = labels)
        }
    invisible()
    }
