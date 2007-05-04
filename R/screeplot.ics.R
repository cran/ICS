`screeplot.ics` <-
    function(x,index=NULL,type="barplot",
               main = deparse(substitute(x)),ylab="kurtosis measure",xlab= "component",...)
    {
    if(class(x)!="ics") stop("'x' must be of class ics")
    type<-match.arg(type,c("barplot","lines"))
    if (is.null(index)) index=1:length(x@Kurt)
    if (type=="barplot")
        {
        barplot(x@Kurt[index],names.arg=index,ylab=ylab,xlab=xlab,main=main,...)
        }
    else
        {
        plot(index,x@Kurt[index],type="b",ylab=ylab,xlab=xlab,axes=F,main=main,...)
        axis(2)
        axis(1, at = index, labels = T)
        }
    invisible()
    }
