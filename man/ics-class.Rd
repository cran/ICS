\name{ics-class}
\docType{class}
\alias{ics-class}

\title{Class ICS}
\description{A S4 class to store results from an invariant coordinate system transformation or independent component computation based on two scatter matrices.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ics", ...)}. But usual objects are created by the function \code{\link{ics}}.
     
}
\section{Slots}{
     \describe{
    \item{\code{Kurt}:}{Object of class \code{"numeric"}. Gives the kurtosis measures of the components}
    \item{\code{UnMix}:}{Object of class \code{"matrix"}. The unmixing matrix.}
    \item{\code{S1}:}{Object of class \code{"character"}. Name of the first scatter matrix.}
    \item{\code{S2}:}{Object of class \code{"character"}. Name of the second scatter matrix.}
    \item{\code{Scores}:}{Object of class \code{"data.frame"}. The underlying components in the invariant coordinate system.}
    \item{\code{DataNames}:}{Object of class \code{"character"}. Names of the original variables.}
    \item{\code{Standardize}:}{Object of class \code{"character"}. Name standardization method for UnMix.}
  }
}
\section{Methods}{
For this class the following generic functions are available: \code{\link{print.ics}}, \code{\link{summary.ics}}, \code{\link{coef.ics}}, \code{\link{fitted.ics}} and \code{\link{plot.ics}}
}

\section{Note}{
In case no extractor function for the slots exists, the component can be extracted the usual way using '@'.
}

\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}



\seealso{
\code{\link{ics}}
}

\keyword{classes}
