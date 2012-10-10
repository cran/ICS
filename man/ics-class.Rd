\name{ics-class}
\docType{class}
\alias{ics-class}

\title{Class ICS}
\description{A S4 class to store results from an invariant coordinate system transformation or independent component computation based on two scatter matrices.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ics", ...)}. But usually objects are created by the function \code{\link{ics}}.
     
}
\section{Slots}{
     \describe{
    \item{\code{gKurt}:}{Object of class \code{"numeric"}. Gives the generalized kurtosis measures of the components}
    \item{\code{UnMix}:}{Object of class \code{"matrix"}. The unmixing matrix.}
    \item{\code{S1}:}{Object of class \code{"matrix"}. The first scatter matrix.}
    \item{\code{S2}:}{Object of class \code{"matrix"}. The second scatter matrix.}
    \item{\code{S1name}:}{Object of class \code{"character"}. Name of the first scatter matrix.}
    \item{\code{S2name}:}{Object of class \code{"character"}. Name of the second scatter matrix.}
    \item{\code{Scores}:}{Object of class \code{"data.frame"}. The underlying components in the invariant coordinate system.}
    \item{\code{DataNames}:}{Object of class \code{"character"}. Names of the original variables.}
    \item{\code{StandardizeB}:}{Object of class \code{"character"}. Names standardization method for UnMix.}
    \item{\code{StandardizegKurt}:}{Object of class \code{"logical"}. States wether the generalized kurtosis is standardized or not.}
  }
}
\section{Methods}{
For this class the following generic functions are available: \code{\link{print.ics}}, \code{\link{summary.ics}}, \code{\link{coef.ics}}, \code{\link{fitted.ics}} and \code{\link{plot.ics}}
}

\section{Note}{
In case no extractor function for the slots exists, the component can be extracted the usual way using '@'.
}

\author{Klaus Nordhausen}



\seealso{
\code{\link{ics}}
}

\keyword{classes}
