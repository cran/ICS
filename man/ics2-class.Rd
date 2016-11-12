\name{ics2-class}
\docType{class}
\alias{ics2-class}

\title{Class ICS2}
\description{A S4 class to store results from an invariant coordinate system transformation or independent component computation based on two scatter matrices and two location vectors.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ics2", ...)}. But usually objects are created by the function \code{\link{ics2}}.
The Class inherits from the \code{ics} class.
     
}
\section{Slots}{
     \describe{
    \item{\code{gSkew}:}{Object of class \code{"numeric"}. Gives the generalized skewness measures of the components}   
    \item{\code{gKurt}:}{Object of class \code{"numeric"}. Gives the generalized kurtosis measures of the components}
    \item{\code{UnMix}:}{Object of class \code{"matrix"}. The unmixing matrix.}
    \item{\code{S1}:}{Object of class \code{"matrix"}. The first scatter matrix.}
    \item{\code{S2}:}{Object of class \code{"matrix"}. The second scatter matrix.}
    \item{\code{T1}:}{Object of class \code{"numeric"}. The first location vector.}
    \item{\code{T2}:}{Object of class \code{"numeric"}. The second location vector.}
    \item{\code{S1name}:}{Object of class \code{"character"}. Name of the first scatter matrix.}
    \item{\code{S2name}:}{Object of class \code{"character"}. Name of the second scatter matrix.}
    \item{\code{S1args}:}{Object of class \code{"list"}. Additional arguments needed when calling function S1.}
    \item{\code{S2args}:}{Object of class \code{"list"}. Additional arguments needed when calling function S2.}
    \item{\code{Scores}:}{Object of class \code{"data.frame"}. The underlying components in the invariant coordinate system.}
    \item{\code{DataNames}:}{Object of class \code{"character"}. Names of the original variables.}
    \item{\code{StandardizeB}:}{Object of class \code{"character"}. Names standardization method for UnMix.}
    \item{\code{StandardizegKurt}:}{Object of class \code{"logical"}. States wether the generalized kurtosis is standardized or not.}
  }
}
\section{Methods}{
For this class the following generic functions are available: \code{\link{print.ics2}}, \code{\link{summary.ics2}}
But naturally the other methods like plot, coef, fitted and so from class ics work via inheritance.
}

\section{Note}{
In case no extractor function for the slots exists, the component can be extracted the usual way using '@'.
}

\author{Klaus Nordhausen}



\seealso{
\code{\link{ics2}}
}

\keyword{classes}
