\name{ICS-package}
\alias{ICS-package}

\docType{package}
\title{
Tools for Exploring Multivariate Data via ICS/ICA
}
\description{
This package implements the 2 scatter matrix transformation to obtain an invariant coordinate sytem or independent components, depending on the underlying assumptions. 
The result of the transformation is an object of the S4 class \code{ics} which is provided by this package.
Besides generic functions to create and work with an \code{ics} object the package contains also some scatter matrices. 
}
\details{
\tabular{ll}{
Package: \tab ICS\cr
Type: \tab Package\cr
Version: \tab 1.1-1\cr
Date: \tab 2008-02-18\cr
License: \tab GPL (>= 2)\cr
}
Some multivariate tests and estimates are not affine equivariant by nature. A possible remedy for the lack of that property is to transform the data points to an invariant coordinate system,
construct tests and estimates from the transformed data, and if needed, retransform the estimates back. The use of two different
scatter matrices to obtain invariant coordinates is implemeted in this package by the function \code{ics}. For an invariant coordinate selection no 
assumptions are made about the data or the scatter matrices and it can be seen as a data transformation method. If the data come, however, from a so called independent component model 
the \code{ics} function can recover the independent components and estimate the mixing matrix under general assumptions. 
Besides the function \code{ics}, which creates an S4 object of class \code{ics}, provides this package tools to work with objects of that class and some
scatter matrices which can be used in the \code{ics} function. Furthermore are there also two tests for multinormality.  
}
\author{
Klaus Nordhausen, Hannu Oja and David E. Tyler

Maintainer: 
Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}
}

\keyword{package}
