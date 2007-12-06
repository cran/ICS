\name{screeplot.ics}
\alias{screeplot.ics}

\title{Screeplot for an ICS Object }
\description{
  Plots the kurtosis measures of an \code{ics} object against its index number. Two versions of this screeplot are available.
}
\usage{
\method{screeplot}{ics}(x, index = NULL, type = "barplot", 
          main = deparse(substitute(x)), ylab = "generalized kurtosis", 
          xlab = "component", ...)
}

\arguments{
  \item{x}{ object of class \code{ics}.}
  \item{index}{index of the components to be plottes. If NULL all components are used.}
  \item{type}{ \code{barplot} if a barplot or \code{lines} if a line plot is preferred.}
  \item{main}{ main title of the plot. }
  \item{ylab}{ y-axis label. }
  \item{xlab}{ x-axis label. }
  \item{\dots}{other arguments for the plotting functions. }
}


\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}

\seealso{\code{\link{plot.ics}}, \code{\link{ics-class}} and \code{\link{ics}}}
\examples{
set.seed(654321)
A <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2),ncol=3)
eigen.A <- eigen(A)
sqrt.A <- eigen.A$vectors \%*\% (diag(eigen.A$values))^0.5 \%*\% t(eigen.A$vectors)
normal.ic <- cbind(rnorm(800), rnorm(800), rnorm(800))
mix.ic <- cbind(rt(800,4), rnorm(800), runif(800,-2,2))

data.normal <- normal.ic \%*\% t(sqrt.A)
data.mix <- mix.ic \%*\% t(sqrt.A)

par(mfrow=c(1,2))
screeplot(ics(data.normal))
screeplot(ics(data.mix), type="lines")
par(mfrow=c(1,1))
rm(.Random.seed)
}

\keyword{ hplot }

