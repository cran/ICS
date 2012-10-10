\name{scovq}
\alias{scovq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Supervised scatter matrix based on quantiles
}
\description{
Function for a supervised scatter matrix that is the weighted
covariance matrix of \code{x} with weights 1/(\code{q2-q1}) if \code{y} is between the
lower (\code{q1}) and upper (\code{q2}) quantile and 0 otherwise (or vice versa). 
}

\usage{
scovq(x, y, q1 = 0, q2 = 0.5, pos = TRUE, type = 7, 
      method = "unbiased", na.action = na.fail, 
      check = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric data matrix with at least two columns.
}
  \item{y}{
numerical vector specifying the dependent variable.
}
  \item{q1}{
percentage for lower quantile of \code{y}. With 0 <= \code{q1} < \code{q2}. See details.
}
  \item{q2}{
percentage for upper quantile of \code{y}. With \code{q1} < \code{q2} <= 1. See details.
}
  \item{pos}{
logical. If TRUE then the weights are 1/(\code{q2-q1}) if \code{y} is between the \code{q1}- and \code{q2}-
     quantiles and 0 othervise. If FALSE then the weights are 0 if \code{y}
     between \code{q1}- and \code{q2}-quantiles  and 1/(\code{1-q2+q1}) otherwise.
}
  \item{type}{
passed on to function \code{\link{quantile}}. 
}
  \item{method}{
passed on to function \code{\link{cov.wt}}. 
}
  \item{na.action}{
a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.
}
  \item{check}{
logical. Checks if the input should be checked for consistency. If not needed setting
it to FALSE might save some time.
}
}
\details{

The weights for this supervised scatter matrix for \code{pos=TRUE} are
\eqn{w(y) = I(q1-quantile < y < q2-quantile)/(q2-q1)}. Then \code{scovq} is calculated as
\deqn{scovq = \sum w(y) (x-\bar{x}_w)'(x-\bar{x}_w).}{%
  scovq = \sum w(y) (x-x_w_bar)'(x-x_w_bar).}  
where \eqn{\bar{x}_w = \sum w(y) x}{ x_w_bar = sum w(y)x}.

To see how this function can be used in the context of supervised invariant coordinate selection
see the example below.
}
\value{
a matrix.
}
\references{
\cite{Liski, E., Nordhausen, K. and Oja, H. (2013), Supervised invariant coordinate selection, \emph{To appear in Statistics}, \bold{???},  ???--???.}
}
\author{Klaus Nordhausen}
\seealso{ \code{\link{cov.wt}} and \code{\link{ics}}}
\examples{
# Creating some data

# The number of explaining variables
p <- 10
# The number of observations
n <- 400
# The error variance
sigma <- 0.5
# The explaining variables 
X <- matrix(rnorm(p*n),n,p)
# The error term
epsilon <- rnorm(n, sd = sigma)
# The response
y <- X[,1]^2 + X[,2]^2*epsilon


# SICS with ics

X.centered <- sweep(X,2,colMeans(X),"-")
SICS <- ics(X.centered, S1=cov, S2=scovq, S2args=list(y=y, q1=0.25, 
        q2=0.75, pos=FALSE), stdKurt=FALSE, stdB="Z")

# Assuming it is known that k=2, then the two directions 
# of interest are choosen as:

k <- 2
KURTS <- SICS@gKurt 
KURTS.max <- ifelse(KURTS >= 1, KURTS, 1/KURTS)
ordKM <- order(KURTS.max, decreasing = TRUE)

indKM <- ordKM[1:k]

# The two variables of interest
Zk <- ics.components(SICS)[,indKM]

# The correspondings transformation matrix
Bk <- coef(SICS)[indKM,]

# The corresponding projection matrix
Pk <- t(Bk) \%*\% solve(Bk \%*\% t(Bk)) \%*\% Bk

# Visualization
pairs(cbind(y,Zk))

# checking the subspace difference

# true projection

B0 <- rbind(rep(c(1,0),c(1,p-1)),rep(c(0,1,0),c(1,1,p-2)))
P0 <- t(B0) \%*\% solve(B0 \%*\% t(B0)) \%*\% B0

# crone and crosby subspace distance measure, should be small
k - sum(diag(P0 \%*\% Pk))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
