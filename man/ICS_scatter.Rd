% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICS_S3.R
\name{ICS_scatter}
\alias{ICS_scatter}
\alias{ICS_cov}
\alias{ICS_cov4}
\alias{ICS_covW}
\alias{ICS_covAxis}
\alias{ICS_tM}
\alias{ICS_scovq}
\title{Location and Scatter Estimates for ICS}
\usage{
ICS_cov(x, location = TRUE)

ICS_cov4(x, location = c("mean", "mean3", "none"))

ICS_covW(x, location = TRUE, alpha = 1, cf = 1)

ICS_covAxis(x, location = TRUE)

ICS_tM(x, location = TRUE, df = 1, ...)

ICS_scovq(x, y, ...)
}
\arguments{
\item{x}{a numeric matrix or data frame.}

\item{location}{for \code{ICS_cov()}, \code{ICS_cov4()}, \code{ICS_covW()},
and \code{ICS_covAxis()}, a logical indicating whether to include the sample
mean as location estimate (default to \code{TRUE}).  For \code{ICS_cov4()},
alternatively a character string specifying the location estimate can be
supplied.  Possible values are \code{"mean"} for the sample mean (the
default), \code{"mean3"} for a location estimate based on third moments,
or \code{"none"} to not include a location estimate.  For \code{ICS_tM()}
a logical inficating whether to include the M-estimate of location
(default to \code{TRUE}).}

\item{alpha}{parameter of the one-step M-estimator (default to 1).}

\item{cf}{consistency factor of the one-step M-estimator (default to 1).}

\item{df}{assumed degrees of freedom of the t-distribution (default to 1,
which corresponds to the Cauchy distribution).}

\item{\dots}{additional arguments to be passed down to \code{\link{scovq}()}.}

\item{y}{numerical vector specifying the dependent variable.}
}
\value{
An object of class \code{"ICS_scatter"} with the following
components:
\item{location}{if requested, a numeric vector giving the location
estimate.}
\item{scatter}{a numeric matrix giving the estimate of the scatter matrix.}
\item{label}{a character string providing a label for the scatter matrix.}
}
\description{
Computes a scatter matrix and an optional location vector to be used in
transforming the data to an invariant coordinate system or independent
components.
}
\details{
\code{ICS_cov()} is a wrapper for the sample covariance matrix as computed
by \code{\link[stats]{cov}()}.

\code{ICS_cov4()} is a wrapper for the scatter matrix based on fourth
moments as computed by \code{\link{cov4}()}. Note that the scatter matrix
is always computed with respect to the sample mean, even though the returned
location component can be specified to be based on third moments as computed
by \code{\link{mean3}()}.  Setting a location component other than the
sample mean can be used to fix the signs of the invariant coordinates in
\code{\link{ICS}()} based on generalized skewness values, for instance
when using the scatter pair \code{ICS_cov()} and \code{ICS_cov4()}.

\code{ICS_covW()} is a wrapper for the one-step M-estimator of scatter as
computed by \code{\link{covW}()}.

\code{ICS_covAxis()} is a wrapper for the one-step Tyler shape matrix as
computed by \code{\link{covAxis}()}, which is can be used to perform
Principal Axis Analysis.

\code{ICS_tM()} is a wrapper for the M-estimator of location and scatter
for a multivariate t-distribution, as computed by \code{\link{tM}()}.

\code{ICS_scovq()} is a wrapper for the supervised scatter matrix based
on quantiles scatter, as computed by \code{\link{scovq}()}.
}
\examples{
data("iris")
X <- iris[,-5]
ICS_cov(X)
ICS_cov4(X)
ICS_covW(X, alpha = 1, cf = 1/(ncol(X)+2))
ICS_covAxis(X)
ICS_tM(X)


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
ICS_scovq(X, y = y)

}
\references{
Arslan, O., Constable, P.D.L. and Kent, J.T. (1995) Convergence behaviour of
the EM algorithm for the multivariate t-distribution, \emph{Communications
in Statistics, Theory and Methods}, \bold{24}(12), 2981--3000.
\doi{10.1080/03610929508831664}.

Critchley, F., Pires, A. and Amado, C. (2006) Principal Axis Analysis.
Technical Report, \bold{06/14}. The Open University, Milton Keynes.

Kent, J.T., Tyler, D.E. and Vardi, Y. (1994) A curious likelihood identity
for the multivariate t-distribution, \emph{Communications in Statistics,
Simulation and Computation}, \bold{23}(2), 441--453.
\doi{10.1080/03610919408813180}.

Oja, H., Sirkia, S. and Eriksson, J. (2006) Scatter Matrices and Independent
Component Analysis. \emph{Austrian Journal of Statistics}, \bold{35}(2&3),
175-189.

Tyler, D.E., Critchley, F., Duembgen, L. and Oja, H. (2009) Invariant
Co-ordinate Selection. \emph{Journal of the Royal Statistical Society,
Series B}, \bold{71}(3), 549--592. \doi{10.1111/j.1467-9868.2009.00706.x}.
}
\seealso{
\code{\link{ICS}()}

\code{\link[base:colSums]{colMeans}()}, \code{\link{mean3}()}

\code{\link[stats]{cov}()}, \code{\link{cov4}()}, \code{\link{covW}()},
\code{\link{covAxis}()}, \code{\link{tM}()}, \code{\link{scovq}()}
}
\author{
Andreas Alfons and Aurore Archimbaud
}
