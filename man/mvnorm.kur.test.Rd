\name{mvnorm.kur.test}
\alias{mvnorm.kur.test}
\title{Test of Multivariate Normality Based on Kurtosis}
\description{
  Test for multivariate normality which uses as criterion the kurtosis measured by the ratio of regular covariance matrix and matrix of
  fourth moments. 
}
\usage{
mvnorm.kur.test(X, method = "integration", n.simu = 1000, 
                na.action = na.fail)
}
\arguments{
  \item{X}{a numeric data frame or matrix. }
  \item{method}{defines the method used for the computation of the p-value. The possibilites are
        "integration" (default), "satterthwaite" or "simulation". Details below.}
  \item{n.simu}{if '\code{method}=simulation' this specifies the number of replications in the simulation. }
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
  }
\details{
  This test implements the multivariate normality test based on kurtosis measured by two different scatter estimates
  as described in Kankainen, Taskinen and Oja. The choice here is based on the regular covariance matrix and matrix of
  fourth moments (\code{\link{cov4}}). 
  The limiting distribution of the test statistic W is a linear combination of independent chi-square variables with different degrees of freedom.
  Exact limiting p-values or approximated p-values are obtained by using the function \code{\link[survey]{pchisqsum}}. However Kankainen et al.
  mention that even for n = 200 the convergence can be poor, therefore also p-values simulated under the NULL can be obtained.  
  
}
\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the test statistic W.}
  \item{parameter}{the degrees of freedom for the test statistic W with their weights or the number of replications depending on the chosen method.}
  \item{p.value}{the p-value for the test.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data.}
}
\references{\cite{Kankainen, A., Taskinen, S. and Oja, H. (2007),Tests of multinormality based on location vectors and scatter matrices, \emph{Statistical Methods and Applications}, \bold{16},  357--379.}}
\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}
\seealso{ \code{\link{mvnorm.skew.test}} }
\examples{
X<-rmvnorm(100, c(2, 4, 5))
mvnorm.kur.test(X)
mvnorm.kur.test(X, method = "satt")
mvnorm.kur.test(X, method = "simu")
}
\keyword{htest}
\keyword{ multivariate }
