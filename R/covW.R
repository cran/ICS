#'  One-step M-estimator covW
#'
#' @param X numeric nxp data matrix or dataframe.
#' @param na.action a function which indicates what should happen when the data contain 'NA's. Default is to fail.
#' @param alpha parameter of the one-step M-estimator. By default equals to 1.
#' @param cf consistency factor of the one-step M-estimator. By default equals to 1
#'
#' @return
#' @export
#'
#' @examples
covW <- function (X, na.action = na.fail, alpha = 1, cf = 1 )
{
  X <- na.action(X)
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (p < 2)
    stop("'X' must be at least bivariate")
  Xmeans <- colMeans(X)
  di <- mahalanobis(X, Xmeans, cov(X))
  X.centered <- sweep(X, 2, Xmeans)
  v.tilde <- 1/n*cf * t( X.centered) %*% diag(di^alpha) %*%  X.centered

  return(v.tilde)
}
