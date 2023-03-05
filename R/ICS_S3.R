# Scatter functions returning class "ICS_scatter" -----

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
ICS_cov <- function(x, location = c("mean", "none")) {
  # initializations
  location <- match.arg(location)
  # compute center and scatter estimates
  location <- if (location == "mean") colMeans(x)
  out <- list(location = location, scatter = cov(x), label = "COV")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
## TODO: Do we need to allow passing other arguments to cov4? Then we also need
##       to figure out what to do with the location argument of cov4().
ICS_cov4 <- function(x, location = c("mean3", "mean", "none")) {
  # initializations
  location <- match.arg(location)
  # compute center and scatter estimates
  location <- switch(location, "mean3" = mean3(x), "mean" = colMeans(x))
  out <- list(location = location, scatter = cov4(x), label = "COV4")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}


#' Title
#'
#' @param x
#' @param na.action
#' @param alpha
#' @param cf
#'
#' @return
#' @export
#'
#' @examples
ICS_covW <- function(x, na.action = na.fail, alpha = 1, cf = 1) {
  # compute center and scatter estimates
  scatter <- covW(x, na.action = na.fail, alpha = alpha, cf = cf)
  out <- list(location = colMeans(x), scatter = scatter,
              label = get_covW_label())
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}


# Main function to compute ICS -----

#' Title
#'
#' @param X
#' @param S1
#' @param S2
#' @param S1_args
#' @param S2_args
#' @param QR
#' @param whiten
#' @param center logical indicating whether to center the ICS coordinates (scores)
#' @param fix_signs character string specifying how to fix_signs the ICS
#                 coordinates, either 'scores' or 'eigenvectors'
#' @param scale_lambda
#' @param na.action
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ICS <- function(X, S1 = ICS_cov, S2 = ICS_cov4, S1_args = list(),
                S2_args =  list(), QR = NULL, whiten = NULL,
                center = FALSE, scale_lambda = FALSE,
                fix_signs = c("scores", "W"),
                na.action = na.fail, ...) {

  # make sure we have a suitable data matrix
  X <- na.action(X)
  X <- as.matrix(X)
  p <- ncol(X)
  if (p < 2L) stop("'X' must be at least bivariate")

  # obtain default labels for scatter matrices
  S1_label <- deparse(substitute(S1))
  S2_label <- deparse(substitute(S2))

  # TODO: we first check 'QR' and 'whiten', and 'center' is only applicable if
  #       both are FALSE.  That is, if a user supplies a function for S2 and
  #       wants to center, simply setting 'center = TRUE' is not enough, they
  #       also have to set 'whiten = FALSE' and/or 'QR = FALSE'. Is this the
  #       behavior we want? We could check which arguments are supplied, and
  #       set defaults for others, although that becomes rather cumbersome to
  #       implement and document. My vote is to leave it as is, in particular
  #       because I think that we can use centering when whitening.

  # check argument for QR algorithm
  have_QR <- !is.null(QR)
  if (have_QR) QR <- isTRUE(QR)
  # QR algorithm requires a certain class of scatter pairs supplied as functions
  have_cov <- identical(S1, cov) || identical(S1, ICS_cov)
  have_covW <- identical(S2, covW) || identical(S2, ICS_covW)
  if (have_cov && have_covW) {
    # use QR algorithm by default
    if (!have_QR) QR <- TRUE
  } else {
    # QR algorithm not applicable for other scatter pairs
    if (have_QR && QR) {
      warning("QR algorithm is not applicable; proceeding without it")
    }
    QR <- FALSE
  }

  # check argument for whitening
  have_whiten <- !is.null(whiten)
  if (have_whiten) whiten <- isTRUE(whiten)
  # whitening requires S2 to be a function
  if (QR) {
    # whitening is not applicable when we have the QR algorithm
    if (have_whiten && whiten) {
      warning("whitening is not applicable for QR algorithm; ",
              "proceeding without whitening")
    }
    whiten <- FALSE
  } else if (is.function(S2)) {
    # use whitening by default
    if (!have_whiten) whiten <- TRUE
  } else {
    # whitening not applicable
    if (have_whiten && whiten) {
      warning("whitening requires 'S2' to be a function; ",
              "proceeding without whitening")
    }
    whiten <- FALSE
  }

  # check argument for centering
  # TODO: see below, I think we can use centering when whitening
  # if ((QR || whiten) && center) {
  #   warning("centering the data is not applicable when QR algorithm ",
  #           "or whitening are used; proceeding without centering")
  #   center <- FALSE
  # }

  # check remaining arguments
  center <- isTRUE(center)
  scale_lambda <- isTRUE(scale_lambda)
  fix_signs <- match.arg(fix_signs)

  # obtain first scatter matrix
  if (is.function(S1)) {
    S1_X <- get_scatter(X, fun = S1, args = S1_args, label = S1_label)
  } else S1_X <- to_ICS_scatter(S1, label = S1_label)
  # update label for first scatter matrix
  S1_label <- S1_X$label

  # obtain second scatter matrix
  if (QR) {

    ## perform numerically stable algorithm based on QR decomposition (second
    ## scatter matrix is specified as a function for a one-step M-estimator)

    # further initializations
    n <- nrow(X)

    # center the columns of data matrix by the mean
    X_centered = sweep(X, 2L, colMeans(X), "-")

    # reorder rows by decreasing by infinity norm (maximum in absolute value)
    norm_inf <- apply(abs(X_centered), 1L, max)
    order_rows <- order(norm_inf, decreasing = TRUE)
    X_reordered <- X_centered[order_rows, ]

    # compute QR decomposition with column pivoting from LAPACK
    qr_X <- qr(X_reordered / sqrt(n-1), LAPACK = TRUE)

    # extract components of the QR decomposition
    Q <- qr.Q(qr_X)  # should be nxp
    R <- qr.R(qr_X)  # should be pxp

    # compute squared Mahalanobis distances (leverage scores)
    d <- (n-1) * rowSums(Q^2)

    # extract relevant arguments for one-step M-estimator
    alpha <- S2_args$alpha
    if (is.null(alpha)) alpha <- formals(covW)$alpha
    cf <- S2_args$cf
    if (is.null(cf)) cf <- formals(covW)$cf

    # compute the second scatter matrix: this works for one-step M-estimators
    S2_Y <- cf*(n-1)/n * crossprod(sweep(Q, 1L, d^alpha, "*"), Q)

    # convert second scatter matrix to class "ICS_scatter" and update label
    S2_Y <- to_ICS_scatter(S2_Y, label = get_covW_label())
    S2_label <- S2_Y$label

    # reorder columns of the centered data matrix
    X_centered <-  X_centered[, qr_X$pivot]

    # we don't have a location estimate for S2, so centering is not applicable
    # TODO: Actually, we do center the scores, but we don't fix the signs based
    #       on the generalized skewness values. Hence it is rather confusing
    #       that 'center' is set to FALSE. Perhaps we need to change argument
    #       'fix_signs' to take values "location" (center = TRUE), "skewness"
    #       (currently "scores"), and "W" (same as currently)? Then argument
    #       'center' would only control whether the scores should be computed
    #       from the centered data matrix (so that it should be set to TRUE
    #       when the QR algorithm is used).
    if (center) {
      warning("centering the data is not applicable for QR algorithm; ",
              "proceeding without centering")
      center <- FALSE
    }

  } else {

    # compute inverse of the square root of the first scatter matrix
    W1 <- mat_sqrt(S1_X$scatter, inverse = TRUE)

    # obtain second scatter matrix
    if (whiten) {
      # whiten the data with respect to the first scatter matrix
      Y <- X %*% W1
      # compute second scatter matrix on whitened data and update label
      S2_Y <- get_scatter(Y, fun = S2, args = S2_args, label = S2_label)
      S2_label <- S2_Y$label
      # check if we have location component
      is_null_T2 <- is.null(S2_Y$location)
    } else {
      # obtain second scatter matrix on original data matrix
      if (is.function(S2)) {
        S2_X <- get_scatter(X, fun = S2, args = S2_args, label = S2_label)
      } else S2_X <- to_ICS_scatter(S2, label = S2_label)
      # update label for second scatter matrix
      S2_label <- S2_X$label
      # transform second scatter matrix
      S2_Y <- to_ICS_scatter(W1 %*% S2_X$scatter %*% W1, label = S2_label)
      # check if we have location component
      is_null_T2 <- is.null(S2_X$location)
    }

    # if centering is requested, check if we have location components
    if (center && (is.null(S1_X$location) || is_null_T2)) {
      warning("location component in 'S1' and 'S2' required for centering ",
              "the data; proceeding without centering")
      center <- FALSE
    }

  }

  # compute eigendecomposition of second scatter matrix
  S2_Y_eigen <- eigen(S2_Y$scatter, symmetric = TRUE)
  lambda <- S2_Y_eigen$values
  if (QR) W <- t(qr.solve(R, S2_Y_eigen$vectors))
  else W <- crossprod(S2_Y_eigen$vectors, W1)

  # scale generalized kurtosis values and fix signs in the 'W' matrix
  if (center) {

    # overrride values of related arguments
    # TODO: Would it make sense to allow scaling the generalized kurtosis
    #       values, or does the scale of the generalized kurtosis values
    #       matter for the interpretation of the distances? I don't see
    #       how it would matter since below the scores are unchanged when
    #       we scale the generalized kurtosis values.
    if (scale_lambda) {
      warning("when centering the data, scaling the generalized kurtosis ",
              "values is not applicable; setting 'scale_lambda' to FALSE")
      scale_lambda <- FALSE
    }
    if (fix_signs != "scores") {
      warning("when centering the data, signs of 'W' matrix are set based on ",
              "generalized skewness of components; setting 'fix_signs' to ",
              '"scores"')
      fix_signs <- "scores"
    }

    # compute location estimates of the scores
    T1_X <- S1_X$location
    T1_Z <- T1_X %*% W
    # T2_Z <- S2_X$location %*% W
    # TODO: Aurore, in an earlier version you mentioned that it's necessary
    #       in case of whitening to apply S2 again for computing the location
    #       estimate on the original data matrix X. But that is not true if
    #       we have an affine equivariant scatter matrix, then we can simply
    #       backtransform the location estimate for Y.
    T2_X <- if (whiten) S2_Y$location %*% solve(W1) else S2_X$location
    T2_Z <- T2_X %*% W
    # center the columns of data matrix by the first location estimate
    X_centered <- sweep(X, 2, T1_X, "-")

  }

  # if requested, scale generalized kurtosis values
  if (scale_lambda) lambda <- lambda / prod(lambda)^(1/p)

  # fix the signs in the 'W' matrix of coefficients
  if (fix_signs == "scores") {

    if (center) {
      # compute generalized skewness values of each component
      gamma <- T1_Z - T2_Z
    } else {
      # compute scores for initial 'W' matrix of coefficients
      Z <- if (QR) tcrossprod(X_centered, W) else tcrossprod(X, W)
      # compute skewness values of each component
      gamma <- colMeans(Z) - apply(Z, 2L, median)
    }

    # compute signs of (generalized) skewness values for each component
    # TODO: For centering, we used '>=' below, but otherwise we used '>'.
    #       I put >= to have the same in both cases, so that the sign is only
    #       changed if the (generalized) skewness is negative. In practice,
    #       this shouldn't really matter, as the skewness value will never be
    #       exactly 0.
    skewness_signs <- ifelse(gamma >= 0, 1, -1)
    # fix signs in 'W' matrix so that generalized skewness values are positive
    gamma <- as.vector(skewness_signs * gamma)
    W_final <- sweep(W, 1L, skewness_signs, "*")

  } else {
    # fix signs in 'W' matrix so that the maximum element per row is positive
    # and that each row has norm 1
    row_signs <- apply(W, 1L, .sign.max)
    row_norms <- sqrt(rowSums(W^2))
    W_final <- sweep(W, 1L, row_norms * row_signs, "/")
    # we don't have (generalized) skewness values in this case
    gamma <- NULL
  }

  # compute the scores
  # TODO: perhaps we want to enable that we return the non-centered scores
  #       when we use the QR algorithm (see also above)?
  if (QR || center) Z_final <- tcrossprod(X_centered, W_final)
  else Z_final <- tcrossprod(X, W_final)
  colnames(Z_final) <- paste("IC", 1:p, sep = ".")

  # construct object to be returned
  res <- list(lambda = lambda, W = W_final, scores = Z_final, gamma = gamma,
              S1_label = S1_label, S2_label = S2_label, S1_args = S1_args,
              S2_args = S2_args, QR = QR, whiten = whiten, center = center,
              scale_lambda = scale_lambda, fix_signs = fix_signs)
  class(res) <- "ICS"
  res

}


# Methods for class "ICS" -----

#' @method summary ICS
#' @export
summary.ICS <- function(object) {
  return(object)
}

#' @method coef ICS
#' @export
coef.ICS <- function(object) {
  return(object$W)
}

#' @method fitted ICS
#' @export
fitted.ICS <- function(object, index=NULL) {
  p <- ncol(object$scores)
  if (is.null(index) == FALSE && max(index)>p) stop("undefined columns selected")
  Mix <- solve(object$W)
  if (is.null(index)) index=1:p
  fits <- tcrossprod(as.matrix(object$scores[,index]), Mix[,index])
  return(as.data.frame(fits))
}


#' @method plot ICS
#' @export
#' @importFrom graphics pairs
plot.ICS <- function(x,index=NULL,...){
  p<-ncol(x$W)
  if (is.null(index) & p<=6) pairs(x$scores,...)
  if (is.null(index) & p>6) pairs(x$scores[,c(1:3,p-2:0)],...)
  if (length(index)==1) stop("index must be NULL or at least a vector of length 2")
  if (length(index)>1) pairs(x$scores[,index],...)
}

#' @method print ICS
#' @export
print.ICS <- function(object, digits = 4){
  cat("\nICS with the following parameters: \n")
  cat("S1:", object$S1_label)
  if (length(object$S1_args) > 0){
    sapply(1:length(object$S1_args), function(i)
      cat("\n", paste0(names(object$S1_args)[i], ":"), object$S1_args[[i]]))
  }
  cat("\nS2:", object$S2_label)
  if (length(object$S2_args)>0){
    sapply(1:length(object$S2_args), function(i)
      cat("\n", paste0(names(object$S2_args)[i], ":"), object$S2_args[[i]]))
  }
  cat("\nQR:", object$QR)
  cat("\nwhiten:", object$whiten)
  cat("\nscale_lambda:", object$scale_lambda)
  cat("\nfix_signs:", object$fix_signs)
  cat("\ncenter:", object$center)
  cat("\n")
  cat("\nThe generalized kurtosis measures (lambda) of the components are:\n")
  print(format(round(object$lambda, digits)), quote = FALSE)
  cat("\n")
  cat("\nThe W matrix of coefficients is:\n")
  print(round(object$W, digits))
  invisible(object)
}


# Internal functions -----


## apply a scatter function to the data matrix
# X ......... data matrix
# fun ....... function to compute a scatter matrix
# args ...... list of additional arguments to be passed to the function
# convert ... logical indicating whether the scatter matrix should be converted
#             to class "ICS_scatter"
# label ..... typically constructed beforehand via deparse(substitute())
# TODO: perhaps argument 'convert' is unnecessary?
get_scatter <- function(X, fun = cov, args = list(), convert = TRUE, label) {
  if (length(args) == 0) scatter <- fun(X)
  else {
    # there should be a more efficient way of doing this, for example via call()
    args <- c(list(X), args)
    scatter <- do.call(fun, args)
  }
  # if requested, convert to class "ICS_scatter": use the function name as
  # label if the function does not return an object of this class already
  if (convert) scatter <- to_ICS_scatter(scatter, label = label)
  # return scatter matrix
  scatter
}


## convert scatter matrices to the required class

to_ICS_scatter <- function(object, ...) UseMethod("to_ICS_scatter")

# object ... a scatter matrix
# label .... typically constructed beforehand via deparse(substitute())
to_ICS_scatter.matrix <- function(object, label, ...) {
  # convert to class "ICS_scatter" with empty location estimate
  out <- list(location = NULL, scatter = object, label = label)
  class(out) <- "ICS_scatter"
  out
}

# object ... typically a list with components 'location' and 'scatter', but it
#            can also have a component 'label'
# label .... typically constructed beforehand via deparse(substitute()) and
#            ignored if the list already has a component 'label'
to_ICS_scatter.list <- function(object, label) {
  # check that list has a component 'scatter'
  scatter <- object$scatter
  if (is.null(scatter)) {
    stop("list should have components 'scatter' (required) as well as ",
         "'location' and 'label' (optional)")
  }
  # check if there already is a component 'label' in the list
  if (!is.null(object$label)) label <- object$label
  # convert to class "ICS_scatter"
  out <- list(location = object$location, scatter = scatter, label = label)
  class(out) <- "ICS_scatter"
  out
}

to_ICS_scatter.ICS_scatter <- function(object, ...) object


# compute the matrix square root (or the inverse thereof) of a symmetric matrix
# TODO: there is already an internal function mat.sqrt() in the package, but it
#       doesn't assume symmetry and doesn't allow to compute the inverse
mat_sqrt <- function(A, inverse = FALSE) {
  # initializations
  power <- if (inverse) -0.5 else 0.5
  # compute eigendecomposition
  eigen_A <- eigen(A, symmetric = TRUE)
  # compute matrix square root or inverse
  eigen_A$vectors %*% tcrossprod(diag(eigen_A$values^power), eigen_A$vectors)
}

# in QR algorithm, the function for the second scatter matrix is not actually
# applied, so we need a different way to get the label for covW()
get_covW_label <- function() "COVW"
