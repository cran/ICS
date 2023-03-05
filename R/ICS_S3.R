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
                S2_args =  list(), QR = NULL,  whiten = NULL,
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
  #       implement and document.

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
  if (is.function(S2)) {
    if (QR) {
      # whitening is not applicable when we have the QR algorithm
      if (have_whiten && whiten) {
        warning("whitening is not applicable for QR algorithm; ",
              "proceeding without whitening")
      }
      whiten <- FALSE
    } else if (!have_whiten) {
      # use whitening by default
      whiten <- TRUE
    }
  } else {
    # whitening not applicable
    if (have_whiten && whiten) {
      warning("whitening requires 'S2' to be a function; ",
              "proceeding without whitening")
    }
    whiten <- FALSE
  }

  # check argument for centering
  center <- isTRUE(center)
  if ((QR || whiten) && center) {
    warning("centering the data is not applicable when QR algorithm ",
            "or whitening are used; proceeding without centering")
    center <- FALSE
  }

  # check remaining arguments
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
    cf <- S2_args$cf
    if (is.null(cf)) cf <- formals(covW)$cf
    alpha <- S2_args$alpha
    if (is.null(alpha)) alpha <- formals(covW)$alpha

    # compute the second scatter matrix: this works for one-step M-estimators
    S2_Y <- cf*(n-1)/n * crossprod(sweep(Q, 1L, d^alpha, "*"), Q)

    # convert second scatter matrix to class "ICS_scatter" and extract label
    S2_Y <- to_ICS_scatter(S2_Y, label = get_covW_label())
    S2_label <- S2_Y$label

    # reorder columns of the centered data matrix
    X_centered <-  X_centered[, qr_X$pivot]

  } else if (is.function(S2)) {

    ## compute second scatter matrix: if requested, first whiten the data
    ## with respect to the first scatter matrix

    # compute B1 = S1_X^-1/2
    S1_X.eigen <- eigen(S1_X$scatter, symmetric=TRUE)
    B1 <- S1_X.eigen$vectors %*% tcrossprod(diag(S1_X.eigen$values^(-0.5)), S1_X.eigen$vectors)

    # obtain second scatter matrix and whiten the data if requested
    # convert scatter matrices to class "ICS_scatter"
    if (whiten) {
      Y <- X %*% B1
      S2_Y <- get_scatter(Y, fun = S2, args = S2_args, label = S2_label)
    } else {
      # or compute S1_X^-1/2*S2_Y*S1_X^-1/2
      S2.X <- get_scatter(X, fun = S2, args = S2_args)
      S2_Y <- to_ICS_scatter(B1 %*% S2.X$scatter %*% B1, label = S2_label)
    }

  } else {

    ## second scatter matrix supplied by user

    # perform 1st step of ICS
    # compute B1 = S1_X^-1/2
    S1_X.eigen <- eigen(S1_X$scatter, symmetric=TRUE)
    B1 <- S1_X.eigen$vectors %*% tcrossprod(diag(S1_X.eigen$values^(-0.5)), S1_X.eigen$vectors)

    # obtain second scatter matrix
    # or compute S1_X^-1/2*S2_Y*S1_X^-1/2
    # convert scatter matrices to class "ICS_scatter"
    S2.X <- to_ICS_scatter(S2, label = S2_label)
    S2_Y <- B1 %*% S2.X$scatter %*% B1
    S2_Y <- to_ICS_scatter(S2_Y, label = S2_label)

  }

  # compute eigendecomposition of second scatter matrix
  S2_Y_eigen <- eigen(S2_Y$scatter, symmetric = TRUE)
  lambda <- S2_Y_eigen$values
  if (QR) W <- t(qr.solve(R, S2_Y_eigen$vectors))
  else W <- crossprod(S2_Y_eigen$vectors, B1)

  # centering only makes sense if the scatter objects have location components
  if (center && (is.null(S1_X$location) || is.null(S2.X$location))) {
    warning("location component in 'S1' and 'S2' required for centering ",
            "the data; proceeding without centering")
    center <- FALSE
  }

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
    T1_Z <- S1_X$location %*% W
    T2_Z <- S2.X$location %*% W
    # compute generalized skewness values and signs of each component
    gamma <- T1_Z - T2_Z
    skewness_signs <- ifelse(gamma >= 0, 1, -1)
    # fix signs in 'W' matrix so that generalized skewness values are positive
    gamma <- as.vector(skewness_signs * gamma)
    W_final <- sweep(W, 1L, skewness_signs, "*")
    # center the columns of data matrix by the first location estimate
    X_centered <- sweep(X, 2, S1_X$location, "-")

  } else {

    # if requested, scale generalized kurtosis values
    if (scale_lambda) lambda <- lambda / prod(lambda)^(1/p)

    # fix the signs in the 'W' matrix of coefficients
    if (fix_signs == "scores") {
      # compute scores for initial 'W' matrix of coefficients
      Z <- if (QR) tcrossprod(X_centered, W) else tcrossprod(X, W)
      # compute skewness values and signs of each component
      gamma <- colMeans(Z) - apply(Z, 2L, median)
      skewness_signs <- ifelse(gamma > 0, 1, -1)
      # fix signs in 'W' matrix so that all components are right-skewed
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

  }

  # compute the scores
  if (QR || center) Z_final <- tcrossprod(X_centered, W_final)
  else Z_final <- tcrossprod(X, W_final)
  colnames(Z_final) <- paste("IC", 1:p, sep = ".")

  # construct object to be returned
  res <- list(lambda = lambda,
              W = W_final,
              scores = Z_final,
              gamma = gamma,
              S1_label = S1_label,
              S2_label = S2_label,
              S1_args = S1_args,
              S2_args = S2_args,
              QR = QR,
              whiten = whiten,
              center = center,
              scale_lambda = scale_lambda,
              fix_signs = fix_signs)
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


# compute the matrix square root or the inverse thereof
# TODO: there is already an internal function mat.sqrt() in the package, but
#       that one doesn't allow to compute the inverse
mat_sqrt <- function(A, inv = FALSE) {
  # initializations
  power <- if (inv) -0.5 else 0.5
  # compute eigendecomposition
  eigen_A <- eigen(A)
  # compute matrix square root or inverse
  eigen_A$vectors %*% tcrossprod(diag(eigen_A$values^power), eigen_A$vectors)
}

# in QR algorithm, the function for the second scatter matrix is not actually
# applied, so we need a different way to get the label for covW()
get_covW_label <- function() "COVW"
