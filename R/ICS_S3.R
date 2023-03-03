

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
ICS_cov <- function(x) {
  # compute center and scatter estimates
  out <- list(location = colMeans(x), scatter = cov(x), label = "COV")
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
ICS_cov4 <- function(x) {
  # compute center and scatter estimates
  out <- list(location = colMeans(x), scatter = cov4(x), label = "COV4")
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
ICS_Mean3Cov4 <- function(x) {
  # compute center and scatter estimates
  out <- list(location = mean3(x), scatter = cov4(x), label = "COV4")
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
ICS_covW <- function(x, na.action = na.fail, alpha = 1, cf = 1 ) {
  # compute center and scatter estimates
  out <- list(location = colMeans(x), scatter = covW(x, na.action = na.fail, alpha = alpha, cf = cf), label = "COVW")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}


# S2_Y <- function(X, S1, S2, ...) UseMethod("S2_Y", S2)
#
# S2_Y.function <- function(X, S1.X, S2, S2_args = list(), QR = FALSE,
#                           whiten = !QR, S1_label, S2_label, ...) {
#   # Initialization
#   S2.X <- S2
#
#   if (QR) {
#
#     # Details
#     n <- nrow(X)
#     p <- ncol(X)
#
#     # perform algorithm based on QR decomposition
#     # centering by the mean
#     X.c = sweep(X, 2, colMeans(X), "-")
#
#     # Permutation by rows
#     # decreasing by infinity norm:  absolute maximum
#     norm_inf <- apply(X.c, 1, function(x) max(abs(x)))
#     order_rows <- order(norm_inf, decreasing = TRUE)
#     X_row_per <- X.c[order_rows,]
#
#     # QR decomposition of X with column pivoting from LAPACK
#     qr_X <- qr( 1/sqrt(n-1)*X_row_per, LAPACK = TRUE)
#
#     # R should be pxp
#     R <- qr.R(qr_X)
#
#     # Q should be nxp
#     Q <- qr.Q(qr_X)
#
#     # computation of d_i
#     d_i <- (n-1)*apply(Q, 1, function(x) sum(x^2))
#
#     # Spectral decomposition (SEP) of S2.Y
#     # Computation of S2 - with general one-step M-estimators
#     S2.Y <- S2_args$cf*1/n*(n-1)*t(Q) %*% diag(d_i^S2_args$alpha) %*% Q
#     # convert scatter matrices to class "ICS_scatter"
#     S2.Y <- to_ICS_scatter(S2.Y, label = S2_label)
#
#     # Reorder the rows and cols of X
#     X_Z <-  X_row_per[order(order_rows),qr_X$pivot]
#     B1 <- NULL
#
#   } else {
#
#     # compute B1 = S1.X^-1/2
#     S1.X.eigen <- eigen(S1.X$scatter, symmetric=TRUE)
#     B1 <- S1.X.eigen$vectors %*% tcrossprod(diag(S1.X.eigen$values^(-0.5)), S1.X.eigen$vectors)
#
#     # obtain second scatter matrix and whiten the data if requested
#     # convert scatter matrices to class "ICS_scatter"
#     if (whiten) {
#       Y <- X %*% B1
#       S2.Y <- get_scatter(Y, fun = S2, args = S2_args, label = S2_label)
#     } else {
#       # or compute S1.X^-1/2*S2.Y*S1.X^-1/2
#       S2.X <- get_scatter(X, fun = S2, args = S2_args)
#       S2.Y <- to_ICS_scatter(B1 %*% S2.X$scatter %*% B1, label = S2_label)
#     }
#     X_Z <- X
#     R <- NULL
#
#   }
#
#   # return list of results
#   list(X, S2.X = S2.X, S2.Y = S2.Y, B1 = B1, R = R, X_Z = X_Z, ...)
# }
#
#
# S2_Y.matrix <- function(X, S1.X, S2, S1_label, S2_label, ...) {
#
#   # perform 1st step of ICS
#   # compute B1 = S1.X^-1/2
#   S1.X.eigen <- eigen(S1.X$scatter, symmetric=TRUE)
#   B1 <- S1.X.eigen$vectors %*% tcrossprod(diag(S1.X.eigen$values^(-0.5)), S1.X.eigen$vectors)
#
#   # obtain second scatter matrix
#   # or compute S1.X^-1/2*S2.Y*S1.X^-1/2
#   # convert scatter matrices to class "ICS_scatter"
#   S2.X <- to_ICS_scatter(S2, label = S2_label)
#   S2.Y <- B1 %*% S2.X$scatter %*% B1
#   S2.Y <- to_ICS_scatter(S2.Y, label = S2_label)
#
#   # return list of results
#   list(S2.X = S2.X, S2.Y = S2.Y, B1 = B1, X_Z = X, ...)
# }
#
# S2_Y.ICS_scatter <- S2_Y.matrix




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
      warning("whitening not applicable for QR algorithm; ",
              "proceeding without whitening")
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

  # TODO: checks for other input arguments
  # E.g., 'center' should only be TRUE if 'QR' and 'whiten' are both FALSE

  # obtain first scatter matrix or convert to class "ICS_scatter"
  if (is.function(S1)) {
    S1.X <- get_scatter(X, fun = S1, args = S1_args, label = S1_label)
  } else S1.X <- to_ICS_scatter(S1, label = S1_label)
  # update label for first scatter matrix
  S1_label <- S1.X$label

  # # compute S2_Y
  # S2_Y <- S2_Y(X, S1.X = S1.X, S2 = S2, S2_args = S2_args,
  #              QR = QR, whiten = whiten, S1_label, S2_label, ... )
  # X_Z <- S2_Y$X_Z

  # TODO: move computation of S2_Y back to ICS() instead of separate function
  if (is.function(S2)) {

    # Initialization
    # TODO: check if this can be avoided
    S2.X <- S2

    if (QR) {

      # Details
      n <- nrow(X)
      p <- ncol(X)

      # perform algorithm based on QR decomposition
      # centering by the mean
      X.c = sweep(X, 2, colMeans(X), "-")

      # Permutation by rows
      # decreasing by infinity norm:  absolute maximum
      norm_inf <- apply(X.c, 1, function(x) max(abs(x)))
      order_rows <- order(norm_inf, decreasing = TRUE)
      X_row_per <- X.c[order_rows,]

      # QR decomposition of X with column pivoting from LAPACK
      qr_X <- qr( 1/sqrt(n-1)*X_row_per, LAPACK = TRUE)

      # R should be pxp
      R <- qr.R(qr_X)

      # Q should be nxp
      Q <- qr.Q(qr_X)

      # computation of d_i
      d_i <- (n-1)*apply(Q, 1, function(x) sum(x^2))

      # Spectral decomposition (SEP) of S2.Y
      # Computation of S2 - with general one-step M-estimators
      S2.Y <- S2_args$cf*1/n*(n-1)*t(Q) %*% diag(d_i^S2_args$alpha) %*% Q
      # convert scatter matrices to class "ICS_scatter"
      S2.Y <- to_ICS_scatter(S2.Y, label = S2_label)

      # Reorder the rows and cols of X
      X_Z <-  X_row_per[order(order_rows),qr_X$pivot]

      # TODO: I don't think that this is necessary
      # B1 <- NULL

    } else {

      # compute B1 = S1.X^-1/2
      S1.X.eigen <- eigen(S1.X$scatter, symmetric=TRUE)
      B1 <- S1.X.eigen$vectors %*% tcrossprod(diag(S1.X.eigen$values^(-0.5)), S1.X.eigen$vectors)

      # obtain second scatter matrix and whiten the data if requested
      # convert scatter matrices to class "ICS_scatter"
      if (whiten) {
        Y <- X %*% B1
        S2.Y <- get_scatter(Y, fun = S2, args = S2_args, label = S2_label)
      } else {
        # or compute S1.X^-1/2*S2.Y*S1.X^-1/2
        S2.X <- get_scatter(X, fun = S2, args = S2_args)
        S2.Y <- to_ICS_scatter(B1 %*% S2.X$scatter %*% B1, label = S2_label)
      }

      # copy data matrix
      # TODO: check if this can be avoided
      X_Z <- X

      # TODO: I don't think that this is necessary
      # R <- NULL

    }

  } else {

    # perform 1st step of ICS
    # compute B1 = S1.X^-1/2
    S1.X.eigen <- eigen(S1.X$scatter, symmetric=TRUE)
    B1 <- S1.X.eigen$vectors %*% tcrossprod(diag(S1.X.eigen$values^(-0.5)), S1.X.eigen$vectors)

    # obtain second scatter matrix
    # or compute S1.X^-1/2*S2.Y*S1.X^-1/2
    # convert scatter matrices to class "ICS_scatter"
    S2.X <- to_ICS_scatter(S2, label = S2_label)
    S2.Y <- B1 %*% S2.X$scatter %*% B1
    S2.Y <- to_ICS_scatter(S2.Y, label = S2_label)

    # copy data matrix
    # TODO: check if this can be avoided
    X_Z <- X

  }

  # decomposition of S2
  S2.Y.eigen <- eigen(S2.Y$scatter, symmetric=TRUE)
  U2 <- S2.Y.eigen$vectors
  lambda <- S2.Y.eigen$values
  if (QR) {
    # Eigenvectors by rows
    Rinv <- qr.solve(R)
    B <- t(U2) %*% t(Rinv)

  } else {
    B <- crossprod(U2, B1)
  }

  # choosing the signs of B
  # centering only makes sense if the scatter objects have location components
  center <- isTRUE(center)
  if (center){
    if (is.function(S2)) {
      # TODO: The first part is unnecessary with checks for argument 'center'
      #       (see also TODO comment below). If S2 is a function, we have S2.X
      #       in most cases so that we can check if we have a location
      #       component. As far as I can see, we only don't have an
      #       "ICS_scatter" object S2_X if the QR algorithm or whitening have
      #       been used.
      if (is.null(S1.X$location) && is.function(S2)) {
        center <- FALSE
        warning("'S1' needs to have a location component for centering ",
                "the data; proceeding without centering")
      }
    } else {
      if (is.null(S1.X$location) || is.null(S2.X$location)) {
        center <- FALSE
        warning("'S1' and 'S2' need to have a location component for ",
                "centering the data; proceeding without centering")
      }
    }
  }


  # TODO: add more comments for this next part
  if (center) {

    # TODO: with proper checks for argument 'center', this part can be removed
    #       (i.e., 'center' is only applicable when we don't have QR algorithm
    #       or whitening, which can be checked before computing S2_Y)
    if (is.function(S2)) {
      S2.X <- get_scatter(X, fun = S2, args = S2_args)
    } else {
      S2.X <- to_ICS_scatter(S2, label = deparse(substitute(S2)))
    }

    T1.Z <- S1.X$location %*% B
    T2.Z <- S2.X$location %*% B

    gamma <- T1.Z - T2.Z
    skew.signs <- ifelse(gamma >= 0, 1, -1)
    gamma <- as.vector(skew.signs*gamma)

    B.res <- sweep(B, 1, skew.signs, "*")
    X_Z <- sweep(X_Z, 2, S1.X$location, "-")
    # What should it be

    fix_signs <- "scores"
    scale_lambda <- FALSE

  } else {

    # further initializations
    fix_signs <- match.arg(fix_signs)
    scale_lambda <- isTRUE(scale_lambda)

    if (scale_lambda) lambda <- lambda/prod(lambda)^(1/p)
    if (fix_signs == "W") {
      row.signs <- apply(B, 1, .sign.max)
      row.norms <- sqrt(rowSums((B)^2))
      B.res <- sweep(B, 1, row.norms * row.signs, "/")
      gamma <- NULL
    }
    if (fix_signs == "scores" && !center) {
      Z1 <- tcrossprod(X_Z, B)
      gamma <- colMeans(Z1) - apply(Z1, 2, median)
      skew.signs <- ifelse(gamma > 0, 1, -1)
      gamma <- as.vector(skew.signs*gamma)
      B.res <- sweep(B, 1, skew.signs, "*")
    }

  }

  # compute the scores
  Z <- tcrossprod(X_Z, B.res)
  colnames(Z) <- paste(rep("IC", p), 1:p, sep = ".")

  # construct object to be returned
  res <- list(lambda = lambda,
              W = B.res,
              scores = Z,
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



## internal function to convert scatter matrices to the required class
## internal function to apply a scatter function to the data matrix
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
