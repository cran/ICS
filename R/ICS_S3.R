# Scatter functions returning class "ICS_scatter" -----

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
ICS_cov <- function(x, location = TRUE) {
  # initializations
  x <- as.matrix(x)
  location <- isTRUE(location)
  # compute location and scatter estimates
  location <- if (location) colMeans(x)
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
## TODO: Should the default location estimate be "mean3" or "mean"? I think
##       that mean3 was only needed for fixing the signs in ics2? But it's
##       more natural to use the mean as the default location estimate.
## TODO: Do we need to allow passing other arguments to cov4? Then we also need
##       to figure out what to do with the location argument of cov4().
ICS_cov4 <- function(x, location = c("mean3", "mean", "none")) {
  # initializations
  x <- as.matrix(x)
  if (is.character(location)) location <- match.arg(location)
  else if (isTRUE(location)) location <- "mean3"
  else location <- "none"
  # compute location and scatter estimates
  location <- switch(location, "mean3" = mean3(x), "mean" = colMeans(x))
  out <- list(location = location, scatter = cov4(x), label = get_cov4_label())
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}


#' Title
#'
#' @param x
#' @param alpha
#' @param cf
#'
#' @return
#' @export
#'
#' @examples
ICS_covW <- function(x, location = TRUE, alpha = 1, cf = 1) {
  # initializations
  x <- as.matrix(x)
  location <- isTRUE(location)
  # compute location and scatter estimates
  location <- if (location) colMeans(x)
  scatter <- covW(x, alpha = alpha, cf = cf)
  out <- list(location = location, scatter = scatter, label = get_covW_label())
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
ICS_covAxis <- function(x, location = TRUE) {
  # initializations
  x <- as.matrix(x)
  location <- isTRUE(location)
  # compute location and scatter estimates
  location <- if (location) colMeans(x)
  out <- list(location = location, scatter = covAxis(x),
              label = get_covAxis_label())
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
#                 coordinates, either "scores" or "W"
#' @param scale_gen_kurtosis
#' @param na.action
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ICS <- function(X, S1 = ICS_cov, S2 = ICS_cov4, S1_args = list(),
                S2_args = list(), QR = NULL, whiten = NULL,
                center = FALSE, scale_gen_kurtosis = FALSE,
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
  have_cov4 <- identical(S2, cov4) || identical(S2, ICS_cov4)
  have_covW <- identical(S2, covW) || identical(S2, ICS_covW)
  have_covAxis <- identical(S2, covAxis) || identical(S2, ICS_covAxis)
  if (have_cov && (have_cov4 || have_covW || have_covAxis)) {
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

  # check remaining arguments
  center <- isTRUE(center)
  scale_gen_kurtosis <- isTRUE(scale_gen_kurtosis)
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

    # obtain column means: if ICS_cov() is used, we may already have them in
    # location component of S1_X (S1 is either cov() or ICS_cov())
    T1_X <- S1_X$location
    if (is.null(T1_X)) T1_X <- colMeans(X)
    # center the columns of data matrix by the mean
    X_centered <- sweep(X, 2L, T1_X, "-")

    # reorder rows by decreasing by infinity norm (maximum in absolute value)
    norm_inf <- apply(abs(X_centered), 1L, max)
    order_rows <- order(norm_inf, decreasing = TRUE)
    X_reordered <- X_centered[order_rows, ]

    # compute QR decomposition with column pivoting from LAPACK: note that this
    # changes the order of the columns internally, which we need to take into
    # account when returning the matrix W of coefficients (or other output)
    qr_X <- qr(X_reordered / sqrt(n-1), LAPACK = TRUE)

    # extract components of the QR decomposition
    Q <- qr.Q(qr_X)  # should be nxp
    R <- qr.R(qr_X)  # should be pxp

    # compute squared Mahalanobis distances (leverage scores)
    d <- (n-1) * rowSums(Q^2)

    # obtain arguments for one-step M-estimator and update label for S2
    if (have_cov4) {
      alpha <- 1
      cf <- 1 / (p+2)
      S2_label <- get_cov4_label()
    } else if (have_covAxis) {
      alpha <- -1
      cf <- p
      S2_label <- get_covAxis_label()
    } else {
      # COVW: get arguments from supplied list or function defaults
      alpha <- S2_args$alpha
      if (is.null(alpha)) alpha <- formals(covW)$alpha
      cf <- S2_args$cf
      if (is.null(cf)) cf <- formals(covW)$cf
      S2_label <- get_covW_label()
    }

    # compute the second scatter matrix: this works for one-step M-estimators
    S2_Y <- cf*(n-1)/n * crossprod(sweep(Q, 1L, d^alpha, "*"), Q)

    # convert second scatter matrix to class "ICS_scatter"
    S2_Y <- to_ICS_scatter(S2_Y, label = S2_label)

    # set flag that we don't have location component in S2_X for fixing
    # signs in W (since we don't have S2_X)
    missing_T2_X <- TRUE

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
      # set flag that we don't have location component in S2_X for fixing
      # signs in W (since we don't have S2_X)
      missing_T2_X <- TRUE
    } else {
      # obtain second scatter matrix on original data matrix
      if (is.function(S2)) {
        S2_X <- get_scatter(X, fun = S2, args = S2_args, label = S2_label)
      } else S2_X <- to_ICS_scatter(S2, label = S2_label)
      # update label for second scatter matrix
      S2_label <- S2_X$label
      # transform second scatter matrix
      S2_Y <- to_ICS_scatter(W1 %*% S2_X$scatter %*% W1, label = S2_label)
      # check if we have location component in S2_X for fixing signs in W
      T2_X <- S2_X$location
      missing_T2_X <- is.null(T2_X)
    }

    # if requested, center the columns of the data matrix
    if (center) {
      # center by the location estimate from S1_X or give warning
      # if S1_X doesn't have a location component
      T1_X <- S1_X$location
      if (is.null(T1_X)) {
        warning("location component in 'S1' required for centering the data; ",
                "proceeding without centering")
        center <- FALSE
      } else X_centered <- sweep(X, 2L, T1_X, "-")
    }

  }

  # compute eigendecomposition of second scatter matrix
  S2_Y_eigen <- eigen(S2_Y$scatter, symmetric = TRUE)
  # extract generalized kurtosis values
  gen_kurtosis <- S2_Y_eigen$values
  if (!all(is.finite(gen_kurtosis))) {
    if (scale_gen_kurtosis) {
      suffix <- "; proceeding without scaling them"
      scale_gen_kurtosis <- FALSE
    } else suffix <- ""
    warning("some generalized kurtosis values are non-finite", suffix)
  }
  # obtain coefficient matrix of the linear transformation
  if (QR) {
    # obtain matrix of coefficients in internal order from column pivoting
    W <- t(qr.solve(R, S2_Y_eigen$vectors))
    # reorder columns of W to correspond to original order of columns in X
    W <- W[, order(qr_X$pivot)]
  } else W <- crossprod(S2_Y_eigen$vectors, W1)

  # if requested, scale generalized kurtosis values
  if (scale_gen_kurtosis) gen_kurtosis <- gen_kurtosis / prod(gen_kurtosis)^(1/p)

  # fix the signs in matrix W of coefficients
  if (fix_signs == "scores") {

    # the condition is phrased so that the last part is only evaluated when
    # necessary (and it is guaranteed that T1_X and T2_X actually exist)
    if (center && !(missing_T2_X || isTRUE(all.equal(T1_X, T2_X)))) {
      # compute generalized skewness values of each component
      T1_Z <- T1_X %*% W
      T2_Z <- T2_X %*% W
      gen_skewness <- as.vector(T1_Z - T2_Z)
    } else {
      # compute scores for initial matrix W of coefficients
      Z <- if (center) tcrossprod(X_centered, W) else tcrossprod(X, W)
      # compute skewness values of each component
      gen_skewness <- colMeans(Z) - apply(Z, 2L, median)
    }

    # compute signs of (generalized) skewness values for each component
    skewness_signs <- ifelse(gen_skewness >= 0, 1, -1)
    # fix signs in W so that generalized skewness values are positive
    gen_skewness <- skewness_signs * gen_skewness
    W_final <- sweep(W, 1L, skewness_signs, "*")

  } else {
    # fix signs in W so that the maximum element per row is positive
    # and that each row has norm 1
    row_signs <- apply(W, 1L, .sign.max)
    row_norms <- sqrt(rowSums(W^2))
    W_final <- sweep(W, 1L, row_norms * row_signs, "/")
    # we don't have (generalized) skewness values in this case
    gen_skewness <- NULL
  }

  # compute the component scores
  if (center) Z_final <- tcrossprod(X_centered, W_final)
  else Z_final <- tcrossprod(X, W_final)

  # set names for different parts of the output
  IC_names <- paste("IC", seq_len(p), sep = ".")
  names(gen_kurtosis) <- IC_names
  dimnames(W_final) <- list(IC_names, colnames(X))
  dimnames(Z_final) <- list(rownames(X), IC_names)
  if (!is.null(gen_skewness)) names(gen_skewness) <- IC_names


  # construct object to be returned
  res <- list(gen_kurtosis = gen_kurtosis, W = W_final, scores = Z_final,
              gen_skewness = gen_skewness, S1_label = S1_label,
              S2_label = S2_label, S1_args = S1_args, S2_args = S2_args,
              QR = QR, whiten = whiten, center = center,
              scale_gen_kurtosis = scale_gen_kurtosis, fix_signs = fix_signs)
  class(res) <- "ICS"
  res

}


# Methods for class "ICS" -----

#' @method summary ICS
#' @export
summary.ICS <- function(object, ...) {
  # currently doesn't do anything but add a subclass
  class(object) <- c("summary_ICS", class(object))
  object
}

#' @method coef ICS
#' @importFrom stats coef
#' @export
coef.ICS <- function(object, select = NULL, drop = FALSE, index = NULL, ...) {
  # back-compatibility check
  if (missing(select) && !missing(index)) {
    warning("argument 'index' is deprecated, use 'select' instead")
    select <- index
  }
  # extract coefficient matrix
  W <- object$W
  # check if we have argument of components to return
  if (!is.null(select)) {
    # check argument specifying components
    if (check_undefined(select, max = nrow(W), names = rownames(W))) {
      stop("undefined components selected")
    }
    # select components
    W <- W[select, , drop = drop]
  }
  # return coefficient matrix for selected components
  W
}

#' @method fitted ICS
#' @export
fitted.ICS <- function(object, select = NULL, index = NULL, ...) {
  # back-compatibility check
  if (missing(select) && !missing(index)) {
    warning("argument 'index' is deprecated, use 'select' instead")
    select <- index
  }
  # initializations
  scores <- object$scores
  p <- ncol(scores)
  if (is.null(select)) select <- seq_len(p)
  else if (check_insufficient(select, target = 1L)) {
    stop("no components selected")
  } else if (check_undefined(select, max = p, names = colnames(scores))) {
    stop("undefined components selected")
  }
  # compute reconstructions from selected components
  # ('drop = FALSE' preserves row and column names)
  W_inverse <- solve(object$W)
  tcrossprod(scores[, select, drop = FALSE], W_inverse[, select, drop = FALSE])
}

#' @export
components <- function(object, ...) UseMethod("components")

#' @export
components.ICS <- function(object, select = NULL, drop = FALSE, index = NULL,
                           ...) {
  # back-compatibility check
  if (missing(select) && !missing(index)) {
    warning("argument 'index' is deprecated, use 'select' instead")
    select <- index
  }
  # extract scores
  scores <- object$scores
  # check if we have argument of components to return
  if (!is.null(select)) {
    # check argument specifying components
    if (check_undefined(select, max = ncol(scores), names = colnames(scores))) {
      stop("undefined components selected")
    }
    # select components
    scores <- scores[, select, drop = drop]
  }
  # return matrix of scores for selected components
  scores
}

#' @export
gen_kurtosis <- function(object, ...) UseMethod("gen_kurtosis")

#' @export
gen_kurtosis.ICS <- function(object, select = NULL, index = NULL, ...) {
  # back-compatibility check
  if (missing(select) && !missing(index)) {
    warning("argument 'index' is deprecated, use 'select' instead")
    select <- index
  }
  # extract generalized kurtosis values
  gen_kurtosis <- object$gen_kurtosis
  # check if we have argument of components to return
  if (!is.null(select)) {
    # check argument specifying components
    if (check_undefined(select, max = length(gen_kurtosis),
                        names = names(gen_kurtosis))) {
      stop("undefined components selected")
    }
    # select components
    gen_kurtosis <- gen_kurtosis[select]
  }
  # return generalized kurtosis values for selected components
  gen_kurtosis
}

#' @method plot ICS
#' @export
#' @importFrom graphics pairs
plot.ICS <- function(x, select = NULL, index = NULL, ...) {
  # back-compatibility check
  if (missing(select) && !missing(index)) {
    warning("argument 'index' is deprecated, use 'select' instead")
    select <- index
  }
  # initializations
  scores <- x$scores
  p <- ncol(scores)
  # create scatterplot matrix
  if (is.null(select)) {
    # not specified which components to plot, use defaults
    if (p <= 6L) pairs(scores, ...)
    else pairs(scores[, c(1:3, p-2:0)], ...)
  } else {
    # check argument specifying components
    if (check_insufficient(select, target = 2L)) {
      stop("'select' must specify at least two components")
    } else if (check_undefined(select, max = p, names = colnames(scores))) {
      stop("undefined components selected")
    }
    # create scatterplot matrix of selected components
    pairs(scores[, select], ...)
  }
}

#' @method print ICS
#' @export
print.ICS <- function(x, info = FALSE, digits = 4L, ...){
  # print information on scatter matrices
  # Andreas: I uncommented the parameter values of the scatters because this
  #          creates problems if the arguments are not single character,
  #          numeric, or logical values
  cat("\nICS based on two scatter matrices")
  cat("\nS1:", x$S1_label)
  # sapply(seq_along(x$S1_args), function(i) {
  #   cat("\n", paste0(" ", names(x$S1_args)[i], ":"), x$S1_args[[i]])
  # })
  cat("\nS2:", x$S2_label)
  # sapply(seq_along(x$S2_args), function(i) {
  #   cat("\n", paste0(" ", names(x$S2_args)[i], ":"), x$S2_args[[i]])
  # })
  # if requested, print information on additional arguments
  if (isTRUE(info)) {
    cat("\n\nInformation on the algorithm:")
    cat("\nQR:", x$QR)
    cat("\nwhiten:", x$whiten)
    cat("\nscale_gen_kurtosis:", x$scale_gen_kurtosis)
    cat("\nfix_signs:", x$fix_signs)
    cat("\ncenter:", x$center)
  }
  # print generalized kurtosis measures and coefficient matrix
  cat("\n\nThe generalized kurtosis measures of the components are:\n")
  # print(formatC(x$gen_kurtosis, digits = digits, format = "f"), quote = FALSE)
  print(x$gen_kurtosis, digits = digits, ...)
  cat("\nThe coefficient matrix of the linear transformation is:\n")
  # print(formatC(x$W, digits = digits, format = "f", flag = " "), quote = FALSE)
  print(x$W, digits = digits, ...)
  # return object invisibly
  invisible(x)
}

#' @method print summary_ICS
#' @export
print.summary_ICS <- function(x, info = TRUE, digits = 4L, ...) {
  # call method for class "ICS" with default for printing additional information
  print.ICS(x, info = info, digits = digits, ...)
}


# Internal functions -----


# check if argument 'select' specifies undefined components
check_undefined <- function(select, max, names) {
  (is.numeric(select) && length(select) > 0L &&
     (min(select) < 1L || max(select) > max)) ||
    (is.character(select) && !all(select %in% names))
}

# check if argument 'select' specifies an insufficient number of components
check_insufficient <- function(select, target = 2L) {
  ((is.numeric(select) || is.character(select)) && length(select) < target) ||
    (is.logical(select) && sum(select) < target)
}


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
# applied, so we need a different way to get the label for those scatters
get_cov4_label <- function() "COV4"
get_covW_label <- function() "COVW"
get_covAxis_label <- function() "COVAxis"
