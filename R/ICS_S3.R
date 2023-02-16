
#' @export
ICS_cov <- function(x) {
  # compute center and scatter estimates
  out <- list(location = colMeans(x), scatter = cov(x), label = "COV")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}

#' @export
ICS_cov4 <- function(x) {
  # compute center and scatter estimates
  out <- list(location = colMeans(x), scatter = cov4(x), label = "COV4")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}

#' @export
ICS_Mean3Cov4 <- function(x) {
  # compute center and scatter estimates
  out <- list(location = mean3(x), scatter = cov4(x), label = "COV4")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}

#' @export
ICS_covW <- function(x, na.action = na.fail, alpha = 1, cf = 1 ) {
  # compute center and scatter estimates
  out <- list(location = colMeans(x), scatter = covW(x, na.action = na.fail, alpha = alpha, cf = cf), label = "COVW")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}



#' @export
ICS <- function(X, S1, S2, ...) UseMethod("ICS", S2)

#' @export
ICS.function <- function(X, S1 = cov, S2 = cov4, S1_args = list(),
                         S2_args = list(), whiten = TRUE, QR = FALSE,                           na.action = na.fail, ...) {
  # initializations
  whiten <- isTRUE(whiten)
  X <- na.action(X)
  X <- as.matrix(X)
  S1_label <- deparse(substitute(S1))
  S2_label <- deparse(substitute(S2))
  S2.X <- S2

  # obtain first scatter matrix
  if (is.function(S1)) {
    S1.X <- get_scatter(X, fun = S1, args = S1_args, label = S1_label)
  } else S1.X <- to_ICS_scatter(S1, label = S1_label)

  # if (QR) {
  #   # use algorithm based on QR decomposition by default if we have scatter
  #   # Aurore: why do we want that?
  #   # pair COV-COV4 and we should not whiten the data with respect to S1 before
  #   # computing S2
  #   QR <- isTRUE(identical(S1, cov) && identical(S2, cov4) && !whiten)
  # }
  # perform ICS
  if (QR) {
    whiten <- FALSE
    if (!(S1_label %in% c("cov", "ICS_cov") &  S2_label %in% c("covW", "ICS_covW"))){

      warning("'QR' is possible only for 'S1' being 'cov' or 'ICS_cov' and ",
              "'S2' being 'covW' or 'ICS_covW'")
      break

    }
    # Aurore: do we want to change cov4 to covW with correct arguments?
    # same for CovAxis?

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
    S2.Y <- to_ICS_scatter(S2_args$cf*1/n*(n-1)*t(Q) %*% diag(d_i^S2_args$alpha) %*% Q, label = S2_label)

    # Reorder the rows and cols of X
    X_Z <-  X_row_per[order(order_rows),qr_X$pivot]
    B1 <- NULL
  } else {
    # compute B1 = S1.X^-1/2
    S1.X.eigen <- eigen(S1.X$scatter, symmetric=TRUE)
    B1 <- S1.X.eigen$vectors %*% tcrossprod(diag(S1.X.eigen$values^(-0.5)), S1.X.eigen$vectors)

    # obtain second scatter matrix
    # whiten the data if requested - only possible if S2 is a function
    if (whiten & is.function(S2)) {
      Y <- X %*% B1
      S2.Y <- get_scatter(Y, fun = S2, args = S2_args, label = S2_label)
    }else{
      # or compute S1.X^-1/2*S2.Y*S1.X^-1/2
      S2.X <- get_scatter(X, fun = S2, args = S2_args)
      S2.Y <- to_ICS_scatter(B1 %*% S2.X$scatter %*% B1, label = S2_label)
    }
    X_Z <- X
    R <- NULL
  }
  # call the default method
  ICS(X, S1 = S1.X, S2 = S2.Y, QR = QR, B1 = B1, R = R, X_Z = X_Z,
      S2.X =  S2.X, S2_args = S2_args,
      ...)
}


#' @export
ICS.matrix <- function(X, S1 = cov, S2 = cov4, S1_args = list(),
                       S2_args = list(), whiten = FALSE, QR = FALSE,
                       na.action = na.fail, ...) {
  # initializations
  X <- na.action(X)
  X <- as.matrix(X)
  S1_label <- deparse(substitute(S1))
  S2_label <- deparse(substitute(S2))

  # obtain first scatter matrix
  if (is.function(S1)) {
    S1.X <- get_scatter(X, fun = S1, args = S1_args, label = S1_label)
  } else S1.X <- to_ICS_scatter(S1, label = S1_label)


  # perform ICS
  # compute B1 = S1.X^-1/2
  S1.X.eigen <- eigen(S1.X$scatter, symmetric=TRUE)
  B1 <- S1.X.eigen$vectors %*% tcrossprod(diag(S1.X.eigen$values^(-0.5)), S1.X.eigen$vectors)

  # obtain second scatter matrix
  # whiten the data if requested - only possible if S2 is a function
  # or compute S1.X^-1/2*S2.Y*S1.X^-1/2
  S2.X <- to_ICS_scatter(S2, label = S2_label)
  S2.Y <- to_ICS_scatter(B1 %*% S2.X$scatter %*% B1, label = S2_label)


  # call the default method
  ICS(X, S1 = S1.X, S2 = S2.Y, QR = QR, B1 = B1, R = NULL, X_Z = X,
      S2.X =  S2.X,
      ...)
}


#' @export
# center ........ logical indicating whether to center the ICS coordinates (scores)
# standardize ... character string specifying how to standardize the ICS
#                 coordinates
ICS.default <- function(X, S1, S2, center = FALSE, QR = FALSE,
                        standardize = c("scores", "eigenvalues"),
                        stdKurt=FALSE,
                        B1 = NULL, R = NULL, X_Z = NULL,  S2.X = NULL,
                        S2_args = list(),
                        ...) {
  # initializations
  p <- dim(X)[2]
  if (p < 2)
    stop("'X' must be at least bivariate")

  center <- isTRUE(center)
  standardize <- match.arg(standardize)
  # convert scatter matrices to class "ICS_scatter"
  # (Andreas: in principle, we could allow that the first scatter matrix is a
  # function, but I'm not sure whether this is useful)
  #S1 <- to_ICS_scatter(S1, label = deparse(substitute(S1)))
  # S2 <- to_ICS_scatter(S2, label = deparse(substitute(S2)))

  # centering only makes sense if the scatter objects have location components
  if (center && (is.null(S1$location) || is.null(S2$location))) {
    center <- FALSE
    warning("'S1' and 'S2' need to have a location component for centering ",
            "the data; proceeding without centering")
  }

  # decomposition of S2
  S2.Y.eigen <- eigen(S2$scatter, symmetric=TRUE)
  U2 <- S2.Y.eigen$vectors
  gKurt <- S2.Y.eigen$values
  if (QR) {
    # Eigenvectors by rows
    Rinv <- qr.solve(R)
    B <- t(U2) %*% t(Rinv)

  }else{
    B <- crossprod(U2, B1)
  }


  # choosing the signs of B
  if (center){
    if (is.function(S2.X)) {
      S2.X <- get_scatter(X, fun = S2.X, args = S2_args)
    }else{
      S2.X <- to_ICS_scatter(S2.X, label = deparse(substitute(S2)))
    }

    T1.Z <- S1$location %*% B
    T2.Z <- S2.X$location %*% B

    gSkew <- T1.Z - T2.Z
    skew.signs <- ifelse(gSkew >= 0, 1, -1)

    B.res <- sweep(B, 1, skew.signs, "*")
    X_Z <- sweep(X_Z, 2, S1$location, "-")
    # What should it be

    standardize <- "scores"
    stdKurt <- FALSE
    stdB = "Z"
  }else{
    if (stdKurt == TRUE) gKurt <- gKurt/prod(gKurt)^(1/p)
    if (standardize == "eigenvalues") {
      stdB = "B"
      row.signs <- apply(B, 1, .sign.max)
      row.norms <- sqrt(rowSums((B)^2))
      B.res <- sweep(B, 1, row.norms * row.signs, "/")
    }
    if (standardize == "scores" & center == FALSE) {
      stdB = "Z"
      Z1 <- tcrossprod(X_Z, B)
      skewness <- colMeans(Z1) - apply(Z1, 2, median)
      skew.signs <- ifelse(skewness > 0, 1, -1)
      B.res <- sweep(B, 1, skew.signs, "*")
    }
  }

  # compute the scores
  Z <- as.data.frame(tcrossprod(X_Z, B.res))
  names(Z) <- paste(rep("IC", p), 1:p, sep = ".")

  if (is.null(colnames(X)) == TRUE)
    names.X <- paste(rep("X", p), 1:p, sep = ".")
  else names.X <- colnames(X)

  res <- list(gKurt = gKurt, UnMix = B.res, S1 = S1, S2 = S2, S1name = S1$label,
              S2name = S2$label, Scores = Z, DataNames = names.X,
              StandardizeB = stdB, StandardizegKurt = stdKurt,
              standardize = standardize)
  class(res) <- "ICS"
  res

}


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


## internal function to convert scatter matrices to the required class

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
