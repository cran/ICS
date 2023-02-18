

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


S2_Y <- function(X, S1, S2, ...) UseMethod("S2_Y", S2)

S2_Y.function <- function(X, S1.X, S2, S2_args = list(),
                          whiten = NULL,
                          QR = NULL, S1_label, S2_label, ...) {
  # Initialization
  S2.X <- S2

  if (isTRUE(QR)) {
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
    # convert scatter matrices to class "ICS_scatter"
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
    # convert scatter matrices to class "ICS_scatter"
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
  list(X, S2.X = S2.X, S2.Y = S2.Y, B1 = B1, R = R, X_Z = X_Z,
       ...)
}


S2_Y.matrix <- function(X, S1.X, S2, S1_label, S2_label, ...) {

  # perform 1st step of ICS
  # compute B1 = S1.X^-1/2
  S1.X.eigen <- eigen(S1.X$scatter, symmetric=TRUE)
  B1 <- S1.X.eigen$vectors %*% tcrossprod(diag(S1.X.eigen$values^(-0.5)), S1.X.eigen$vectors)

  # obtain second scatter matrix
  # or compute S1.X^-1/2*S2.Y*S1.X^-1/2
  # convert scatter matrices to class "ICS_scatter"
  S2.X <- to_ICS_scatter(S2, label = S2_label)
  S2.Y <- to_ICS_scatter(B1 %*% S2.X$scatter %*% B1, label = S2_label)

  list(S2.X = S2.X, S2.Y = S2.Y, B1 = B1, X_Z = X, ...)
}

S2_Y.ICS_scatter <- S2_Y.matrix




#' Title
#'
#' @param X
#' @param S1
#' @param S2
#' @param S1_args
#' @param S2_args
#' @param center logical indicating whether to center the ICS coordinates (scores)
#' @param QR
#' @param fix_signs character string specifying how to fix_signs the ICS
#                 coordinates, either 'scores' or 'eigenvectors'
#' @param scale_lambda
#' @param na.action
#' @param whiten
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ICS <- function(X, S1 = ICS_cov, S2 = ICS_cov4, S1_args = list(), S2_args =  list(),
                center = FALSE,
                QR = NULL,
                fix_signs = c("scores", "W"),
                scale_lambda = FALSE,
                whiten = NULL,
                na.action = na.fail,
                ...) {
  # initializations - data matrix
  X <- na.action(X)
  X <- as.matrix(X)
  p <- dim(X)[2]
  if (p < 2)
    stop("'X' must be at least bivariate")

  # obtain first scatter matrix
  S1_label <- deparse(substitute(S1))
  S2_label <- deparse(substitute(S2))

  # convert scatter matrices to class "ICS_scatter"
  if (is.function(S1)) {
    S1.X <- get_scatter(X, fun = S1, args = S1_args, label = S1_label)
  } else S1.X <- to_ICS_scatter(S1, label = S1_label)


  # compute S2_Y
  # determine which should be the best default value of whiten
  if(is.null(whiten) & is.function(S2)){
    whiten <- TRUE
  }
  # determine which should be the best default value of QR
  if(is.null(QR) & is.function(S2)){
    if (S1_label %in% c("cov", "ICS_cov") &
        S2_label %in% c("covW", "ICS_covW")){
      QR <- TRUE

    }
  }
  S2_Y <- S2_Y(X, S1.X = S1.X, S2 = S2, S2_args = S2_args,
               whiten = whiten, QR = QR, S1_label, S2_label, ... )
  X_Z <- S2_Y$X_Z
  center <- isTRUE(center)
  fix_signs <- match.arg(fix_signs)

  # decomposition of S2
  S2.Y.eigen <- eigen(S2_Y$S2.Y$scatter, symmetric=TRUE)
  U2 <- S2.Y.eigen$vectors
  lambda <- S2.Y.eigen$values
  if (isTRUE(QR)) {
    # Eigenvectors by rows
    Rinv <- qr.solve(S2_Y$R)
    B <- t(U2) %*% t(Rinv)

  }else{
    B <- crossprod(U2, S2_Y$B1)
  }

  # choosing the signs of B
  # centering only makes sense if the scatter objects have location components
  if (center){
    if (is.function(S2)){
      if(is.null(S1.X$location) & is.function(S2)) {
        center <- FALSE
        warning("'S1' need to have a location component for centering ",
                "the data; proceeding without centering")
      }
    }else{
      if(is.null(S1.X$location) ||
         is.null(S2_Y$S2.X$location)){
        center <- FALSE
        warning("'S1' and 'S2' need to have a location component for centering ",
                "the data; proceeding without centering")
      }
    }
  }



  if (center){
    if (is.function(S2)) {
      S2.X <- get_scatter(X, fun = S2, args = S2_args)
    }else{
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
    stdB = "Z"
  }else{
    if (scale_lambda == TRUE) lambda <- lambda/prod(lambda)^(1/p)
    if (fix_signs == "W") {
      row.signs <- apply(B, 1, .sign.max)
      row.norms <- sqrt(rowSums((B)^2))
      B.res <- sweep(B, 1, row.norms * row.signs, "/")
      gamma <- NULL
    }
    if (fix_signs == "scores" & center == FALSE) {
      Z1 <- tcrossprod(X_Z, B)
      gamma <- colMeans(Z1) - apply(Z1, 2, median)
      skew.signs <- ifelse(gamma > 0, 1, -1)
      gamma <- as.vector(skew.signs*gamma)
      B.res <- sweep(B, 1, skew.signs, "*")
    }
  }

  # compute the scores
  Z <- as.data.frame(tcrossprod(X_Z, B.res))
  names(Z) <- paste(rep("IC", p), 1:p, sep = ".")

  if (is.null(colnames(X)) == TRUE)
    names.X <- paste(rep("X", p), 1:p, sep = ".")
  else names.X <- colnames(X)

  res <- list(lambda = lambda, W = B.res,
              scores = Z,
              gamma = gamma,
              S1_label = S1_label,
              S2_label = S2_label,
              S1_args = S1_args,
              S2_args = S2_args,
              center = center,
              QR = ifelse(is.null(QR), FALSE, QR),
              fix_signs = fix_signs,
              scale_lambda = scale_lambda,
              whiten = ifelse(is.null(whiten), FALSE, whiten),
              X = X) # it is required if we want to plot X.

  class(res) <- "ICS"
  res

}




#' @method summary ICS
#' @export
summary.ICS <- function(object,digits=4)
{
  cat("\nICS with the following parameters: \n")
  cat("S1:", object$S1_label)
  if (length(object$S1_args)>0){
    sapply(1:length(object$S1_args), function(i)
      cat("\n", paste0(names(object$S1_args)[i],":"), object$S1_args[[i]]))
  }
  if (length(object$S2_args)>0){
    cat("\nS2:", object$S2_label)
    sapply(1:length(object$S2_args), function(i)
      cat("\n", paste0(names(object$S2_args)[i],":"), object$S2_args[[i]]))
  }
  cat("\nQR:", object$QR)
  cat("\nwhiten:", object$whiten)
  cat("\nscale_lambda:", object$scale_lambda)
  cat("\nfix_signs:", object$fix_signs)
  cat("\ncenter:", object$center)
  cat("\n")
  cat("\nThe generalized kurtosis measures (lambda) of the components are:\n")
  print(format(round(object$lambda,digits)),quote=F)
  cat("\n")
  cat("\nThe W matrix is:\n")
  print(round(object$W,digits))
  invisible(object)
}

#' @method coef ICS
#' @export
coef.ICS <- function(object){
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
print.ICS <- function(object){
  tmp <- list(lambda=object$lambda,
              W=object$W)
  print(tmp, quote = FALSE)
  invisible(tmp)
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
