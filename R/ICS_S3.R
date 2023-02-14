
#' @export
ICS_cov <- function(x) {
  # compute center and scatter estimates
  out <- list(location = colMeans(x), scatter = cov(x), label = "COV")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}


#' @export
ICS <- function(X, S1, S2, ...) UseMethod("ICS", S2)

#' @export
ICS.function <- function(X, S1 = cov, S2 = cov4, S1_args = list(),
                         S2_args = list(), whiten = TRUE, qr = NULL,
                         ...) {
  # initializations
  whiten <- isTRUE(whiten)
  if (is.null(qr)) {
    # use algorithm based on QR decomposition by default if we have scatter
    # pair COV-COV4 and we should not whiten the data with respect to S1 before
    # computing S2
    qr <- isTRUE(identical(S1, cov) && identical(S2, cov4) && !whiten)
  }
  # perform ICS
  if (qr) {
    # TODO: perform algorithm based on QR decomposition
  } else {
    # obtain first scatter matrix
    if (is.function(S1)) S1 <- get_scatter(X, fun = S1, args = S1_args)
    else S1 <- to_ICS_scatter(S1)
    # whiten the data if requested
    if (whiten) {
      # TODO
    }
    # obtain second scatter matrix
    S2 <- get_scatter(X, fun = S2, args = S2_args)
    # TODO: can we call the default method now?
  }
}

#' @export
# center ........ logical indicating whether to center the data prior to ICS
# standardize ... character string specifying how to standardize the ICS
#                 coordinates
ICS.default <- function(X, S1, S2, center = FALSE,
                        standardize = c("scores", "eigenvalues"),
                        ...) {
  # initializations
  center <- isTRUE(center)
  standardize <- match.arg(standardize)
  # convert scatter matrices to class "ICS_scatter"
  # (Andreas: in principle, we could allow that the first scatter matrix is a
  # function, but I'm not sure whether this is useful)
  S1 <- to_ICS_scatter(S1)
  S2 <- to_ICS_scatter(S2)
  # centering only makes sense if the scatter objects have location components
  if (center && (is.null(S1$location) || is.null(S2$location))) {
      center <- FALSE
      warning("'S1' and 'S2' need to have a location component for centering ",
              "the data; proceeding without centering")
  }
  # TODO
}


## FIXME: default labels in the internal functions below don't work properly
##        when called from function ICS(). We need to pass down the correct
##        environment so that deparse(substitute()) works as intended.


## internal function to apply a scatter function to the data matrix
# X ......... data matrix
# fun ....... function to compute a scatter matrix
# args ...... list of additional arguments to be passed to the function
# convert ... logical indicating whether the scatter matrix should be converted
#             to class "ICS_scatter". This is useful for the algorithm based on
#             the QR decomposition, which needs functions that return just the
#             scatter matrix.
get_scatter <- function(X, fun = cov, args = list(), convert = TRUE) {
  if (length(args) == 0) scatter <- fun(X)
  else {
    # there should be a more efficient way of doing this, for example via call()
    args <- c(list(X), args)
    scatter <- do.call(fun, args)
  }
  # if requested, convert to class "ICS_scatter": use the function name as
  # label if the function does not return an object of this class already
  if (convert) {
    default_label = deparse(substitute(fun))
    scatter <- to_ICS_scatter(scatter, label = default_label)
  }
  # return scatter matrix
  scatter
}


## internal function to convert scatter matrices to the required class

to_ICS_scatter <- function(object, ...) UseMethod("to_ICS_scatter")

to_ICS_scatter.matrix <- function(object, label = NULL, ...) {
  # if not supplied, use object name as label
  if (is.null(label)) label <- deparse(substitute(object))
  # convert to class "ICS_scatter" with empty location estimate
  out <- list(location = NULL, scatter = object, label = label)
  class(out) <- "ICS_scatter"
  out
}

to_ICS_scatter.list <- function(object, label = NULL) {
  # check that list has a component 'scatter'
  scatter <- object$scatter
  if (is.null(scatter)) {
    stop("list should have components 'scatter' (required) as well as ",
         "'location' and 'label' (optional)")
  }
  # check if there already  is a component 'label', otherwise use the supplied
  # label, or use the object name if no label is supplied
  if (!is.null(object$label)) label <- object$label
  if (is.null(label)) label <- deparse(substitute(object))
  # convert to class "ICS_scatter"
  out <- list(location = object$location, scatter = scatter, label = label)
  class(out) <- "ICS_scatter"
  out
}

to_ICS_scatter.ICS_scatter <- function(object, ...) object
