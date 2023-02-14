
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
    S1_label <- deparse(substitute(S1))
    if (is.function(S1)) {
      S1 <- get_scatter(X, fun = S1, args = S1_args, label = S1_label)
    } else S1 <- to_ICS_scatter(S1, label = S1_label)
    # whiten the data if requested
    if (whiten) {
      # TODO
    }
    # obtain second scatter matrix
    S2_label <- deparse(substitute(S2))
    S2 <- get_scatter(X, fun = S2, args = S2_args, label = S2_label)
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
  S1 <- to_ICS_scatter(S1, label = deparse(substitute(S1)))
  S2 <- to_ICS_scatter(S2, label = deparse(substitute(S2)))
  # centering only makes sense if the scatter objects have location components
  if (center && (is.null(S1$location) || is.null(S2$location))) {
      center <- FALSE
      warning("'S1' and 'S2' need to have a location component for centering ",
              "the data; proceeding without centering")
  }
  # TODO
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
