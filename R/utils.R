#' Print method for ContRespPP class
#'
#' @param x A list of class ContRespPP
#' @param ... Other arguments passed to \code{print}
#'
#' @return Prints the predictive or posterior probability result from ContRespPP class
#' @export
print.ContRespPP <- function(x, ...) {
  # Define printing method for predictive prob results
  if(! is.null(x$pp)){
    cat("PP Result: ", x$pp)
  }
  # Define printing method for posterior prob results
  if(! is.null(x$probability)){
    cat("Posterior Probability: ", x$probability)
  }
  invisible(x)
}

#' Summary method for ContRespPP class
#'
#' @param object A list of class ContRespPP
#' @param ... Other arguments passed to \code{summary}
#'
#' @return Prints a basic summary of posterior distribution from ContRespPP class
#' @export
summary.ContRespPP <- function(object, ...) {
  # Define printing method for results
  print(summary(object$posterior))
  invisible(object)
}

parameter.data <- function(i, Y, X, full.design) {
  # This function subsets data based on parameter / new responses
  full.data.set <- cbind(Y, X)
  # only keep data that involves the parameter (i)
  subset.data <- full.data.set[(full.data.set[, i] == 1), ]
  #remove p column and parameter column (i)
  subset.data <- subset.data[, -i]
}
