parameter.data <- function(i, Y, X, full.design) {
  # This function subsets data based on parameter / new reponses
  full.data.set <- cbind(Y, X)
  # only keep data that involves the parameter (i)
  subset.data <- full.data.set[(full.data.set[, i] == 1), ]
  #remove p column and parameter column (i)
  subset.data <- subset.data[, -i]
}

print.ContRespPP <- function(x) {
  # Define printing method for results
  cat("PP Result: ", x$pp)
  invisible(x)
}
