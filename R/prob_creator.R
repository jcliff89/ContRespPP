#' Probability Matrix Creator.
#'
#' \code{prob.creator} creates a matrix of probabilities of encountering factor levels to support
#' construction of input for `gibbs.sampler`, which uses the Bayesian Mission Mean approach to analysis.
#'
#' The ANOVA model includes main effects and two-way interactions. Priors
#' on model parameters are assumed to be independent of each other; beta is
#' then defined as the set of model parameters, which is multivariate normal.
#'
#' @param num.factors Number of factors in the model (e.g., factor \eqn{\alpha_i} with
#'   levels i=1,2 is 1 factor). Input for `num.factors` should be a single number.
#' @param num.factor.levels Number of levels for each factor (e.g., factor \eqn{\alpha_i}
#'   with i=1,2 has 2 levels). Input for `num.factor.levels` may be a vector, a matrix, or a dataframe.
#' @param likelihood.encountering The probability of seeing each level of each factor (e.g., if the factor levels for
#'   \eqn{\alpha_i} are equally likely, then the probabilities would be c(1/2, 1/2)).
#'   The probabilities for each factor should sum to one. Input for `likelihood.encountering` may be a vector, a matrix,
#'   or a dataframe.
#' @param print.result Displays final probability matrix.
#' @return Returns a matrix with two columns, one with the factor number and the other with the likelihoods of encountering.
#' @export
#'
#'

prob.creator <- function(num.factors, num.factor.levels, likelihood.encountering, print.result = FALSE) {

  # Convert to matrix if needed
  if(! any(class(likelihood.encountering) == "matrix")){
    likelihood.encountering <- as.matrix(likelihood.encountering, ncol = 1)
  }
  if(! any(class(num.factor.levels) == "matrix")){
    num.factor.levels <- as.matrix(num.factor.levels, ncol = 1)
  }

  # ERROR: length of num.factor.levels should be equal to num.factors
  if(nrow(num.factor.levels) != num.factors) {
    stop("Number of factor levels provided do not match number of factors provided.")
  }

  # ERROR: length of likelihood.encountering should be equal to sum of num.factor.levels
  if(length(likelihood.encountering) != sum(num.factor.levels)){
    stop("Total num.factor.levels does not match length of likelihood.encountering.")
  }

  # Construct prob matrix
  prob <- matrix(NA, nrow = nrow(likelihood.encountering), ncol = 2)
  prob[,2] <- likelihood.encountering
  levels.new <- NULL
  for(i in 1:num.factors){
    levels.new <- c(levels.new, rep(i, num.factor.levels[i, 1]))
  }
  prob[,1] <- levels.new

  # ERROR: prob object probability values should sum to 1 within groups
  if(! all(aggregate(prob[,2], by = list(prob[,1]), sum)[,2] == 1)) {
    prob_sumNot1 <- aggregate(prob[,2], by = list(prob[,1]), sum)[aggregate(prob[,2], by = list(prob[,1]), sum)[,2] != 1,]
    colnames(prob_sumNot1) <- c("Factor", "Sum")
    prob_sumNot1$Label <- paste("Factor", prob_sumNot1$Factor, "sums to", prob_sumNot1$Sum)
    prob_sumNot1_errorMsg <- paste(
      "prob object contains factor levels with probability values that do not sum to 1: ",
      paste(prob_sumNot1$Label, collapse = ", ")
    )
    stop(prob_sumNot1_errorMsg)
  }

  colnames(prob) <- c("Factor", "Likelihood of Encountering")

  if(print.result) { print(prob) }

  return(prob)
}

