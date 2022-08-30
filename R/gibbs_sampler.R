#' Continuous Response Predictive Probability.
#'
#' \code{gibbs.sampler} runs the estimation of the continuous response ANOVA predictive probability.
#'
#' The ANOVA model includes main effects and two-way interactions. Priors
#' on model parameters are assumed to be independent of each other; beta is
#' then defined as the set of model parameters, which is multivariate normal.
#'
#' @param X Design matrix for the test (matrix of indicator functions defining
#'   which model parameters are active in each test event).
#' @param Y A vector of the responses from the test.
#' @param n.seen Number of test events already observed (i.e., the number of rows
#'   of the design matrix, X, that have been observed).
#' @param beta.mean Mean vector of the multivariate normal distribution
#'   (or the mean of each of the priors on the model parameters), ordered
#'   the same as the columns of the design matrix, X. It also serves as
#'   the initialization for the model parameters.
#' @param beta.precision Precisions of the multivariate normal distribution
#'   (precision of each of the priors on the model parameters), corresponding
#'   to the beta.mean values.
#' @param precision.a Hyperparameter alpha for gamma prior on the precision of the ANOVA model, tau.
#' @param precision.b Hyperparameter beta for gamma prior on the precision of the ANOVA model, tau.
#' @param n.sim Number of non-conditional posterior draws (i.e., number of
#'   draws that will be returned to the user from the function after burn-in draws for the
#'   non-conditional draws are removed).
#' @param y.burnin Number of burn-in samples for the non-conditional posterior.
#' @param b.sim Number of conditional posterior draws used in analysis
#'   for each non-conditional draw.
#' @param b.burnin Number of burn-in samples for the conditional posterior.
#' @param phi.0 Threshold value the parameter of interest (BMM) must obtain
#'   (i.e., BBM > \code{phi.0}).
#' @param theta.t Certainty threshold for the conditional posterior probability of the
#'   parameter of interest (the Bayesian mission mean, "BMM") obtaining \code{phi.0}
#'   (i.e., BMM > \code{phi.0}) that the conditional posterior probability must obtain
#'   (the certainty threshold for conditional P(BMM > \code{phi.0}) ) must obtain for the
#'   question of interest to be evaluated as successfully passing the test.
#' @param prob Matrix or dataframe of the "likelihood of encountering" (or probability of seeing a
#'   factor level); it is a two column matrix (or dataframe), where the first column identifies the
#'   factor numerically and the second column defines the probability of seeing each
#'   factor level.
#' @param factor.no.2way Optional vector of model parameters (as defined by prob)
#'   that are not incorporated in the two way interactions for the model.
#' @param colnames.pick Optional vector of model parameter names in the same order
#'   as in the design matrix to label the returned dataframe columns.
#' @return Returns a list with three elements:
#'   \describe{
#'     \item{\code{pp}}{The predicted probability of the test ending in a successful evaluation of the question of interest}
#'     \item{\code{posterior}}{The full dataframe of non-conditional posterior draws}
#'     \item{\code{indicator}}{The vector of test success results for each posterior draw}
#'   }
#'   Printing the result object will display the predicted probability result.
#' @importFrom stats aggregate rgamma rmultinom rnorm
#' @export
gibbs.sampler <- function(X, Y, n.seen, beta.mean, beta.precision, precision.a, precision.b,
                          n.sim, y.burnin, b.sim, b.burnin,
                          phi.0, theta.t, prob, factor.no.2way = NA, colnames.pick = NA) {

  # Check X is not NULL
  if(is.null(X) | is.na(X)) { stop("X must not be NULL or NA.") }

  ###### Run predictive function or posterior function, depending on n.seen vs design matrix

  if(n.seen == nrow(X)) {
    warning(
      "Since n.seen is the same as the number of rows of the design matrix,
       note you are calculating the posterior probability and not the predictive probability.",
      immediate. = TRUE
    )
    output <- gibbs.sampler.posterior(
      X, Y, beta.mean, beta.precision,
      precision.a, precision.b, b.sim, b.burnin,
      phi.0, prob, factor.no.2way, colnames.pick
    )
  } else {
    output <- gibbs.sampler.predictive(
      X, Y, n.seen, beta.mean, beta.precision, precision.a, precision.b,
      n.sim, y.burnin, b.sim, b.burnin,
      phi.0, theta.t, prob, factor.no.2way, colnames.pick
    )
  }

  return(output)
}
