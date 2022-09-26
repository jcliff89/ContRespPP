#' Continuous Response Posterior Probability.
#'
#' \code{gibbs.sampler.posterior.rjags} obtains the ContRespPP method for obtaining the posterior (i.e., all
#' observations have been seen, and this reverts to a traditional Bayesian analysis) using rjags to
#' generate samples instead of base R functions, drastically decreasing the computational time required
#' to obtain predictive draws.
#'
#' @param X Design matrix for the test (matrix of indicator functions defining
#'   which model parameters are active in each test event).
#' @param Y A vector of the responses from the test.
#' @param beta.mean Mean vector of the multivariate normal distribution
#'   (or the mean of each of the priors on the model parameters), ordered
#'   the same as the columns of the design matrix, X. It also serves as
#'   the initialization for the model parameters.
#' @param beta.precision Precisions of the multivariate normal distribution
#'   (precision of each of the priors on the model parameters), corresponding
#'   to the beta.mean values.
#' @param shape Hyperparameter alpha for gamma prior on the precision of the ANOVA model, tau.
#' @param rate Hyperparameter beta for gamma prior on the precision of the ANOVA model, tau.
#' @param b.sim Number of conditional posterior draws used in analysis
#'   for each non-conditional draw.
#' @param b.burnin Number of burn-in samples for the conditional posterior.
#' @param phi.0 Threshold value the parameter of interest (BMM) must obtain
#'   (i.e., BBM > \code{phi.0}).
#' @param prob Matrix or dataframe of the "likelihood of encountering" (or probability of seeing a
#'   factor level); it is a two column matrix (or dataframe), where the first column identifies the
#'   factor numerically and the second column defines the probability of seeing each
#'   factor level.
#' @param factor.no.2way Optional vector of model parameters (as defined by prob)
#'   that are not incorporated in the two way interactions for the model.
#' @param colnames.pick Optional vector of model parameter names in the same order
#'   as in the design matrix to label the returned dataframe columns.
#' @param seed Optional selection which will create a reproducible result from the function.
#' @param verbose Allows suppression of sampler progress printing in console.
#' @return Returns a list with three elements:
#'   \describe{
#'     \item{\code{pp}}{This value will be NA since this function only calculates the posterior}
#'     \item{\code{posterior}}{The full dataframe of non-conditional posterior draws}
#'     \item{\code{indicator}}{This value will be NA since this function only calculates the posterior}
#'   }
#'   Printing the result object will display the predicted probability result.
#' @importFrom stats aggregate rgamma rmultinom rnorm
#' @export
gibbs.sampler.posterior.rjags <- function(X, Y, beta.mean, beta.precision, shape, rate,
                                          b.sim, b.burnin, phi.0, prob, factor.no.2way = NA, colnames.pick = NA,
                                          seed = NA, verbose = TRUE) {

  # Convert non-matrix inputs to matrix for remainder of function to run smoothly
  if(any(class(X) == "data.frame")){ X <- as.matrix(X) }
  if(any(class(Y) == "data.frame")){ Y <- as.matrix(Y) }
  if(! any(class(Y) == "matrix")){ Y <- matrix(Y, ncol = 1) }
  if(! any(class(beta.mean) == "matrix")){ beta.mean <- matrix(beta.mean, ncol = 1) }
  if(! any(class(beta.precision) == "matrix")){ beta.precision <- matrix(beta.precision, ncol = 1) }
  if(any(class(prob) == "data.frame")){ prob <- as.matrix(prob) }

  ##### Validation checks to make sure function will run properly

  # ERROR: Numeric inputs expected
  if(! is.numeric(X)) { stop("X must be a numeric type.") }
  if(! is.numeric(Y)) { stop("Y must be a numeric type.") }
  if(! is.numeric(beta.mean)) { stop("beta.mean must be a numeric type.") }
  if(! is.numeric(shape)) { stop("shape must be a numeric type.") }
  if(! is.numeric(rate)) { stop("rate must be a numeric type.") }
  if(! is.numeric(b.sim)) { stop("b.sim must be a numeric type.") }
  if(! is.numeric(b.burnin)) { stop("b.burnin must be a numeric type.") }
  if(! is.numeric(phi.0)) { stop("phi.0 must be a numeric type.") }
  if(! is.numeric(prob)) { stop("prob must be a numeric type.") }
  if(! is.na(seed) & ! is.numeric(seed)) { stop("seed must be a numeric type.") }

  # ERROR: Design matrix should be same size as number of priors
  if( ncol(X) != length(beta.mean) ) {
    stop("Design matrix X does not have the same number of parameters (columns) as model
         initialization vector beta.mean.")
  }

  # ERROR: prob object probability values should all be 0-1 and should sum to 1 within groups
  if(! all(prob[,2] >= 0 & prob[,2] <= 1)) {
    stop("prob object contains probability values not within [0,1]")
  }

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

  # ERROR: Design matrix values should be all 0/1
  if(! all(X == 0 | X == 1)) { stop("Not all values in design matrix X are 0 or 1") }

  # ERROR: Design matrix should have intercept
  if (! (all(X[,1] == 1))) { stop("Model must have an intercept, where the first column of the design matrix contains all 1 values; not all values in the first column of the design matrix are 1") }

  # ERROR: gamma prior on tau parameters must be positive
  if(! (shape > 0)) { stop("Gamma parameters must be greater than 0 (shape is not)") }
  if(! (rate > 0)) { stop("Gamma parameters must be greater than 0 (rate is not)") }


  ##### Warnings to user when function will run alright but something is strange

  if(b.burnin > b.sim) {
    warning(
      "Burn-in draws (b.burnin) higher than non-conditional posterior draws (b.sim).",
      immediate. = TRUE
    )
  }

  ##### The user provides b.sim, expecting those to be the number of posterior draws returned; however code uses
  # b.sim as the total number sampled. So here we re-define b.sim to be b.sim provided by the user plus b.burnin.
  b.sim <- b.sim + b.burnin

  ##### Create a matrix with the response in the first column, and the remaining columns as the design matrix
  full.design <- cbind(Y, X)

  ##### Get the number of model parameters -- not including tau
  num.param <- length(beta.mean)

  ##### Get column name to assign to results (if not provided, use colnames of design matrix)
  if(is.na(colnames.pick[1])) {colnames.pick <- c(colnames(X), "tau")}



  ############################################################################################################
  #############################################   MISSION SETS   #############################################
  ############################################################################################################

  ##### We need a mission set matrix that is num.param columns (1 for each parameter in the model)
  # and b.sim - b.burnin rows (1 for each posterior draw that we will keep after burn in).
  ##### This grabs the mission for each column of the matrix, e.g. the first one grabs nothing but 1's because
  # the model parameter eta is the intercept essentially. The rest are generated from binomial or multinomial distributions,
  # as defined by the likelihood of encountering the parameter. Once combined, this essentially will create a new design matrix
  # of possible scenarios we could see in the field and evlaute one conditional posterior draw against one of the rows in the mission sets.

  if(is.na(seed) == FALSE){
    set.seed(seed)
  }

  num.missions <- b.sim - b.burnin
  num.maineffects <- length(prob[, 1]) - length(unique(prob[, 1]))

  missions <-
    matrix(NA, ncol = (1 + num.maineffects), nrow = num.missions)
  missions[, 1] <- 1
  for (f in 1:max(prob[, 1])) {
    sample.mission.set <-
      rmultinom(num.missions, 1, prob[prob[, 1] == f, 2])
    last.notna.col <- sum(!is.na(missions)) / num.missions
    ##### This takes the number of rows sampled (again, in case it isn't a two level factor),and puts it in
    # the next columns filled with NA. If there is a two level factor it will take the second row of the multinomial
    # vector and put it in one column. If it is a three-level factor it will take the second and third row of the multinomial vector
    # and put it in one column. So on and so forth.
    missions[, (last.notna.col + 1):(last.notna.col + nrow(sample.mission.set) -
                                       1)] <- t(sample.mission.set)[, 2:ncol(t(sample.mission.set))]
  }

  ##### Label each row by the factor number from the object prob, so that main effects that are not included in the
  # two-way interactions can be removed by their factor number value.
  colnames(missions) <- c(0, prob[duplicated(prob[, 1]), 1])

  ##### Remove the factors not needed for two-way interactions
  if (!is.na(factor.no.2way)) {
    sub.missions <- missions[, colnames(missions) != factor.no.2way]
  }
  ##### Remove the intercept / reference cell number
  sub.missions <- sub.missions[, -1]

  ##### Get the number of two-way interactions needed from the remaining values: (n * (n-1))/2
  num.int <- ((ncol(sub.missions)) * (ncol(sub.missions) - 1)) / 2
  ##### Create a matrix to store the two-way interactions in
  missions.int <- matrix(NA, ncol = num.int, nrow = num.missions)

  ##### Get two-way interactions by multiplying each main effect column (element-wise) by the main effect columns that follow it, one column at a time.
  for (t in 1:(ncol(sub.missions) - 1)) {
    for (u in (t + 1):ncol(sub.missions)) {
      missions.int[, ((sum(!is.na(missions.int)) / num.missions) + 1)] <-
        sub.missions[, t] * sub.missions[, u]
    }
  }

  ##### Create the full mission set
  mission.set <- cbind(missions, missions.int)
  mission.set <- t(mission.set)


  #### Calculate initial value for tau as mean of gamma distribution
  tau.int <- shape / rate


  ##############################################################################################################
  #############################################   Posterior Calc   #############################################
  ##############################################################################################################

  parameters <- c( "beta", "tau" )
  data <- list( "X"=X, "n"=nrow(X), "Y"=Y, "beta.mean"=c(beta.mean), "beta.precision"=c(beta.precision), "m"=length(beta.mean), "shape"=shape, "rate"=rate)
  if(is.na(seed) == FALSE){
    inits <- list( ".RNG.name"= "base::Wichmann-Hill", ".RNG.seed"= seed, beta=c(beta.mean), tau=tau.int )
  } else{
    inits <- list( beta=c(beta.mean), tau=tau.int )
  }

  modelstring <- "model
        	{
        		for( i in 1:n ){
        			Y[i,] ~ dnorm(mu[i], tau)
        			mu[i]<- inprod(beta[],X[i,])
        		}

        	## Define Priors
        		tau ~ dgamma(shape, rate)
        		for( j in 1:m){
          		beta[j]~dnorm(beta.mean[j],beta.precision[j])
        		}
        	}"



  post.m <- rjags::jags.model( data=data, inits=inits, file=textConnection(modelstring), n.chains=1 , quiet = TRUE)
  post.wbi <- rjags::coda.samples(post.m, parameters, n.iter=b.sim, thin=1, n.burn=0, progress.bar = NULL)
  post.wbi <- as.matrix(post.wbi)
  post.wbi <- cbind(post.wbi, NA)


  ####################################################################################################################################
  ####################################################################################################################################
  #############################################								           #############################################
  #############################################   SAVE VARIOUS OUTPUTS / RETURN VALUES   #############################################
  #############################################					       ` 			   #############################################
  ####################################################################################################################################
  ####################################################################################################################################

  ##### Remove Burn in from posterior draws
  post.wobi <- post.wbi[-(1:b.burnin),]

  ##### Calculate mission mean using the conditional posterior draws and mission set matrix. To do this, first row of
  # post.wobi multiplied by first column of mission.set, second row of post.wobi multiplied by second column of
  # mission.set, etc (alternatively, you could see this as the diagonal of the posterior draws matrix multiplied by the
  # mission set matrix). Save the result in the last column of post.wobi.
  for (s in 1:nrow(post.wobi)) {
    post.wobi[s, (num.param + 2)] <-
      post.wobi[s, 1:num.param] %*% mission.set[, s]
  }

  colnames(post.wobi) <- c(colnames.pick, "m")

  ###### Calculate the conditional probability of the test being a success at the end of the test; this is, the number
  # of mission means that were above phi.0 divided by the total number of mission means to get Pr.
  post.prob <- sum(post.wobi[, (num.param + 2)] > phi.0) / nrow(post.wobi)

  if(verbose){cat("\n Simulation complete \n")}

  output <- structure(
    list(
      probability = post.prob,
      posterior = post.wobi
    ),
    class = "ContRespPP"
  )

  return(output)

}
