#' Continuous Response Predictive Probability.
#'
#' \code{gibbs.sampler.predictive.jags} does the predictive ContRespPP method using rjags to generate
#' samples instead of base R functions, drastically decreasing the computational time required to obtain
#' predictive draws.
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
#' @param shape Hyperparameter alpha for gamma prior on the precision of the ANOVA model, tau.
#' @param rate Hyperparameter beta for gamma prior on the precision of the ANOVA model, tau.
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
#' @param seed Optional selection which will create a reproducible result from the function.
#' @param verbose Allows suppression of sampler progress printing in console.
#' @return Returns a list with three elements:
#'   \describe{
#'     \item{\code{pp}}{The predicted probability of the test ending in a successful evaluation of the question of interest}
#'     \item{\code{posterior}}{The full dataframe of non-conditional posterior draws}
#'     \item{\code{indicator}}{The vector of test success results for each posterior draw}
#'   }
#'   Printing the result object will display the predicted probability result.
#' @importFrom stats aggregate rgamma rmultinom rnorm
#' @export
gibbs.sampler.predictive.rjags <- function(X, Y, n.seen, beta.mean, beta.precision, shape, rate,
                                           n.sim, y.burnin, b.sim, b.burnin,
                                           phi.0, theta.t, prob, factor.no.2way = NA, colnames.pick = NA,
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
  if(! is.numeric(n.seen)) { stop("n.seen must be a numeric type.") }
  if(! is.numeric(beta.mean)) { stop("beta.mean must be a numeric type.") }
  if(! is.numeric(shape)) { stop("shape must be a numeric type.") }
  if(! is.numeric(rate)) { stop("rate must be a numeric type.") }
  if(! is.numeric(n.sim)) { stop("n.sim must be a numeric type.") }
  if(! is.numeric(y.burnin)) { stop("y.burnin must be a numeric type.") }
  if(! is.numeric(b.sim)) { stop("b.sim must be a numeric type.") }
  if(! is.numeric(b.burnin)) { stop("b.burnin must be a numeric type.") }
  if(! is.numeric(phi.0)) { stop("phi.0 must be a numeric type.") }
  if(! is.numeric(theta.t)) { stop("theta.t must be a numeric type.") }
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

  # ERROR: Theta.t should be 0-1
  if(! (theta.t >= 0 & theta.t <= 1)) { stop("theta.t is not a probability in [0,1]") }

  # ERROR: n.seen should be less than size of design matrix X
  if(n.seen > nrow(X)) { stop("n.seen value must be less than or equal to the number of rows available in the design matrix X") }

  # ERROR: n.seen should be less than number of test responses Y
  if(n.seen > length(Y)) { stop("n.seen value must be less than or equal to the number of rows available in the response vector Y") }

  # ERROR: gamma prior on tau parameters must be positive
  if(! (shape > 0)) { stop("Gamma parameters must be greater than 0 (shape is not)") }
  if(! (rate > 0)) { stop("Gamma parameters must be greater than 0 (rate is not)") }

  # ERROR: Posterior not predictive probability
  if(n.seen == nrow(X)) {
    stop(
      "Since n.seen is the same as the number of rows of the design matrix,
      you should be using the gibbs.sampler.posterior function."
    )
  }

  ##### Warnings to user when function will run alright but something is strange

  # WARNING: Burn-in has more draws than actual simulation
  if(y.burnin > n.sim) {
    warning(
      "Burn-in draws (y.burnin) higher than non-conditional posterior draws (n.sim).",
      immediate. = TRUE
    )
  }
  if(b.burnin > b.sim) {
    warning(
      "Burn-in draws (b.burnin) higher than non-conditional posterior draws (b.sim).",
      immediate. = TRUE
    )
  }

  # WARNING: Unobtainable certainty thresholds for successful passing of test
  if(theta.t == 0 | theta.t == 1) {
    warning(
      "A theta.t equal to 0 or 1 are not obtainable, so the conclusion with these values will always be
       that the test was or was not successfully passed.",
      immediate. = TRUE
    )
  }

  ##### The user provides b.sim and n.sim, expecting those to be the number of conditional posterior draws or
  # (non-conditional) posterior draws returned; however code uses b.sim and n.sim as the total number sampled
  # in the inner and outer loop. So here we re-define b.sim to be b.sim provided by the user plus b.burnin and
  # we re-define n.sim to be n.sim provided by the user plus y.burnin.
  ##### Note, n.sim needs an additional "+ 1" because the inner loop defines the indicator function for the
  # previous outer draw. So to get N observations from the outer loop, you need N + 1 draws, so you can drop
  # the last line that will be NA because it didn't have an inner draw.
  n.sim <- n.sim + y.burnin + 1
  b.sim <- b.sim + b.burnin

  # If needed, make the length of Y match X for creating the full design matrix
  if(length(Y) < nrow(X)){ Y <- matrix(c( Y, rep(NA_real_, nrow(X) - length(Y)) ), ncol = 1) }

  ##### Create a matrix with the response in the first column, and the remaining columns as the design matrix
  full.design <- cbind(Y, X)

  ##### Get the number of model parameters -- not including tau
  num.param <- length(beta.mean)

  ##### Create "success.ind" to keep track of indicators of a test being successful at the end of the test or not
  success.ind <- matrix(NA_real_, nrow=(n.sim - 1), ncol=1)

  ##### Get column name to assign to results (if not provided, use colnames of design matrix)
  if(is.na(colnames.pick[1])) {colnames.pick <- c(colnames(X), "tau")}

  ##### Create an object "posterior.results" to keep track of posterior draws and return from function
  posterior.results <- matrix(NA_real_, nrow=n.sim, ncol=(num.param + 1), dimnames=list(NULL, colnames.pick))

  ##### Calculate initial value for tau as mean of gamma distribution
  tau.int <- shape / rate

  ##### Store the initial values for the model parameters in the matrix that was just created
  posterior.results[1, ] <- c(beta.mean, tau.int)

  ##### For predictive probability calculations, you will not have seen every observation yet;
  # this "NA"s out the rows after the last "n.seen" observation
  full.design[(n.seen + 1):nrow(full.design), 1] <- NA_real_




  ############################################################################################################
  #############################################   OUTSIDE LOOP   #############################################
  ############################################################################################################

  for(o in 2:n.sim) {

    # Progress printing
    if(verbose) {
      if(o - 1 <= y.burnin) {
        cat("\r", paste("Running Burn-in", o - 1, "of", y.burnin))
      } else {
        cat("\r", paste("Running Simulation", o - y.burnin - 1, "of", n.sim - y.burnin - 1))
      }
    }

    ##### Generate the remaining data yet to be seen (n.all - n.seen left to observe) and put them into the full design
    # (i.e., augment the full design with the generated observations)

    # Set seed for generating future observations for reproduceability
    if(is.na(seed) == FALSE){
      set.seed(seed + o - 2)
    }

    for(g in (n.seen + 1):nrow(full.design)) {
      ##### Get the true mean based on the parameter estimates.
      ##### Data is in the first column of the object "full design", so full.design[g, 2:ncol(full.design)%*%posterior.results[o-1, 1:num.param] is the X\beta matrix
      dist.mean <- full.design[g, 2:ncol(full.design)] %*% posterior.results[o-1, 1:num.param]
      ##### Now that we have the true mean based on the parameter estimates, use that to draw random draws from a normal distribution
      # centered at the mean; posterior.results[o-1, (num.param + 1)] is tau, so the std is sqrt(1/posterior.results[o-1, (num.param + 1)])
      full.design[g, 1] <- rnorm(1, dist.mean, sqrt(1/posterior.results[o-1, (num.param + 1)]))
    }

    # Reassign Y to include the newly imputed Y values that have not been observed.
    Y <- as.matrix(full.design[,1])


    ############################################################################################################
    #############################################   MISSION SETS   #############################################
    ############################################################################################################

    ##### For each inner loop, we need a mission set matrix that is num.param columns (1 for each parameter in the model)
    # and b.sim - b.burnin rows (1 for each conditional posterior draw that we will keep after burn in).
    ##### This grabs the mission for each column of the matrix, e.g. the first one grabs nothing but 1's because
    # the model parameter eta is the intercept essentially. The rest are generated from binomial or multinomial distributions,
    # as defined by the likelihood of encountering the parameter. Once combined, this essentially will create a new design matrix
    # of possible scenarios we could see in the field and evlaute one conditional posterior draw against one of the rows in the mission sets.

    num.missions <- b.sim-b.burnin
    num.maineffects <- length(prob[, 1]) - length(unique(prob[, 1]))

    missions <- matrix(NA, ncol=(1 + num.maineffects), nrow=num.missions)
    missions[, 1] <- 1
    for(f in 1:max(prob[, 1])) {
      sample.mission.set <- rmultinom(num.missions, 1, prob[prob[, 1] == f, 2])
      last.notna.col <- sum(!is.na(missions)) / num.missions
      ##### This takes the number of rows sampled (again, in case it isn't a two level factor),and puts it in
      # the next columns filled with NA. If there is a two level factor it will take the second row of the multinomial
      # vector and put it in one column. If it is a three-level factor it will take the second and third row of the multinomial vector
      # and put it in one column. So on and so forth.
      missions[, (last.notna.col + 1):(last.notna.col + nrow(sample.mission.set)-1)] <- t(sample.mission.set)[, 2:ncol(t(sample.mission.set))]
    }

    ##### Label each row by the factor number from the object prob, so that main effects that are not included in the
    # two-way interactions can be removed by their factor number value.
    colnames(missions) <- c(0, prob[duplicated(prob[, 1]), 1])

    ##### Remove the factors not needed for two-way interactions
    if(! is.na(factor.no.2way)) {
      sub.missions <- missions[, colnames(missions) != factor.no.2way] # This currently drops by name
    } else {
      sub.missions <- missions # I added this since factor.no.2way = NA meant it couldn't find a sub.missions object
    }

    ##### Remove the intercept / reference cell number
    sub.missions <- sub.missions[,-1]

    ##### Get the number of two-way interactions needed from the remaining values:
    # (n * (n-1))/2 gives you the total number of interactions that are possible, with no constraints
    # length(prob[duplicated(prob[, 1]), 1][duplicated(prob[duplicated(prob[, 1]), 1])]) can be broken down as follows:
    #     prob[duplicated(prob[, 1]), 1] tells you what values are duplicated (all factors will have at least one duplicate; the number of duplicates
    #         are the number of levels for each factor minus 1)
    #     prob[duplicated(prob[, 1]), 1][duplicated(prob[duplicated(prob[, 1]), 1])] tells you what values are duplicated of the duplicated; any factor that has only two
    #         levels will not be presented here, only factors that have more than two levels. For example, if factor 1 had four levels, the previous commented code would
    #         result in "1 1 1" while this commented code would result in "1 1".
    #     length(prob[duplicated(prob[, 1]), 1][duplicated(prob[duplicated(prob[, 1]), 1])]) then takes the length of the previous commented code (in our example 2);
    # This is subtracted off of the first piece becuase these are the number of columns that cannot have two way interactions with themselves. (i.e., factor levels are
    #    independent of themselves and will have no interactoin).

    num.levels.repeat <- sum( prob[duplicated(prob[, 1]), 1][duplicated(prob[duplicated(prob[, 1]), 1])] != factor.no.2way )
    num.int <- ((ncol(sub.missions)) * (ncol(sub.missions) - 1)) / 2 - num.levels.repeat

    ##### Create a matrix to store the two-way interactions in
    missions.int <- matrix(NA, ncol = num.int, nrow = num.missions)

    ##### Get two-way interactions by multiplying each main effect column (element-wise) by the main effect columns that follow it, one column at a time.
    # If the main effect column that follows at any point is the same as the main effect column name that is being used then the element-wise muliplication is skipped
    #   and the next column is considered. This avoids having two way interactions within a factor.
    #  If the column name is the same, the
    for(t in 1:(ncol(sub.missions) - 1)) {
      for(u in (t + 1):ncol(sub.missions)) {
        if( colnames(sub.missions)[t] != colnames(sub.missions)[u] ) {
          missions.int[, ((sum(! is.na(missions.int)) / num.missions) + 1)] <- sub.missions[, t] * sub.missions[, u]
        }
      }
    }

    ##### Create the full mission set
    mission.set <- cbind(missions, missions.int)
    mission.set <- t(mission.set)

    ##### Create an object "post" for saving individual conditional posterior draws for the inner loop
    post <- matrix(NA_real_, nrow = b.sim, ncol = (num.param + 2), dimnames = list(NULL, c(colnames.pick, "m")))

    ##### Store initial values in the first row for inner loop
    post[1,] <- c(posterior.results[o-1, 1:(num.param + 1)], NA_real_)


    ############################################################################################################################
    #############################################   INSIDE LOOP - RJAGS APPROACH   #############################################
    ############################################################################################################################

    parameters <- c( "beta", "tau" )
    data <- list( "X"=X, "n"=nrow(X), "Y"=Y, "beta.mean"=c(beta.mean), "beta.precision"=c(beta.precision), "m"=length(beta.mean), "shape"=shape, "rate"=rate)
    if(is.na(seed) == FALSE){
      inits <- list( ".RNG.name"= "base::Wichmann-Hill", ".RNG.seed"= (seed + o - 2), beta=c(beta.mean), tau=tau.int )
    } else{
      inits <- list( beta=c(beta.mean), tau=tau.int )
    }


    modelstring<-"model
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

    ##### Remove Burn in from conditional posterior draws (burnin function within rjags not used above, as it did not properly
    # demonstrate random draws from when used)
    post.wobi <- post.wbi[-(1:b.burnin), ]

    ##### Calculate mission mean using the conditional posterior draws and mission set matrix. To do this, first row of
    # post.wobi multiplied by first column of mission.set, second row of post.wobi multiplied by second column of
    # mission.set, etc (alternatively, you could see this as the diagonal of the posterior draws matrix multiplied by the
    # mission set matrix). Save the result in the last column of post.wobi.
    for(s in 1:nrow(post.wobi)) {
      post.wobi[s, (num.param + 2)] <- post.wobi[s, 1:num.param] %*% mission.set[, s]
    }

    ###### Calculate the conditional probability of the test being a success at the end of the test; this is, the number
    # of mission means that were above phi.0 divided by the total number of mission means to get Pr.
    cond.prob <- sum(post.wobi[, (num.param + 2)] > phi.0) / nrow(post.wobi)

    ###### If the conditional probability is greater than theta.t, assign a 1 for success; if not, assign a 0 for failure
    if(cond.prob >= theta.t){ success.ind[o-1, 1] <- 1 } else{ success.ind[o-1, 1] <- 0 }

    ##### Add draw results to final table
    posterior.results[o, 1:(num.param + 1)] <- post.wobi[1, 1:(num.param + 1)]

    if(verbose & o == n.sim) cat("\n Simulation complete \n")


  }
  ##### End of outer loop

  ##### Remove burn-in observations from the outer draws ("y.burnin") to return only the user defined n.sim requested values
  posterior.results.wobi<-posterior.results[-(1:y.burnin), ]
  success.ind <- success.ind[-(1:y.burnin), ]

  ##### Remove the last outer draw because it doesn't have a value for the indicator function--the information we care about.
  posterior.results.wobi<-posterior.results.wobi[-( nrow(posterior.results.wobi) ), ]

  ##### Convert the posterior draws matrix to a data.frame and indicator matrix of successes to a vector for return purposes
  posterior.results.wobi <- as.data.frame(posterior.results.wobi)
  success.ind <- as.vector(success.ind)

  ##### Calculate the predictive probability of the test being a success at the end
  pp <- sum(success.ind)/( length(success.ind) )

  ##### Return the predictive probability, the posterior draws from the outer loop, and the vector of indicators defining
  # whether a draw is a success (1) or not (0) from the function.
  # Set a new class so that when users print results, the console only shows pp

  output <- structure(
    list(
      pp = pp,
      posterior = posterior.results.wobi,
      indicator = success.ind
    ),
    class = "ContRespPP"
  )

  return(output)

}
