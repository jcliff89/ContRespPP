#' Continuous Response Posterior Probability.
#'
#' \code{gibbs.sampler.posterior} obtains the ContRespPP method for obtaining the posterior (i.e., all
#' observations have been seen, and this reverts to a traditional Bayesian analysis) using base R
#' functions, drastically increasing the computational time required to obtain predictive draws.
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
gibbs.sampler.posterior <- function(X, Y, beta.mean, beta.precision, shape, rate,
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
  if(! is.numeric(beta.precision)) { stop("beta.precision must be a numeric type.") }
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

  ##### Get the number of observations for each model parameter
  n.obs.per.modparam <- as.matrix(apply(X, 2, sum), ncol=1)

  ##### Get the total number of observations (which is also the number of observations for the reference parameter)
  n.all <- n.obs.per.modparam[1,1]

  ##### Create a matrix that stores the prior mean, prior precision, number of observations on each model parameter,
  # and the column number it is in the design matrix
  prior.tab <- cbind(beta.precision, beta.mean, n.obs.per.modparam, matrix(1:length(n.obs.per.modparam)))

  ##### Reorder rows so the intercept (either a reference cell or a grand mean) is the last thing in the matrix--
  # because it is the last thing to update in the Gibbs sampler.
  prior.tab <- rbind(prior.tab[2:nrow(prior.tab), ], prior.tab[1,])
  colnames(prior.tab) <- c("beta.precision", "beta.mean", "num.obs", "design.matrix.colnum")

  ##### Get column name to assign to results (if not provided, use colnames of design matrix)
  if(is.na(colnames.pick[1])) {colnames.pick <- c(colnames(X), "tau")}

  ##### Set seed for function if provided
  if(is.na(seed) == FALSE){
    set.seed(seed)
  }

  ############################################################################################################
  #############################################   SUBSET DATA    #############################################
  ############################################################################################################

  ##### Subset data based on parameters as discussed in the function section of this code.

  data.list <- vector("list", num.param)
  ##### Intercept or Reference Cell Parameter
  data.list[[num.param]] <- parameter.data(2, Y, X, full.design)
  for (d in 2:num.param) {
    data.list[[(d - 1)]] <- parameter.data(d + 1, Y, X, full.design)
  }


  ############################################################################################################
  #############################################   MISSION SETS   #############################################
  ############################################################################################################

  ##### We need a mission set matrix that is num.param columns (1 for each parameter in the model)
  # and b.sim - b.burnin rows (1 for each posterior draw that we will keep after burn in).
  ##### This grabs the mission for each column of the matrix, e.g. the first one grabs nothing but 1's because
  # the model parameter eta is the intercept essentially. The rest are generated from binomial or multinomial distributions,
  # as defined by the likelihood of encountering the parameter. Once combined, this essentially will create a new design matrix
  # of possible scenarios we could see in the field and evlaute one conditional posterior draw against one of the rows in the mission sets.

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


  ##### Create an object "post" to keep track of posterior draws and return from function
  post <-
    matrix(
      NA_real_,
      nrow = b.sim,
      ncol = (num.param + 2),
      dimnames = list(NULL, c(colnames.pick, "m"))
    )

  #### Calculate initial value for tau as mean of gamma distribution
  tau.int <- shape / rate

  ##### Store initial values in the first row
  post[1,] <- c(beta.mean, tau.int, NA_real_)


  ###########################################################################################################
  #############################################   INSIDE LOOP   #############################################
  ###########################################################################################################

  for (k in 2:b.sim) {
    # Progress printing
    if(verbose) {
      if (k <= b.burnin) {
        cat("\r", paste("Running Burn-in", k, "of", b.burnin))
      } else {
        cat("\r",
            paste("Running Simulation", k - b.burnin, "of", b.sim - b.burnin))
      }
    }

    ##### We start with the entire previous iteration. As model parameters are updated, they are not only saved in post,
    # but they are saved in param.it so they can be used in the smapler. So, for exampe, since the first full conditional
    # will use the entire row of the previous iteration. However, the next model parameter will have only one observation
    # (the first, that is calculated before second), so it will take the current first parameter after updating the first
    # model parameter, and then everything else will be the previous model parameter values. This is to get the right values
    # multipled for each full conditional.

    param.it <- matrix(post[(k - 1), 1:num.param], ncol = 1)


    ##### Conditional Gibbs sampler for model parameters (i.e., not tau)

    for (p in 1:(num.param - 1)) {
      ##### Sample parameter | data and other parameters (k or k-1)

      ##### First, calculate variance for parameter | data and other parameters, called "param.var"
      ##### post[k-1, (num.param + 1)] is the precision from the previous conditional posterior draw; prior.tab[p,3]
      # is the number of observations seen for the p-th parameter; prior.tab[p,1] is the prior precision for the p-th parameter
      param.var <-
        1 / (post[(k - 1), (num.param + 1)] * prior.tab[p, 3] + prior.tab[p, 1])

      ##### Next, calculate the mean for parameter | data and other parameters
      ##### The mean has the form (b/a)+(c/a).
      ##### This calculates the "a" from the previous comment; post[k-1, (num.param + 1)] is the precision from the previous
      # conditional posterior draw; prior.tab[p,3] is the number of observations seen for the p-th parameter; prior.tab[p1,1]
      # is the prior precision for the p-th parameter
      bottom.frac <-
        prior.tab[p, 3] + (prior.tab[p, 1] / post[(k - 1), (num.param + 1)])

      ##### This removes the value of the parameter list that is under consideration now (e.g., the p-th model parameter is
      # column (p+1), so it is removed from the list)
      param.it.mod <- as.matrix(param.it[-(prior.tab[p, 4]), 1])
      ##### This takes the data list, grabs the p-th one, and subtracts off the expected values
      # (based on the param.it.mod multiplication term) so that all you have left is the piece just from the p-th model parameter
      sum.diff <-
        sum(data.list[[p]][, 1] - data.list[[p]][, 2:ncol(data.list[[1]])] %*% param.it.mod)
      ##### Finally, we actually calculate the parameter estimate for the p-th model parameter
      param.hat <-
        (1 / bottom.frac) * sum.diff + ((prior.tab[p, 1] / post[(k - 1), (num.param + 1)]) / bottom.frac) * prior.tab[p, 2]

      ##### We now sample parameter | data and other parameters, based on the mean and variance we just calculated.
      ##### place the p-th estimate in the p+1 spot because the intercept / reference cell is calculated last.
      post[k, (p + 1)] <- rnorm(1, param.hat, sqrt(param.var))
      param.it[(p + 1), 1] <- post[k, (p + 1)]

    }


    ##### INTERCEPT / REFERENCE CELL #####

    #Sample parameter | data and other parameters (k or k-1)
    #calculate variance for parameter | data and other parameters
    param.var <-
      1 / (post[(k - 1), (num.param + 1)] * prior.tab[num.param, 3] + prior.tab[num.param, 1])

    #calculate mean for parameter | data and other parameters
    bottom.frac <-
      prior.tab[num.param, 3] + (prior.tab[num.param, 1] / post[(k - 1), (num.param + 1)])
    param.it.ref <-
      as.matrix(param.it[-(prior.tab[num.param, 4]), 1])

    sum.diff <-
      sum(data.list[[num.param]][, 1] - data.list[[num.param]][, 2:ncol(data.list[[num.param]])] %*% param.it.ref)
    param.hat <-
      (1 / bottom.frac) * sum.diff + ((prior.tab[num.param, 1] / post[(k - 1), (num.param + 1)]) / bottom.frac) * prior.tab[num.param, 2]

    #Sample parameter | data and other parameters
    post[k, 1] <- rnorm(1, param.hat, sqrt(param.var))
    param.it[1, 1] <- post[k, 1]


    ##### PRECISION #####

    #Sample parameter | data and other parameters (k or k-1)
    #a for parameter | data and other parameters
    a.val <- shape + n.all / 2

    #b for parameter | data and other parameters
    tau.sum <-
      sum((full.design[, 1] - full.design[, 2:ncol(full.design)] %*% param.it) ^ 2)
    b.val <- rate + 0.5 * tau.sum

    #Sample parameter | data and other parameters
    post[k, (num.param + 1)] <- rgamma(1, a.val, b.val)

  }

  ##### End of the Inner Loop

  ####################################################################################################################################
  ####################################################################################################################################
  #############################################								           #############################################
  #############################################   SAVE VARIOUS OUTPUTS / RETURN VALUES   #############################################
  #############################################					       ` 			   #############################################
  ####################################################################################################################################
  ####################################################################################################################################

  ##### Remove Burn in from posterior draws
  post.wobi <- post[-(1:b.burnin),]

  ##### Calculate mission mean using the conditional posterior draws and mission set matrix. To do this, first row of
  # post.wobi multiplied by first column of mission.set, second row of post.wobi multiplied by second column of
  # mission.set, etc (alternatively, you could see this as the diagonal of the posterior draws matrix multiplied by the
  # mission set matrix). Save the result in the last column of post.wobi.
  for (s in 1:nrow(post.wobi)) {
    post.wobi[s, (num.param + 2)] <-
      post.wobi[s, 1:num.param] %*% mission.set[, s]
  }

  ###### Calculate the conditional probability of the test being a success at the end of the test; this is, the number
  # of mission means that were above phi.0 divided by the total number of mission means to get Pr.
  post.prob <- sum(post.wobi[, (num.param + 2)] > phi.0) / nrow(post.wobi)

  if (verbose & k == b.sim)
    cat("\n Simulation complete \n")

  output <- structure(
    list(
      probability = post.prob,
      posterior = post.wobi
    ),
    class = "ContRespPP"
  )

  return(output)

}
