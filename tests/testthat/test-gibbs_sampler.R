library(ContRespPP)

### Set parameters
full.data <- ContRespPP::exData

X <- full.data[,c(2:14)]
Y <- as.matrix(full.data[,1], ncol = 1)

theta.t <- 0.8
phi.0 <- 400

beta.mean <- matrix(c(400, 50, 50, -25, -50, 100, 100, rep(0, 6)),ncol=1)
beta.precision <- matrix(c(1/10000, 1/10000, 1/10000, 1/2500, 1/2500, 1/10000, 1/10000, rep(1/10000,6)), ncol=1)
shape <- 0.0001
rate <- 0.0001

prob<-matrix(c(1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 1/2, 1/2, 4/9, 5/9, 1/3, 1/3, 1/3, 1/2, 1/2, 1/2, 1/2), ncol=2, dimnames=list(NULL, c("factor", "probability")))

factor.no.2way<-c(3)

colnames.pick<-c("eta", "alpha", "beta", "omega2", "omega3", "theta", "gamma", "alphabeta", "alphatheta", "alphagamma", "betatheta", "betagamma", "thetagamma", "tau")
colnames(X)<-colnames.pick[1:13]

n.seen<-75

b.sim<-10
b.burnin<-10

n.sim<-10
y.burnin<-10

factor.no.2way<-c(3)

# Run predictive results (base R)
predictive_results <- gibbs.sampler.predictive(
  X = X,
  Y = Y,
  n.seen = n.seen,
  beta.mean = beta.mean,
  beta.precision = beta.precision,
  shape = shape,
  rate = rate,
  n.sim = n.sim,
  y.burnin = y.burnin,
  b.sim = b.sim,
  b.burnin = b.burnin,
  phi.0 = phi.0,
  theta.t = theta.t,
  prob = prob,
  factor.no.2way = factor.no.2way,
  colnames.pick = colnames.pick,
  seed = 512
)

# Run posterior results (base R)
posterior_results <- gibbs.sampler.posterior(
  X = X[1:75,],
  Y = Y[1:75],
  beta.mean = beta.mean,
  beta.precision = beta.precision,
  shape = shape,
  rate = rate,
  phi.0 = phi.0,
  b.sim = b.sim,
  b.burnin = b.burnin,
  prob = prob,
  factor.no.2way = factor.no.2way,
  colnames.pick = colnames.pick,
  seed = 512
)

# Run predictive results (rjags)
predictive_results_jags <- gibbs.sampler.predictive.rjags(
  X = X,
  Y = Y,
  n.seen = n.seen,
  beta.mean = beta.mean,
  beta.precision = beta.precision,
  shape = shape,
  rate = rate,
  n.sim = n.sim,
  y.burnin = y.burnin,
  b.sim = b.sim,
  b.burnin = b.burnin,
  phi.0 = phi.0,
  theta.t = theta.t,
  prob = prob,
  factor.no.2way = factor.no.2way,
  colnames.pick = colnames.pick,
  seed = 512
)

# Run posterior results (rjags)
posterior_results_rjags <- gibbs.sampler.posterior.rjags(
  X = X[1:75,],
  Y = Y[1:75],
  beta.mean = beta.mean,
  beta.precision = beta.precision,
  shape = shape,
  rate = rate,
  phi.0 = phi.0,
  b.sim = b.sim,
  b.burnin = b.burnin,
  prob = prob,
  factor.no.2way = factor.no.2way,
  colnames.pick = colnames.pick,
  seed = 512
)


### Runs tests on resulting objects elements

# Predictive (base R)
testthat::test_that("base R predictive pp works", {
  testthat::expect_equal(predictive_results$pp, 0.7)
})
testthat::test_that("base R predictive posterior works", {
  testthat::expect_equal(
    predictive_results$posterior[10,],
    data.frame(
      eta = 347.54094,
      alpha = 49.919477,
      beta = 24.046309,
      omega2 = -12.109324,
      omega3 = 6.0590025,
      theta = 105.425247,
      gamma = 42.983682,
      alphabeta = -15.8259474,
      alphatheta = 44.147483,
      alphagamma = 17.2665694,
      betatheta = 37.598518,
      betagamma = 44.983309,
      thetagamma = 7.0175283,
      tau = 0.00041002082,
      row.names = as.integer(10)
    )
  )
})
testthat::test_that("base R predictive indicator works", {
  testthat::expect_equal(predictive_results$indicator, c(1, 1, 1, 0, 1, 1, 1, 0, 0, 1))
})

# Posterior (base R)
testthat::test_that("base R posterior probability works", {
  expect_equal(posterior_results$probability, 0.8)
})
testthat::test_that("base R predictive posterior works", {
  testthat::expect_equal(
    posterior_results$posterior[10,],
    c(
      eta = 372.0181307702488311,
      alpha = 72.4334729407908640,
      beta = 5.6437250088257507,
      omega2 = -18.9518715921377172,
      omega3 = -13.0444823033842709,
      theta = 62.9324083169688535,
      gamma = 25.3786141418328413,
      alphabeta = -25.0701648123980192,
      alphatheta = 48.6332235669522959,
      alphagamma = 28.6441252117864664,
      betatheta = 71.03357877,
      betagamma = 37.55655942,
      thetagamma = 25.13565349,
      tau = 0.00042312,
      m = 372.01813077
    )
  )
})

## Only run rjags test on Mac OS since JAGS RNG has differing results for other OS

if(Sys.info()['sysname'] == "Darwin") {
  # Predictive (rjags)
  testthat::test_that("rjags predictive pp works", {
    testthat::expect_equal(predictive_results_jags$pp, 0.7)
  })
  testthat::test_that("rjags predictive posterior works", {
    testthat::expect_equal(
      predictive_results_jags$posterior[10,],
      data.frame(
        eta = 351.12446,
        alpha = 62.189768,
        beta = 41.251491,
        omega2 = -11.5724895,
        omega3 = -17.2152391,
        theta = 83.784724,
        gamma = 31.757688,
        alphabeta = -34.809998,
        alphatheta = 42.67921,
        alphagamma = 48.548013,
        betatheta = 57.428093,
        betagamma = 16.6852106,
        thetagamma = 12.135702,
        tau = 0.00042180727,
        row.names = as.integer(10)
      )
    )
  })
  testthat::test_that("rjags predictive indicator works", {
    testthat::expect_equal(predictive_results_jags$indicator, c(1, 0, 1, 1, 1, 1, 1, 0, 0, 1))
  })

  # Posterior (rjags)
  testthat::test_that("rjags posterior probability works", {
    testthat::expect_equal(posterior_results_rjags$probability, 0.7)
  })
  testthat::test_that("rjags posterior works", {
    testthat::expect_equal(
      posterior_results_rjags$posterior[10,],
      c(
        eta = 376.06796183672054,
        alpha = 41.443110830794026,
        beta = -32.467186112679855,
        omega2 = 25.40386013084774,
        omega3 = -13.109862581334514,
        theta = 56.598592757067088,
        gamma = 28.991130508340987,
        alphabeta = 16.138461933670587,
        alphatheta = 46.697133566645284,
        alphagamma = 26.833356466019538,
        betatheta = 115.942451444620644,
        betagamma = 28.209858257950600,
        thetagamma = 14.022252094204337,
        tau = 0.00028087935,
        m = 376.06796183672
      )
    )
  })
}


# Clean up environment after tests
rm(full.data, X, Y, n.seen, beta.mean, beta.precision, rate, shape,
   n.sim, y.burnin, b.sim, b.burnin, phi.0, theta.t, prob, factor.no.2way, colnames.pick,
   posterior_results, predictive_results, posterior_results_rjags, predictive_results_jags)
