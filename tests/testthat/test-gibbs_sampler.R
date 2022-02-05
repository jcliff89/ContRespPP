library(ContRespPP)

### Set parameters
full.data <- ContRespPP::exData

X <- full.data[,c(2:14)]
Y <- as.matrix(full.data[,1], ncol = 1)

theta.t <- 0.8
phi.0 <- 400

beta.mean <- matrix(c(400, 50, 50, -25, -50, 100, 100, rep(0, 6)),ncol=1)
beta.precision <- matrix(c(1/10000, 1/10000, 1/10000, 1/2500, 1/2500, 1/10000, 1/10000, rep(1/10000,6)), ncol=1)
precision.a <- 0.0001
precision.b <- 0.0001

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

# Run predictive results to compare outputs in tests
set.seed(512)
predictive_results <- gibbs.sampler(
  X = X,
  Y = Y,
  n.seen = n.seen,
  beta.mean = beta.mean,
  beta.precision = beta.precision,
  precision.a = precision.a,
  precision.b = precision.b,
  n.sim = n.sim,
  y.burnin = y.burnin,
  b.sim = b.sim,
  b.burnin = b.burnin,
  phi.0 = phi.0,
  theta.t = theta.t,
  prob = prob,
  factor.no.2way = factor.no.2way,
  colnames.pick = colnames.pick
)

# Run posterior results to compare outputs
set.seed(512)
posterior_results <- gibbs.sampler.posterior(
  X = X[1:75,],
  Y = Y[1:75],
  n.seen = n.seen,
  beta.mean = beta.mean,
  beta.precision = beta.precision,
  precision.a = precision.a,
  precision.b = precision.b,
  b.sim = b.sim,
  b.burnin = b.burnin,
  prob = prob,
  factor.no.2way = factor.no.2way,
  colnames.pick = colnames.pick
)


# Runs tests on resulting objects elements versus the internal results objects
testthat::test_that("predictive pp works", {
  expect_equal(predictive_results$pp, 0.7)
})
test_that("predictive posterior works", {
  expect_equal(
    predictive_results$posterior[10,],
    data.frame(
      eta = 349.11155514645,
      alpha = 39.69934317663,
      beta = 26.35416334985,
      omega2 = -8.96002391443,
      omega3 = 7.04971687008,
      theta = 100.87461920453,
      gamma = 42.17204449009,
      alphabeta = -10.84151893116,
      alphatheta = 37.85817304383,
      alphagamma = 31.35655526784,
      betatheta = 29.272223,
      betagamma = 36.625298,
      thetagamma = 9.5689962,
      tau = 0.00040760768,
      row.names = as.integer(10)
    )
  )
})
test_that("predictive indicator works", {
  expect_equal(predictive_results$indicator, c(1, 1, 1, 0, 1, 1, 1, 0, 0, 1))
})
test_that("predictive posterior works", {
  expect_equal(
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
      betatheta = 71.0335788,
      betagamma = 37.5565594,
      thetagamma = 25.1356535,
      tau = 0.0004231,
      m = 345.4537665
    )
  )
})


# Clean up environment after tests
rm(full.data, X, Y, n.seen, beta.mean, beta.precision, precision.a, precision.b,
   n.sim, y.burnin, b.sim, b.burnin, phi.0, theta.t, prob, factor.no.2way, colnames.pick,
   posterior_results, predictive_results)
