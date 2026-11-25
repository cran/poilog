# General properties of the univariate and bivariate pmf, independent of any
# reference values.

test_that("dpoilog values are probabilities", {
  for (p in list(c(0, 1), c(-3, 2), c(4, 0.3), c(1, 6))) {
    q <- dpoilog(0:200, p[1], p[2])
    expect_true(all(q >= 0 & q <= 1))
    expect_true(sum(q) <= 1 + 1e-12)
    expect_true(all(dpoilog(0:200, p[1], p[2], log = TRUE) <= 0))
  }
})

test_that("dpoilog has the Poisson-lognormal mean and variance", {
  mu <- 0.5; sig <- 0.8; n <- 0:3000
  q <- dpoilog(n, mu, sig)
  m <- exp(mu + sig^2 / 2)
  expect_equal(sum(q), 1, tolerance = 1e-10)
  expect_equal(sum(n * q), m, tolerance = 1e-8)
  expect_equal(sum((n - m)^2 * q), m + m^2 * (exp(sig^2) - 1), tolerance = 1e-8)
})

test_that("dpoilog tends to the Poisson pmf as sig -> 0", {
  n <- 0:30
  expect_equal(dpoilog(n, 1.5, 1e-6), dpois(n, exp(1.5)), tolerance = 1e-5)
  expect_equal(dpoilog(n, 1.5, 1e-6, log = TRUE), dpois(n, exp(1.5), log = TRUE), tolerance = 1e-5)
})

test_that("dpoilog is vectorised over n and handles repeated values", {
  n <- c(7, 0, 3, 3, 120, 0, 1)
  expect_equal(dpoilog(n, 0.3, 1.2), vapply(n, dpoilog, 0, mu = 0.3, sig = 1.2))
  expect_equal(dpoilog(as.integer(n), 0.3, 1.2), dpoilog(n, 0.3, 1.2))
  expect_length(dpoilog(integer(0), 0, 1), 0)
})

test_that("dpoilog is finite on the log scale over a wide parameter range", {
  grid <- expand.grid(n = c(0, 1, 10, 1e4, 1e8), mu = c(-30, 0, 12), sig = c(1e-4, 1, 20))
  lq <- mapply(dpoilog, grid$n, grid$mu, grid$sig, MoreArgs = list(log = TRUE))
  expect_true(all(is.finite(lq) & lq <= 0))
})

test_that("dbipoilog values are probabilities and sum to 1", {
  # sig chosen so that the mass beyond 40 in either sample is below 1e-12
  g <- expand.grid(n1 = 0:40, n2 = 0:40)
  q <- dbipoilog(g$n1, g$n2, 0, 0.3, 0.45, 0.4, 0.4)
  expect_true(all(q >= 0 & q <= 1))
  expect_equal(sum(q), 1, tolerance = 1e-9)
})

test_that("dbipoilog marginalises to dpoilog", {
  # the mass of n2 beyond 200 is below 1e-12 here
  n2 <- 0:200
  for (n1 in c(0, 2, 9)) {
    m <- sum(dbipoilog(rep(n1, length(n2)), n2, 0.5, 0, 0.8, 0.6, -0.6))
    expect_equal(m, dpoilog(n1, 0.5, 0.8), tolerance = 1e-9)
  }
})

test_that("dbipoilog is symmetric under swapping the two samples", {
  a <- dbipoilog(c(0, 3, 12), c(5, 1, 0), 1, -0.5, 1.2, 0.7, 0.3, log = TRUE)
  b <- dbipoilog(c(5, 1, 0), c(0, 3, 12), -0.5, 1, 0.7, 1.2, 0.3, log = TRUE)
  expect_equal(a, b, tolerance = 1e-8)
})

test_that("dbipoilog: positive rho makes equal counts more likely", {
  expect_gt(dbipoilog(5, 5, 1, 1, 1, 1, 0.8), dbipoilog(5, 5, 1, 1, 1, 1, -0.8))
})

test_that("dbipoilog log = TRUE matches log() and is vectorised", {
  n1 <- c(0, 1, 4, 30); n2 <- c(2, 0, 4, 1)
  lq <- dbipoilog(n1, n2, 0, 1, 1, 1.5, 0.2, log = TRUE)
  expect_equal(lq, log(dbipoilog(n1, n2, 0, 1, 1, 1.5, 0.2)), tolerance = 1e-12)
  expect_equal(lq, mapply(dbipoilog, n1, n2, MoreArgs = list(mu1 = 0, mu2 = 1, sig1 = 1,
                          sig2 = 1.5, rho = 0.2, log = TRUE)))
})
