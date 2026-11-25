# rpoilog and rbipoilog. Seeds are fixed; statistical tests use a lenient level
# so that they are not flaky across platforms.

test_that("rpoilog returns the requested number of counts", {
  set.seed(10)
  expect_length(rpoilog(500, 0, 1, keep0 = TRUE), 500)
  x <- rpoilog(500, 0, 1, condS = TRUE)
  expect_length(x, 500)
  expect_true(all(x > 0))
  y <- rpoilog(500, 0, 1)
  expect_true(length(y) <= 500 && all(y > 0))
})

test_that("rpoilog returns non-negative whole numbers", {
  set.seed(11)
  x <- rpoilog(1000, 1, 2, keep0 = TRUE)
  expect_true(all(x >= 0 & x == round(x)))
})

test_that("rpoilog draws from dpoilog (chi-square test)", {
  set.seed(12)
  mu <- 0.5; sig <- 1
  x <- rpoilog(20000, mu, sig, keep0 = TRUE)
  k <- 0:14
  p <- c(dpoilog(k, mu, sig), 1 - sum(dpoilog(k, mu, sig)))
  obs <- c(tabulate(x + 1, nbins = 15), sum(x > 14))
  expect_gt(suppressWarnings(chisq.test(obs, p = p)$p.value), 1e-4)
})

test_that("nu scales the mean like a shift of mu by log(nu)", {
  set.seed(13)
  x <- rpoilog(20000, 0, 0.5, nu = 3, keep0 = TRUE)
  m <- exp(log(3) + 0.5^2 / 2)
  expect_equal(mean(x), m, tolerance = 0.03)
})

test_that("rbipoilog returns the requested shape and drops double zeros", {
  set.seed(14)
  xy <- rbipoilog(400, 0, 0, 1, 1, 0.5, keep0 = TRUE)
  expect_equal(dim(xy), c(400L, 2L))
  xy <- rbipoilog(400, 0, 0, 1, 1, 0.5, condS = TRUE)
  expect_equal(nrow(xy), 400L)
  expect_true(all(rowSums(xy) > 0))
  expect_true(all(xy >= 0 & xy == round(xy)))
})

test_that("rbipoilog draws from dbipoilog (chi-square test on a grid)", {
  set.seed(15)
  xy <- rbipoilog(20000, 0, 0.3, 0.7, 0.6, 0.4, keep0 = TRUE)
  g <- expand.grid(n1 = 0:3, n2 = 0:3)
  p <- dbipoilog(g$n1, g$n2, 0, 0.3, 0.7, 0.6, 0.4)
  cell <- ifelse(xy[, 1] <= 3 & xy[, 2] <= 3, xy[, 1] + 4 * xy[, 2] + 1, 17)
  obs <- tabulate(cell, nbins = 17)
  expect_gt(suppressWarnings(chisq.test(obs, p = c(p, 1 - sum(p)))$p.value), 1e-4)
})

test_that("rbipoilog gives positively correlated counts for rho > 0", {
  set.seed(16)
  xy <- rbipoilog(5000, 1, 1, 1, 1, 0.8, keep0 = TRUE)
  expect_gt(cor(log1p(xy[, 1]), log1p(xy[, 2])), 0.3)
})

test_that("random generators reject invalid input", {
  expect_error(rpoilog(0, 0, 1), "not positive")
  expect_error(rpoilog(10.5, 0, 1), "not an integer")
  expect_error(rpoilog(10, 0, -1), "sig is not positive")
  expect_error(rpoilog(10, 0, 1, nu = -1), "nu is not positive")
  expect_error(rbipoilog(10, 0, 0, -1, 1, 0.5), "sig1 is not positive")
  expect_error(rbipoilog(10, 0, 0, 1, -1, 0.5), "sig2 is not positive")
  expect_error(rbipoilog(10, 0, 0, 1, 1, 0.5, nu1 = -1), "nu1 is not positive")
  expect_error(rbipoilog(10, 0, 0, 1, 1, 0.5, nu2 = -1), "nu2 is not positive")
  expect_error(rpoilog(10, numeric(0), 1), "single number")
  expect_error(rpoilog(10, 0, 1, condS = NA), "TRUE or FALSE")
  expect_error(rbipoilog(10, 0, 0, 1, 1, 2), "rho")
  expect_error(rbipoilog(10, 0, 0, 1, 1, 0.5, keep0 = "yes"), "TRUE or FALSE")
})

test_that("an intensity too large for a double stops instead of giving NA counts", {
  # exp(710) is Inf, for which rpois gives NA (returned as a count by 0.4.2.1)
  expect_error(rpoilog(3, 710, 0), "exceeds the largest double")
  expect_error(rpoilog(3, 710, 0, condS = TRUE), "exceeds the largest double")
  expect_error(rbipoilog(3, 0, 710, 1, 0, 0.5), "exceeds the largest double")
  set.seed(1)
  expect_error(rpoilog(100, 0, 1000), "exceeds the largest double")
})

test_that("impossible conditional simulations stop instead of looping", {
  expect_error(rpoilog(10, 0, 1, nu = 0, condS = TRUE), "impossible")
  expect_error(rbipoilog(10, 0, 0, 1, 1, 0.5, nu1 = 0, nu2 = 0, condS = TRUE), "impossible")
  # with keep0 = TRUE zeros count, so S zeros are the answer, as in poilog 0.4.2.1
  expect_identical(rpoilog(10, 0, 1, nu = 0, condS = TRUE, keep0 = TRUE), rep(0, 10))
  expect_identical(rbipoilog(10, 0, 0, 1, 1, 0.5, nu1 = 0, nu2 = 0, condS = TRUE, keep0 = TRUE),
                   matrix(0L, 10, 2, dimnames = list(NULL, c('x', 'y'))))
  # P(N > 0) is about 1e-26 here
  expect_error(rpoilog(1, -60, 0.5, condS = TRUE), "too small")
})

test_that("seeded draws are identical to those of poilog 0.4.2.1", {
  # drawn with poilog 0.4.2.1 by tests/reference/draws-0.4.2.1.R
  # (compared as doubles: some generators return integer vectors)
  set.seed(101)
  expect_identical(as.numeric(rpoilog(12, 1, 1.5, keep0 = TRUE)),
                   c(4, 8, 0, 3, 4, 17, 6, 1, 7, 1, 8, 0))
  set.seed(102)
  expect_identical(as.numeric(rpoilog(30, 0, 2)),
                   c(3, 1, 56, 14, 16, 6, 2, 65, 19, 1, 35, 2, 11, 1, 1, 1, 13, 53))
  set.seed(103)
  expect_identical(as.numeric(rpoilog(10, -1, 2, nu = 0.5, condS = TRUE)),
                   c(1, 2, 23, 1, 1, 1, 1, 1, 1, 4))
  set.seed(104)
  xy <- rbipoilog(8, 0, 0.5, 1, 1.2, 0.6, keep0 = TRUE)
  expect_identical(unname(xy) + 0,
                   cbind(c(2, 0, 0, 1, 1, 4, 0, 2), c(1, 0, 6, 0, 5, 0, 0, 4)))
  set.seed(105)
  xy <- rbipoilog(10, -1, -0.5, 2, 1.5, -0.3, nu1 = 0.7, condS = TRUE)
  expect_identical(unname(xy) + 0,
                   cbind(c(0, 0, 4, 1, 0, 0, 0, 1, 1, 0), c(3, 14, 0, 4, 1, 4, 1, 0, 1, 1)))
})

test_that("rbipoilog with S = 1 and condS = TRUE returns a one-row matrix", {
  set.seed(17)
  xy <- rbipoilog(1, 0, 0, 1, 1, 0.5, condS = TRUE)
  expect_equal(dim(xy), c(1L, 2L))
  expect_gt(sum(xy), 0)
})
