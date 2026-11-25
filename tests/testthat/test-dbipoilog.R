# Reference values: log q(n1, n2; mu1, mu2, sig1, sig2, rho) and
# log P(N1 + N2 > 0), computed by tests/reference/make-reference.R with an
# implementation that does not use poilog (tests/reference/ref.R).

ref_b <- data.frame(
  n1 = c(3, 0, 0, 0, 500, 4),
  n2 = c(5, 2, 0, 500, 0, 4),
  mu1 = c(1, 0.5, -2, 3, 3, 1),
  mu2 = c(1, 1, -2, 3, 3, 1),
  sig1 = c(1, 1.5, 5, 2, 2, 1),
  sig2 = c(1, 1, 5, 2, 2, 1),
  rho = c(0.5, -0.5, 0.5, 0.95, 0.95, 1),
  lq = c(-4.86111548781883, -3.32697207601703, -0.810534175298127,
         -37.8235590048638, -37.8235590048636, -4.40523370997353)
)

test_that("dbipoilog agrees with the reference values", {
  for (i in seq_len(nrow(ref_b))) {
    r <- ref_b[i, ]
    expect_equal(dbipoilog(r$n1, r$n2, r$mu1, r$mu2, r$sig1, r$sig2, r$rho, log = TRUE),
                 r$lq, tolerance = 1e-9,
                 label = sprintf("log dbipoilog(%g, %g, ..., rho = %g)", r$n1, r$n2, r$rho))
  }
})

test_that("q(0, 0) is accurate at large sig (fixed in 0.5.0)", {
  expect_equal(dbipoilog(0, 0, -2, -2, 5, 5, 0.5), exp(ref_b$lq[3]), tolerance = 1e-9)
})

test_that("very unequal counts at high correlation are symmetric (fixed in 0.5.0)", {
  a <- dbipoilog(0, 500, 3, 3, 2, 2, 0.95, log = TRUE)
  b <- dbipoilog(500, 0, 3, 3, 2, 2, 0.95, log = TRUE)
  expect_equal(a, b, tolerance = 1e-10)
})

test_that("rho = 1 and rho = -1 work (failed in 0.4.2.1)", {
  # at rho = +-1 both log intensities are driven by the same normal deviate u, so the
  # pmf is the single integral of phi(u) dpois(x, exp(mu1 + sig1 u)) dpois(y, exp(mu2 +-
  # sig2 u)), computed here with integrate()
  one <- function(x, y, mu1, mu2, s1, s2, r)
    integrate(function(u) dnorm(u) * dpois(x, exp(mu1 + s1 * u)) * dpois(y, exp(mu2 + r * s2 * u)),
              -12, 12, rel.tol = 1e-12, subdivisions = 1000)$value
  expect_equal(dbipoilog(4, 4, 1, 1, 1, 1, 1), one(4, 4, 1, 1, 1, 1, 1), tolerance = 1e-9)
  expect_equal(dbipoilog(4, 2, 1, 1, 1, 1, -1), one(4, 2, 1, 1, 1, 1, -1), tolerance = 1e-9)
  expect_equal(dbipoilog(0, 7, 0.5, 1.5, 0.8, 1.3, -1), one(0, 7, 0.5, 1.5, 0.8, 1.3, -1), tolerance = 1e-9)
})

test_that("rho = 0 factorises into the univariate pmf", {
  expect_equal(dbipoilog(3, 5, 1, 0, 1, 1.5, 0),
               dpoilog(3, 1, 1) * dpoilog(5, 0, 1.5), tolerance = 1e-9)
})

test_that("log P(N1 + N2 > 0) agrees with the reference values", {
  lpos <- poilog:::lbipoilogpos
  # q(0, 0) below 1/2, where it is computed as log1p(-q(0, 0))
  expect_equal(lpos(-1, -1, 2, 2, 0.5), -0.502297460512365, tolerance = 1e-9)
  expect_equal(lpos(0, 1, 1, 1.5, -0.5), -0.0514817816679212, tolerance = 1e-9)
  # q(0, 0) above 1/2, where it is integrated directly
  expect_equal(lpos(-3, -3, 2, 2, 0.5), -1.35509039070593, tolerance = 1e-9)
  # on both sides of q(0, 0) = 1/2 (0.5047 and 0.4954), where the two forms meet
  expect_equal(lpos(-1.62, -1.5692, 2, 2, 0.5), -0.702544527645374, tolerance = 1e-9)
  expect_equal(lpos(-1.52, -1.5692, 2, 2, 0.5), -0.684007042437691, tolerance = 1e-9)
})

test_that("log P(N1 + N2 > 0) is accurate where q(0, 0) is small and rho is close to -1", {
  # q(0, 0) = 7.0e-7 here, so log P(N1 + N2 > 0) = log1p(-q(0, 0)) is about -7.0e-7, and
  # an error of q(0, 0) relative to q(0, 0) carries over to it unchanged. The reference
  # is log1p(-q(0, 0)) with q(0, 0) from the reference pmf; the bound, 1e-12 absolute,
  # is 1.4e-6 relative to the value.
  lp <- poilog:::lbipoilogpos(-0.6286, 8.1117, 0.6698, 2.1814, -0.99988)
  expect_lt(abs(lp - (-7.0418507082729953e-07)), 1e-12)
  # and it is smooth in rho there: along 201 points 1e-6 apart, the fourth differences
  # of a smooth function are of order step^4 (here far below 1e-20), while a jump of
  # size J in the values gives fourth differences of up to 3 J; so max |fourth
  # difference| / 3 bounds the largest jump.
  rr <- -0.99988 + (-100:100) * 1e-6
  v <- vapply(rr, function(r) poilog:::lbipoilogpos(-0.6286, 8.1117, 0.6698, 2.1814, r), 0)
  expect_lt(max(abs(diff(v, differences = 4))) / 3, 1e-12)
})

test_that("dbipoilog rejects invalid input", {
  expect_error(dbipoilog(1:3, 1:2, 0, 0, 1, 1, 0.5), "equal length")
  expect_error(dbipoilog(c(1, -1), 1:2, 0, 0, 1, 1, 0.5), "negative")
  expect_error(dbipoilog(1, 3e9, 0, 0, 1, 1, 0.5), "integer.max")
  expect_error(dbipoilog(1, 1, 0, 0, 0, 1, 0.5), "sig1")
  expect_error(dbipoilog(1, 1, 0, 0, 1, 1, 1.5), "rho")
  expect_error(dbipoilog(1, 1, numeric(0), 0, 1, 1, 0.5), "single number")
  expect_error(dbipoilog(1, 1, 0, 0, 1, 1, NA_real_), "finite")
})

test_that("sig1 and sig2 up to 1000 are supported, larger values are refused", {
  expect_true(is.finite(dbipoilog(1, 2, 0, 0, 1000, 1000, 0.5, log = TRUE)))
  expect_true(is.finite(poilog:::lbipoilogpos(0, 0, 1000, 1, -1)))
  expect_error(dbipoilog(1, 1, 0, 0, 1, 2000, 0.5), "sig2 larger than 1000 is not supported")
  expect_error(dbipoilog(1, 1, 0, 0, 1e300, 1, 0.5), "sig1 larger than 1000 is not supported")
})
