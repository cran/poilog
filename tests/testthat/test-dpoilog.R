# Reference values: log q(n; mu, sig) and log P(N > 0), computed by
# tests/reference/make-reference.R with an implementation that does not use
# poilog (tests/reference/ref.R; that directory is in the source repository,
# not in the built package).

ref_u <- data.frame(
  n   = c(1, 3, 10, 100, 2, 5, 0, 0, 0, 1000, 1e6),
  mu  = c(0, 1, 2, 4, 3, 10, -2, 0, -8, 0, 10),
  sig = c(1, 1.5, 1, 1, 0.5, 2, 5, 8, 3, 0.5, 1),
  lq  = c(-1.35148288213476, -2.50248456069446, -3.30110818905671,
          -5.70743994059652, -6.4778993298831, -11.8409727627916,
          -0.494829270670746, -0.751010795894318, -0.0131574936959435,
          -102.17207345338, -22.0135008127735)
)
ref_lq_at <- function(n, mu, sig) ref_u$lq[ref_u$n == n & ref_u$mu == mu & ref_u$sig == sig]

test_that("dpoilog agrees with the reference values", {
  for (i in seq_len(nrow(ref_u))) {
    r <- ref_u[i, ]
    expect_equal(dpoilog(r$n, r$mu, r$sig, log = TRUE), r$lq, tolerance = 1e-9,
                 label = sprintf("log dpoilog(%g, %g, %g)", r$n, r$mu, r$sig))
  }
})

test_that("small counts far from the mode are not overestimated (fixed in 0.5.0)", {
  # 0.4.2.1 returned about e times the true value here
  expect_equal(dpoilog(5, 10, 2), exp(ref_lq_at(5, 10, 2)), tolerance = 1e-9)
})

test_that("q(0) is accurate at large sig", {
  expect_equal(dpoilog(0, -2, 5), exp(ref_lq_at(0, -2, 5)), tolerance = 1e-9)
  expect_equal(dpoilog(0, 0, 8), exp(ref_lq_at(0, 0, 8)), tolerance = 1e-9)
})

test_that("log = TRUE matches log() and does not underflow", {
  n <- c(0, 1, 2, 5, 20, 100)
  expect_equal(dpoilog(n, 1, 1.5, log = TRUE), log(dpoilog(n, 1, 1.5)), tolerance = 1e-12)
  # a probability below the smallest positive double (about 4.9e-324, i.e. exp(-744.4))
  # is returned as 0, and only the log scale can represent it
  lq <- dpoilog(1e5, -5, 0.3, log = TRUE)
  expect_true(is.finite(lq) && lq < -745)
  expect_identical(dpoilog(1e5, -5, 0.3), 0)
})

test_that("the probabilities sum to 1", {
  expect_equal(sum(dpoilog(0:500, 0, 1)), 1, tolerance = 1e-9)
})

test_that("log P(N > 0) is accurate along the ridge mu -> -Inf, sig -> Inf (fixed in 0.5.0)", {
  # 0.4.2.1 computed log(1 - q0) by subtraction, which rounds to log(0) here
  expect_equal(poilog:::lpoilogpos(-200, 25), -34.7250063765372, tolerance = 1e-9)
  expect_equal(poilog:::lpoilogpos(0, 1), -0.480872829174658, tolerance = 1e-9)
  expect_equal(poilog:::lpoilogpos(-5, 2), -3.335141361168, tolerance = 1e-9)
})

test_that("dpoilog rejects invalid counts", {
  expect_error(dpoilog(-1, 0, 1), "negative")
  expect_error(dpoilog(1.5, 0, 1), "integers")
  expect_error(dpoilog(c(1, NA), 0, 1), "NA")
  expect_error(dpoilog("3", 0, 1), "numeric")
  expect_error(dpoilog(3e9, 0, 1), "integer.max")
})

test_that("dpoilog rejects invalid parameters", {
  expect_error(dpoilog(1, numeric(0), 1), "single number")
  expect_error(dpoilog(1, c(0, 1), 1), "single number")
  expect_error(dpoilog(1, NA_real_, 1), "finite")
  expect_error(dpoilog(1, 0, NA), "sig must be finite")
  expect_error(dpoilog(1, 0, Inf), "finite")
  expect_error(dpoilog(1, 0, 0), "larger than 0")
  expect_error(dpoilog(1, 0, 1, log = NA), "TRUE or FALSE")
})

test_that("sig up to 1000 is supported, larger values are refused", {
  expect_true(is.finite(dpoilog(1, 0, 1000, log = TRUE)))
  expect_true(is.finite(poilog:::lpoilogpos(0, 1000)))
  expect_error(dpoilog(1, 0, 1000.001), "sig larger than 1000 is not supported")
  expect_error(dpoilog(1, 0, 1e200), "sig larger than 1000 is not supported")
  expect_error(poilog:::lpoilogpos(0, 1e4), "sig larger than 1000 is not supported")
})
