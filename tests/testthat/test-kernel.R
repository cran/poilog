# Regression tests for the integration kernel (src/bipoilog_s_cint.c). Reference values
# come from an independent quadrature, never from poilog; tests/reference/
# kernel-reference.R recomputes them. Error measure: |log q - ref| / max(1, |ref|).

rel <- function(a, r) abs(a - r) / pmax(1, abs(r))

test_that("log q is accurate where two trapezoid levels can agree by coincidence", {
  # At the first three points two successive trapezoid levels agree before the step
  # resolves the integrand. A rule that accepted agreement alone, without the
  # resolution test (a) of the stopping rule, would be off there by 2.0e-8, 2.6e-9 and
  # 2.3e-9 (the error measure above); with (a) they are accurate to 3e-13. So this test
  # fails if (a) is removed (see the header of the C code). The others are among the
  # least accurate points of dense comparisons with the reference.
  n   <- c(0, 1, 0, 0, 0, 1e5)
  mu  <- c(-5.631, -9.63, -6.89, -51.96, -48.1, -56.15)
  sig <- c(2, 2, 2, 10, 15, 30)
  ref <- c(-0.021223526387426972, -7.6512419682254489, -0.006816544828065263,
           -1.8638508636481532e-07, -0.00080146722018181737, -18.376545137676253)
  for (i in seq_along(n))
    expect_lt(rel(dpoilog(n[i], mu[i], sig[i], log = TRUE), ref[i]), 1e-9)
})

test_that("the gradient of log q is accurate where it needs the level-agreement test", {
  # The stopping rule accepts a level only when the step resolves the integrand (a) and
  # two successive levels agree (b). Here (a) alone accepts too early: the gradient is
  # then off by 1e-8 (d/dmu) and 4.5e-8 (d/dsig); with (b) it is accurate to 1e-15. This
  # test fails if (b) is removed, or if its tolerance TOL1 is loosened to 1e-3; with
  # TOL1 up to 1e-4 the gradient here stays accurate to 1e-15 (see the header of the C
  # code).
  r <- .C("poilog1g", n = 3L, mu = -6.4, sig2 = 4, nrN = 1L, val = 0, gr = double(2),
          PACKAGE = "poilog")
  score_ref <- c(1.6046496654837963, 4.7419816903090535)
  expect_lt(max(abs(r$gr - score_ref)), 1e-10)
})

test_that("log P(N > 0) is accurate across the kink where exp(eta) passes 1", {
  mu  <- c(-18.73, 2.73, -18.09, 3.7072671199841, -2.82612498245936, -2.83)
  sig <- c(3, 1, 15, 50, 5, 5)
  ref <- c(-14.231837429182248, -0.0098540067244330398, -2.1007445092993549,
           -0.6271141907547344, -1.1069502075523654, -1.1077721556312157)
  for (i in seq_along(mu))
    expect_lt(rel(poilog:::lpoilogpos(mu[i], sig[i]), ref[i]), 1e-9)
})

test_that("bivariate probabilities are accurate at large sig", {
  lq <- dbipoilog(c(0, 0, 3), c(0, 1, 2), 1, 0.5, 50, 50, 0.5, log = TRUE)
  ref <- c(-1.1308856974604533, -5.5439875326744987, -11.310120323614075)
  expect_true(all(rel(lq, ref) < 1e-9))
  expect_lt(rel(poilog:::lbipoilogpos(1, 0.5, 50, 50, 0.5), -0.38971077196876003), 1e-9)
})

test_that("values and gradients of log P(N > 0) are finite at large sig", {
  # the gradient involves s = lambda / expm1(lambda), which, formed directly, is
  # Inf / Inf (NaN) once lambda = exp(eta) overflows, as it does at large sig
  pos1 <- function(mu, s) { r <- .C("poilogpos1g", mu = mu, sig2 = s^2, val = 0, gr = double(2),
                                    PACKAGE = "poilog"); c(r$val, r$gr) }
  pos2 <- function(p) { r <- .C("poilogpos2g", mu1 = p[1], mu2 = p[2], sig1 = p[3]^2, sig2 = p[4]^2,
                                rho = p[5], c = sqrt((1 - p[5]) * (1 + p[5])), val = 0, gr = double(5),
                                PACKAGE = "poilog"); c(r$val, r$gr) }
  for (s in c(90, 1000)) {
    a <- pos1(1, s)
    expect_true(all(is.finite(a)))
    h <- 1e-5 * s
    expect_equal(a[3], (pos1(1, s + h)[1] - pos1(1, s - h)[1]) / (2 * h), tolerance = 1e-6)
    expect_equal(a[2], (pos1(1 + 1e-5, s)[1] - pos1(1 - 1e-5, s)[1]) / 2e-5, tolerance = 1e-6)
    p <- c(1, 0, s, 0.8 * s, 0.7)
    b <- pos2(p)
    expect_true(all(is.finite(b)))
    e <- c(0, 0, 1e-5 * s, 0, 0)
    expect_equal(b[4], (pos2(p + e)[1] - pos2(p - e)[1]) / (2e-5 * s), tolerance = 1e-6)
  }
})

test_that("rho = 1 with a very negative mu2 gives a finite log probability", {
  # q(0, 1) = exp(-800) int phi(u) e^u exp(-e^u) du to double precision
  expect_lt(rel(dbipoilog(0, 1, 0, -800, 1, 1, 1, log = TRUE), -801.35148288213475), 1e-9)
  # here a step of the mode search lands so far past the mode that exp() overflows; the
  # search then narrows its bracket and bisects (without that, it stops with "mode
  # search failed")
  expect_lt(rel(dbipoilog(0, 500, 0, -20, 2, 2, 1, log = TRUE), -10012.650644221178), 1e-9)
})
