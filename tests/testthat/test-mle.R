# A sparse simulated data set (32 species observed out of 80, mostly singletons)
# on which poilogMLE 0.4.2.1 stopped at a false optimum (mu = -46, sig = 11,
# about 1e6 species). The maximum of the zero-truncated log-likelihood is at
# mu = -4.4943, sig = 2.7565, with log-likelihood -57.4128679001; it was found
# by a profile search with an implementation that does not use poilog
# (tests/reference/make-reference.R).
sparse <- c(rep(1, 18), rep(2, 6), 4, 5, 5, 6, 6, 7, 16, 45)
sparse_max <- -57.4128679001

test_that("poilogMLE reaches the maximum on sparse data (fixed in 0.5.0)", {
  f <- poilogMLE(sparse)
  expect_equal(f$logLval, sparse_max, tolerance = 1e-6)
  expect_equal(unname(f$par), c(-4.4943, 2.7565), tolerance = 1e-3)
  expect_true(length(sparse) / f$p < 400)
})

test_that("poilogMLE's log-likelihood is the zero-truncated log-likelihood", {
  f <- poilogMLE(sparse)
  ll <- sum(dpoilog(sparse, f$par[1], f$par[2], log = TRUE)) -
    length(sparse) * poilog:::lpoilogpos(f$par[1], f$par[2])
  expect_equal(f$logLval, ll, tolerance = 1e-10)
})

test_that("the truncated likelihood is accurate far out on the ridge (fixed in 0.5.0)", {
  # Far out on the ridge mu -> -Inf, sig -> Inf, q0 is so close to 1 that
  # log(1 - q0) computed by subtraction is log(0), and the zero-truncated
  # log-likelihood came out as +Inf, far above the true maximum. The reference value
  # is from tests/reference/make-reference.R.
  tab <- table(sparse)
  ll <- -poilog:::nll1(c(-250, log(25)), as.numeric(names(tab)), as.numeric(tab), TRUE)$value
  expect_equal(ll, -62.6199996667863, tolerance = 1e-9)
  expect_lt(ll, sparse_max)
})

test_that("the analytic gradients agree with finite differences", {
  nd <- function(f, z, h = 1e-5)
    vapply(seq_along(z), function(j) { e <- replace(0 * z, j, h); (f(z + e) - f(z - e)) / (2 * h) }, 0)
  un1 <- c(0, 1, 2, 3, 7, 20); nr1 <- c(5, 9, 4, 3, 2, 1)
  for (zt in c(TRUE, FALSE)) for (z in list(c(0.3, log(1.2)), c(-2, log(3)))) {
    u <- if (zt) un1[-1] else un1; r <- if (zt) nr1[-1] else nr1
    f <- function(z) poilog:::nll1(z, u, r, zt)$value
    expect_equal(poilog:::nll1(z, u, r, zt)$gradient, nd(f, z), tolerance = 1e-6)
  }
  un2 <- cbind(c(0, 1, 0, 2, 5), c(1, 0, 3, 2, 1)); nr2 <- c(4, 3, 2, 2, 1)
  for (zt in c(TRUE, FALSE)) for (z in list(c(0.2, -0.3, log(1.1), log(0.8), 0.6),
                                            c(-1, 0.5, log(2.5), log(1.5), -1.2))) {
    f <- function(z) poilog:::nll2(z, un2, nr2, zt)$value
    expect_equal(poilog:::nll2(z, un2, nr2, zt)$gradient, nd(f, z), tolerance = 1e-6)
  }
})

test_that("the bootstrap returns and gof is a proportion (it never returned in 0.4.2.1)", {
  set.seed(1)
  invisible(utils::capture.output(f <- poilogMLE(sparse, nboot = 5)))
  expect_equal(dim(f$boot), c(5L, 3L))
  expect_named(f$boot, c("mu", "sig", "logLval"))
  expect_length(f$gof, 1)
  expect_true(f$gof >= 0 && f$gof <= 1)
})

test_that("gof is the share of bootstrap log-likelihoods below the observed one", {
  expect_equal(poilog:::boot_gof(5, c(1, 2, 3)), 1)
  expect_equal(poilog:::boot_gof(0, c(1, 2, 3)), 0)
  expect_equal(poilog:::boot_gof(2, c(1, 2, 2, 3)), 0.25)
})

test_that("the bootstrap stops when refits keep failing", {
  expect_error(poilog:::run_boot(3, simulate = function() 1, refit = function(s) NULL,
                                 cols = "a", report = function(...) NULL),
               "only 0 successful refits")
})

test_that("bipoilogMLE fits and its bootstrap returns", {
  set.seed(2)
  xy <- rbipoilog(S = 30, mu1 = 1, mu2 = 1, sig1 = 1, sig2 = 1, rho = 0.5)
  f <- bipoilogMLE(xy)
  expect_true(all(is.finite(f$par)))
  # a start on the other side of rho = 0 reaches the same maximum
  expect_equal(bipoilogMLE(xy, startVals = c(0, 0, 1, 1, -0.5))$logLval, f$logLval, tolerance = 1e-6)
  set.seed(3)
  invisible(utils::capture.output(fb <- bipoilogMLE(xy, nboot = 2)))
  expect_equal(nrow(fb$boot), 2L)
  expect_named(fb$boot, c("mu1", "mu2", "sig1", "sig2", "rho", "logLval"))
  expect_true(fb$gof >= 0 && fb$gof <= 1)
})

test_that("bipoilogMLE's log-likelihood is the zero-truncated log-likelihood", {
  set.seed(4)
  xy <- rbipoilog(S = 25, mu1 = 0, mu2 = 0.5, sig1 = 2, sig2 = 1.5, rho = 0.3)
  f <- bipoilogMLE(xy)
  pr <- f$par
  ll <- sum(dbipoilog(xy[, 1], xy[, 2], pr[["mu1"]], pr[["mu2"]], pr[["sig1"]], pr[["sig2"]],
                      pr[["rho"]], log = TRUE)) -
    nrow(xy) * poilog:::lbipoilogpos(pr[["mu1"]], pr[["mu2"]], pr[["sig1"]], pr[["sig2"]], pr[["rho"]])
  expect_equal(f$logLval, ll, tolerance = 1e-6)
})

# Seven species; the maximum of the zero-truncated log-likelihood is at mu = 0.8629,
# sig = 1.0341, with log-likelihood -17.2185247833 (tests/reference/make-reference.R).
few <- c(1, 1, 2, 3, 5, 8, 13)

test_that("poilogMLE reaches the maximum from a large start value of sig", {
  f <- poilogMLE(few, startVals = c(mu = 1, sig = 90))
  expect_equal(f$logLval, -17.2185247833, tolerance = 1e-7)
  expect_equal(unname(f$par), c(0.8629, 1.0341), tolerance = 1e-3)
})

test_that("start values of sig beyond the supported range are refused", {
  expect_error(poilogMLE(few, startVals = c(mu = 1, sig = 1e150)),
               "start value of sig larger than 1000 is not supported")
  expect_error(bipoilogMLE(few, rev(few), startVals = c(1, 1, 2, 1e4, 0.5)),
               "start value of sig2 larger than 1000 is not supported")
})

test_that("data whose likelihood has no maximum are refused before the search", {
  # every count 1 (zero-truncated) or every count 0 (not truncated): the likelihood
  # increases towards a limit as mu goes to -Inf. With the default control the fit
  # stops with the cause named; with a loose reltol the optimiser stops close to the
  # limit, where the end-point test alone does not recognise it.
  msg <- "the likelihood has no maximum for these data"
  expect_error(poilogMLE(rep(1, 20)), msg)
  expect_error(poilogMLE(rep(1, 20), control = list(reltol = 1e-2)), msg)
  expect_error(poilogMLE(rep(0, 20), zTrunc = FALSE, control = list(reltol = 1e-2)), msg)
  expect_error(bipoilogMLE(rep(0, 5), rep(0, 5), zTrunc = FALSE), msg)
  expect_error(bipoilogMLE(c(1, 0, 1), c(0, 1, 0)), msg)
  # in a bootstrap such a simulated data set is a failed refit, which is replaced
  expect_null(poilog:::boot_refit(c(0, 0), NULL, rep(1, 5), list(), "BFGS",
                                  poilog:::MODEL1, TRUE))
})

test_that("a fit that does not end at a maximum stops with an error", {
  # 68 species simulated with rho = 0.5; from the start rho = -0.99999 BFGS ends near
  # rho = -1 at a point 71.9 below the maximum that is not a maximum (a saddle of the
  # likelihood), and optim reports success there. From rho = 0 the fit reaches the
  # maximum, logL -383.1706.
  sim2 <- function(S, mu1, mu2, s1, s2, rho) {
    z1 <- rnorm(S); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(S)
    xy <- cbind(rpois(S, exp(mu1 + s1 * z1)), rpois(S, exp(mu2 + s2 * z2)))
    xy[rowSums(xy) > 0, , drop = FALSE]
  }
  set.seed(1616)
  xy <- sim2(80, 1, 0.5, 1.5, 1.2, 0.5)
  expect_error(bipoilogMLE(xy, startVals = c(1, 1, 2, 2, -0.99999)),
               "did not end at a maximum of the likelihood")
  expect_equal(bipoilogMLE(xy, startVals = c(1, 1, 2, 2, 0))$logLval, -383.1706, tolerance = 1e-6)
})

test_that("bipoilogMLE reaches rho beyond 0.9999 (rho was limited to 0.9999 before 0.5.0)", {
  # 100 species with large counts and rho = 0.99999: the maximum is at rho = 0.99994,
  # beyond the old limit, where the likelihood used to be flat in rho; two starts
  # reach the same maximum
  set.seed(5)
  z1 <- rnorm(100); z2 <- 0.99999 * z1 + sqrt(1 - 0.99999^2) * rnorm(100)
  xy <- cbind(rpois(100, exp(8 + z1)), rpois(100, exp(8 + z2)))
  f <- bipoilogMLE(xy, startVals = c(8, 8, 1, 1, 0.9))
  expect_gt(f$par[["rho"]], 0.9999)
  expect_equal(bipoilogMLE(xy, startVals = c(8, 8, 1, 1, 0.5))$logLval, f$logLval, tolerance = 1e-9)
})

test_that("a univariate SANN fit gets no gradient and ends near the maximum", {
  # optim's SANN would use a gradient function as its generator of candidate points;
  # given the gradient, the fit ends far from the maximum (and stops with an error)
  set.seed(1)
  f <- poilogMLE(sparse, method = "SANN")
  expect_gt(f$logLval, sparse_max - 0.05)
})

test_that("during the search a point the C code cannot evaluate counts as +Inf", {
  # A stand-in objective in z = (mu, log sig) with its minimum at mu = 1, sig = 2. It
  # fails like the C code (an error whose message starts with "poilog:") for mu > 3,
  # where the first BFGS step from the start lands.
  obj <- function(msg, from = 3) list(
    fn = function(z, nr, un) {
      if (z[1] > from) stop(msg)
      10 * sum((z - c(1, log(2)))^2)
    },
    gr = function(z, nr, un) 20 * (z - c(1, log(2))))
  tab <- list(un = 1, nr = 1)
  fit <- function(ob, z = c(0, 0))
    poilog:::fit_ml(z, ob, tab, list(maxit = 1000), "BFGS", poilog:::MODEL1, zTrunc = FALSE)
  f <- fit(obj("poilog: integration did not terminate"))
  expect_equal(f$par, c(1, log(2)), tolerance = 1e-6)
  # any other error stops the fit
  expect_error(fit(obj("something else")), "something else")
  # and so does a failure of the C code at the start values
  expect_error(fit(obj("poilog: integration did not terminate"), z = c(4, 0)), "^poilog:")
})

test_that("a bootstrap refit that fails is rejected, any other error stops the bootstrap", {
  ob <- poilog:::nll_obj(function(z, nr, un) poilog:::nll1(z, un, nr, TRUE))
  refit <- function(ob, sim, maxit = 1000)
    poilog:::boot_refit(c(0, 0), ob, sim, list(maxit = maxit), "BFGS", poilog:::MODEL1, TRUE)
  # a refit that does not converge (here: one iteration allowed) is rejected
  expect_null(refit(ob, sparse, maxit = 1))
  expect_false(is.null(refit(ob, sparse)))
  bad <- list(fn = function(z, ...) stop("not a fit failure"), gr = function(z, ...) z)
  expect_error(refit(bad, sparse), "not a fit failure")
})

test_that("a bootstrap data set with a count above .Machine$integer.max is redrawn", {
  i <- 0
  simulate <- function() { i <<- i + 1; if (i == 1) c(1, 3e9) else c(1, 2) }
  expect_message(b <- poilog:::run_boot(2, simulate, refit = function(s) sum(s), cols = "a",
                                        report = function(...) NULL),
                 "1 simulated data set(s) had a count above .Machine$integer.max", fixed = TRUE)
  expect_equal(b[, "a"], c(3, 3))
})

test_that("bipoilogMLE stops on unequal lengths (it recycled in 0.4.2.1)", {
  expect_error(bipoilogMLE(1:5, 1:4), "equal length")
})

test_that("poilogMLE recovers the parameters from a large sample", {
  set.seed(20)
  x <- rpoilog(3000, 1, 1, keep0 = TRUE)
  x <- x[x > 0]
  f <- poilogMLE(x)
  expect_equal(unname(f$par), c(1, 1), tolerance = 0.1)
  expect_equal(f$p, 1 - dpoilog(0, 1, 1), tolerance = 0.05)
})

test_that("poilogMLE with zTrunc = FALSE uses the zeros", {
  set.seed(21)
  x <- rpoilog(2000, 0.5, 1, keep0 = TRUE)
  f <- poilogMLE(x, zTrunc = FALSE)
  expect_equal(unname(f$par), c(0.5, 1), tolerance = 0.15)
  expect_equal(f$logLval, sum(dpoilog(x, f$par[1], f$par[2], log = TRUE)), tolerance = 1e-10)
})

test_that("poilogMLE gives the same fit from a distant start, for reordered data and with Nelder-Mead", {
  f1 <- poilogMLE(sparse)
  f2 <- poilogMLE(sparse, startVals = c(mu = 3, sig = 0.5))
  f3 <- poilogMLE(rev(sparse))
  f4 <- poilogMLE(sparse, method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-12))
  expect_equal(f2$logLval, f1$logLval, tolerance = 1e-6)
  expect_equal(f3$logLval, f1$logLval, tolerance = 1e-10)
  expect_equal(f4$logLval, f1$logLval, tolerance = 1e-6)
})

test_that("bipoilogMLE accepts a two-column matrix or two vectors", {
  set.seed(22)
  xy <- rbipoilog(S = 40, mu1 = 1, mu2 = 0.5, sig1 = 1, sig2 = 1, rho = 0.3)
  a <- bipoilogMLE(xy)
  b <- bipoilogMLE(xy[, 1], xy[, 2])
  expect_equal(a$logLval, b$logLval, tolerance = 1e-10)
  expect_equal(a$par, b$par, tolerance = 1e-8)
  # a matrix as n2 is refused (its columns would otherwise be taken as extra rows)
  expect_error(bipoilogMLE(c(xy[, 1], xy[, 1]), cbind(xy[, 2], xy[, 2])), "n2 must be a vector")
  # a matrix as n1 together with n2 is refused (n2 would otherwise be ignored)
  expect_error(bipoilogMLE(xy, xy[, 2]), "n1 is a matrix and n2 is specified")
  expect_error(bipoilogMLE(xy, startVals = c(1, 1, 0, 2, 0.5)), "start value of sig1")
  expect_error(bipoilogMLE(xy, startVals = c(1, 1, 2, -1, 0.5)), "start value of sig2")
})

test_that("the MLE functions reject invalid counts", {
  # invalid counts are rejected before anything reaches the C code
  expect_error(poilogMLE(c(1, 2, -1)), "negative")
  expect_error(poilogMLE(c(1, 2.5, 3)), "integers")
  expect_error(poilogMLE(c(1, NA, 3)), "NA")
  expect_error(poilogMLE(numeric(0)), "empty")
  expect_error(poilogMLE(c(1, 3e9)), "integer.max")
  expect_error(bipoilogMLE(c(1, 2, -1), c(1, 0, 2)), "negative")
  expect_error(bipoilogMLE(c(1, 2, 3), c(1, 0.5, 2)), "integers")
})

test_that("zero-truncated fits reject zero counts", {
  expect_error(poilogMLE(c(0, 1, 2, 3, 5)), "zTrunc = FALSE")
  expect_error(bipoilogMLE(c(0, 1, 2), c(0, 0, 1)), "zTrunc = FALSE")
})

test_that("the MLE functions reject invalid arguments", {
  expect_error(poilogMLE(1:3, zTrunc = NA), "TRUE or FALSE")
  expect_error(poilogMLE(1:3, nboot = -1), "non-negative integer")
  expect_error(poilogMLE(1:3, nboot = 1.5), "non-negative integer")
  expect_error(poilogMLE(1:3, startVals = c(mu = NA, sig = 1)), "finite")
  expect_error(poilogMLE(1:3, startVals = c("1", "2")), "startVals must be numeric")
  expect_error(poilogMLE(1:3, startVals = c(1, 2, 3)), "length of startVals is not 2")
  expect_error(poilogMLE(1:3, startVals = c(mu = 1, sig = 0)), "sig")
  expect_error(poilogMLE(1:3, method = "Brent"), "method must be")
  expect_error(poilogMLE(1:3, method = "L-BFGS-B"), "not supported")
  expect_error(poilogMLE(cbind(1:3, 1:3)), "bipoilogMLE")
  expect_error(bipoilogMLE(1:5), "n2 is missing")
  expect_error(bipoilogMLE(cbind(1:3, 1:3, 1:3)), "2 columns")
  expect_error(bipoilogMLE(1:3, 1:3, startVals = c(1, 1, 1, 1, 1)), "rho")
})

test_that("bipoilogMLE reaches the maximum at large sig (fixed in 0.5.0)", {
  # 34 species simulated with sig1 = sig2 = 4. At such sig, 0.4.2.1 computed q(0, 0) too
  # low, so its log-likelihood was wrong and it underestimated the number of species.
  # The maximum of the zero-truncated log-likelihood, found with an implementation that
  # does not use poilog (tests/reference/make-reference.R), is at mu1 = -0.689029,
  # mu2 = -0.036663, sig1 = 3.767654, sig2 = 3.518454, rho = 0.673021, with
  # log-likelihood -266.8836338857. The likelihood is flat in mu2, so the fit, which
  # stops within optim's default relative tolerance, can differ there by a few 1e-3.
  n1 <- c(6, 0, 455, 1, 33, 51, 0, 0, 34, 1, 102, 1, 0, 4, 2, 580, 4, 0, 0, 0, 1, 6,
          237, 1, 231, 0, 86, 8, 1, 0, 0, 9, 123, 6)
  n2 <- c(0, 3, 5, 0, 123, 0, 21, 11, 1, 1, 2578, 0, 4, 16, 2, 38, 138, 3, 11, 4, 2, 1,
          71, 0, 14, 8, 2059, 438, 0, 1, 2, 23, 31, 2)
  f <- bipoilogMLE(n1, n2)
  expect_equal(f$logLval, -266.8836338857, tolerance = 1e-7)
  expect_lt(max(abs(unname(f$par) - c(-0.689029, -0.036663, 3.767654, 3.518454, 0.673021))), 0.01)
})
