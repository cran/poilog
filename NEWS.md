# poilog 0.5.0

This release improves the accuracy, robustness and speed of the
computations, and makes the bootstrap in `poilogMLE()` and `bipoilogMLE()`
work again.

`dpoilog()` and `dbipoilog()` gain a `log` argument. It returns the log
probability, which stays accurate even where the probability itself is too
small to be stored and would be returned as 0.

In this work we found numerical errors in earlier versions for very sparse
data and at large `sig`; they are listed below. On the four count data sets
in sads 0.6.5 (bci, birds, moths and okland), 0.4.1, 0.4.2.1 and 0.5.0 give
the same estimates of `mu` and `sig` to four decimals.

Some calls that earlier versions accepted now stop with an error; they are
listed below.

## Improvements

* New integration code computes `dpoilog()`, `dbipoilog()` and the
  probability of a zero count used in zero-truncated fits, for all parameter
  values up to `sig = 1000`. The log probabilities agree with an independent
  reference computation to a relative error of about 1e-10. The
  log-likelihood is also smooth in the parameters, with noise below about
  1e-10, which the optimizer and numerical Hessians need.
* `poilogMLE()` and `bipoilogMLE()` use analytic gradients, and return an
  estimate only if the fit ends at a maximum of the likelihood (`?poilogMLE`
  says how this is checked); otherwise they stop with an error. In
  `bipoilogMLE()`, rho can take any value between -1 and 1.
* The bootstrap (`nboot > 0`) works again; in 0.4.2 and 0.4.2.1 it never
  finished (see below). As in 0.4.1, each simulated data set is refitted
  from the estimates, and a refit that fails is replaced by a new simulated
  data set. A refit now also counts as failed when it does not end at a
  maximum of the likelihood, so the bootstrap uses only refits that are
  estimates (`?poilogMLE`).
* Speed, compared with 0.4.2.1:
  * `poilogMLE()` on `rpoilog(100, 0, 2)` is about 10 times as fast;
  * `bipoilogMLE()` on `rbipoilog(100, 0, 0, 2, 2, 0.8)` about 3 times;
  * `dpoilog(rep(0:999, 5), 1, 2)` about 1.8 times.
* An error during integration can no longer leave memory allocated, because
  the C code no longer allocates memory.
* Documentation corrected:
  * `?dpoilog`: the plain-text (console) version of the formula for the
    integrand was wrong; the typeset version was right.
  * `?dpoilog`: the zero-truncated distribution with sampling intensity `nu`
    now has `mu + ln nu` in place of `mu`.
  * `?dpoilog` and `?dbipoilog`: `mu` and `sig` (and `mu1`, `mu2`, `sig1`,
    `sig2`) are described as the mean and standard deviation of the log
    abundances, and `rho` as the correlation of the log abundances in the two
    samples; `mu2` was described as belonging to sample 1.
  * `?poilogMLE`: the bootstrap refers to `rpoilog()`, not `rbipoilog()`.
  * `?bipoilogMLE`: the sampling-intensity example gives the expected
    estimates of `mu1` and `mu2` as 1 + log(0.5) = 0.31, not 1 - log(0.5).

## Errors in earlier versions, corrected in 0.5.0

* Very sparse data, 0.4.2 and 0.4.2.1: for counts from 0 to 7, `dpoilog()`
  returned too high a probability whenever the correct probability was very
  small (below about 2e-4). On very sparse data (few observed species, most
  seen once or twice), this could make `poilogMLE()` stop at a false
  optimum, with an estimated number of species 35 to 3,000 times too large.
* Large `sig`, 0.4.1 and earlier: `dpoilog()` gave too low a probability of a
  zero count, q(0), at `sig` above about 3, and more so the lower `mu`.
  Zero-truncated fits with `poilogMLE()` use q(0) to estimate the number of
  species. On the four sads data sets (estimated `sig` 1.6 to 2.2), 0.4.1,
  0.4.2.1 and 0.5.0 give the same estimates. In simulations of 100 species
  they also agreed at `sig = 2`. They differed only at `sig = 5` with half
  or more of the species unobserved: there 0.4.1 returned a median of about
  10% fewer species than 0.5.0 on the same data.
* Large `sig`, bivariate, 0.4.2.1 and earlier: the probability that a species
  is seen in neither sample, `dbipoilog(0, 0, ...)`, was too low only at
  large `sig1` and `sig2` combined with low `mu1` and `mu2`. The effect on
  `bipoilogMLE()` was small. In simulations of 100 species at
  `sig1 = sig2 = 4` and `rho = 0.8`, with half the species seen in neither
  sample, 0.4.2.1 returned a median of 3% fewer species than 0.5.0 on the
  same data, and less than 1% at higher `mu`.
* In 0.4.2 and 0.4.2.1 the bootstrap in `poilogMLE()` and `bipoilogMLE()`
  never finished, so it gave no results rather than wrong ones.
* Edge cases (rarely relevant for real data):
  * `dbipoilog()` was too low for very unequal counts at high correlation
    (35% too low for counts 0 and 500 at `rho = 0.95`), and failed at
    `rho = 1` or `-1`.
  * `bipoilogMLE()` kept rho between -0.9999 and 0.9999, which mattered only
    when the maximum lay beyond that limit.
  * Fits started from rho beyond -0.9999 or 0.9999 (for example 0.99995)
    could end at a point that is not a maximum.
  * The zero-truncated `poilogMLE()` overstated the log-likelihood where q(0)
    is close to 1 (around `sig = 25`), which could draw a fit there.

## Changes to the bootstrap results

* `gof` is the proportion of bootstrap log-likelihoods below the observed
  one. 0.4.1 instead returned the number of them below the observed one,
  plus 1, divided by `nboot`, which could exceed 1.
* The column of `boot` in `poilogMLE()` that holds `sig` is named `sig`; in
  0.4.1 it was named `sig2`.

## Calls that now stop with an error

Before, some of these calls gave meaningless results, ran for ever, or
stopped with an unclear message.

* Counts that are not numeric (0.4.2.1 took `TRUE` as 1), negative, not
  integers, missing or larger than `.Machine$integer.max`; parameters that
  are missing, not finite or not single numbers; and invalid `nboot`,
  `startVals` or `method`.
* The logical arguments `log`, `zTrunc`, `condS` and `keep0` when they are not
  `TRUE` or `FALSE`; 0.4.2.1 also accepted, for example, 1 for `condS`,
  `keep0` and `zTrunc`.
* `rpoilog()` with a vector `mu`, which 0.4.2.1 recycled over the draws
  without documenting it.
* `dpoilog()` with `sig` above 1000, `dbipoilog()` with `sig1` or `sig2` above
  1000, and `poilogMLE()` and `bipoilogMLE()` with such start values: the
  probabilities are supported up to `sig = 1000`, the range over which the
  numerical integration has been validated.
* `rpoilog()` and `rbipoilog()` when a Poisson intensity they draw exceeds the
  largest double (a log intensity above 709.78), for which 0.4.2.1 returned
  `NA` counts.
* `dbipoilog()` with `sig1 = 0`, which returned NaN.
* `poilogMLE()` and `bipoilogMLE()` with `method = "L-BFGS-B"`. That method is
  for parameters with bounds, and these functions take none; 0.4.2.1 passed
  it on to `optim()`.
* `bipoilogMLE()` with `n1` and `n2` of different lengths, instead of
  recycling the shorter one.
* `poilogMLE()` with `zTrunc = TRUE` and zero counts, and `bipoilogMLE()` with
  `zTrunc = TRUE` and species with zero counts in both samples, which gave a
  meaningless fit.
* `poilogMLE()` and `bipoilogMLE()` on data for which the likelihood has no
  maximum, so that there is no estimate:
  * in `poilogMLE()`, every count equal to 1 with `zTrunc = TRUE`, or every
    count 0 with `zTrunc = FALSE`;
  * in `bipoilogMLE()`, every species seen exactly once in the two samples
    together with `zTrunc = TRUE`, or in neither sample with
    `zTrunc = FALSE`.

  Such data are fitted ever better by ever rarer species (a sample of
  singletons, for example, as `mu` decreases), so they say nothing about the
  number of species or the shape of the distribution. 0.4.2.1 returned an
  arbitrary point far out (for 20 singletons, `mu = -10.8`), or stopped with
  an integration error.
* `poilogMLE()` and `bipoilogMLE()` when the likelihood cannot be evaluated at
  the start values, where 0.4.2.1 stopped with an unclear error from
  `optim()`.
* `poilogMLE()` and `bipoilogMLE()` when the fit does not end at a maximum of
  the likelihood (see above).
* `rpoilog()` and `rbipoilog()` with `condS = TRUE` and `keep0 = FALSE` when
  positive counts are impossible (`nu = 0`, or `nu1 = nu2 = 0`), where they
  looped for ever, or so rare that more than `max(2e7, 100*S)` draws would be
  needed, which could take gigabytes of memory.
* The bootstrap, when fewer than `nboot` refits succeed among `10*nboot + 50`
  simulated data sets.

## New features

* `dpoilog()` and `dbipoilog()` gain a `log = FALSE` argument. With
  `log = TRUE` they return the log probability, which is accurate also where
  the probability itself is too small to represent in double precision
  (below about 5e-324) and would be returned as 0.

## Testing

* The package now has a test suite (testthat). It checks:
  * the probabilities against independently computed reference values and
    against general properties (sums, means and variances, marginals,
    symmetry, the Poisson limit);
  * the random generators against the probability functions, and their
    seeded draws against those of 0.4.2.1;
  * the maximum-likelihood fits against independently computed maxima, their
    gradients and the bootstrap;
  * the handling of invalid input.

  The scripts that compute the reference values are in `tests/reference/` in
  the source repository (<https://github.com/VidarGrotan/poilog>).
