/* Poisson-lognormal (PLN) and bivariate Poisson-lognormal probabilities.
 *
 * For a count x with log intensity my + sd u, u ~ N(0, 1), the pmf is
 *   q(x; my, sd) = int phi(u) dpois(x, exp(my + sd u)) du.
 * The bivariate pmf conditions the second log intensity on the same u: given u it is
 * my2 + ro sd2 u + sd2 c v with v ~ N(0, 1) and c = sqrt(1 - ro^2). So q(x, y) is an
 * outer integral over u of the first Poisson factor times a univariate PLN pmf of the
 * second count (the inner integral, over v). Everything is on the log scale, so tiny
 * probabilities do not underflow.
 *
 * Zero truncation. The zero-truncated fits need log P(N > 0) and log P(N1 + N2 > 0).
 *  - log P(N > 0) is integrated directly, as int phi(u) (1 - exp(-lambda)) du: the form
 *    log(1 - q(0)) cancels when q(0) is close to 1.
 *  - log P(N1 + N2 > 0) is log1p(-q0), q0 = q(0, 0), wherever q0 <= 1/2. There an error
 *    of q0 relative to q0 changes log1p(-q0) by at most 2 q0 times that error, so the
 *    error shrinks with q0. Where q0 > 1/2, log1p(-q0) would cancel, and the direct form
 *      log P(N1 + N2 > 0) = log[ P(N1 > 0) + int phi(u) exp(-lambda1) P(N2 > 0 | u) du ]
 *    is used instead: both terms are positive, so it has no cancellation, and there the
 *    error relative to P(N1 + N2 > 0) < 1/2 is what matters. The direct form is not used
 *    where q0 <= 1/2 for three reasons: its error does not shrink with q0, so where q0 is
 *    small (log P close to 0) it is large relative to log P; where the stopping rule
 *    accepts one level later along a parameter path its value jumps by up to about
 *    1.6e-12 (near rho = -1, with q0 about 7e-7), against 1e-20 for log1p(-q0), and a
 *    likelihood multiplies such a jump by the number of species; and it takes about
 *    2.5 to 3 times as long. At q0 = 1/2 both forms are accurate, so switching between
 *    them adds no visible jump (below 2e-15).
 *  - Both truncation terms (and the q0 of the switch) are integrated one level finer than
 *    the stopping rule requires (see "extra level" below), because they enter a
 *    zero-truncated log-likelihood multiplied by the number of species.
 *
 * Quadrature. Every integrand exp(f(u)) here is log-concave with f'' <= -1, and f is
 * -u^2/2 plus terms of the form x eta - exp(eta) or log(1 - exp(-exp(eta))) in a linear
 * predictor eta = a + b u (or, for the outer integral of the bivariate, a smoothed
 * version of these: the inner integral as a function of its location).
 *  - The mode is found first (univariate pmf: closed form via Lambert W, polished by
 *    Newton; otherwise safeguarded Newton with analytic first and second derivatives).
 *  - Trapezoid rule on nodes centred at the mode, summing outward until the log
 *    integrand is DROP nats below the mode (by log-concavity all further terms are
 *    smaller and decay faster than geometrically), then halving the step, reusing all
 *    nodes. The integrands are entire, so the rule converges geometrically once the
 *    step resolves the integrand's features.
 *  - Change of variables. The exp(eta) terms are sharp, on the scale 1 / |b| in u,
 *    where exp(eta) is of order one (a cliff for a zero count, a kink for P(N > 0)),
 *    while the rest of the integrand varies on the scale of the Gaussian. Uniform steps
 *    would need about 40 |b| nodes, which for large sd (and nested integrals) is too
 *    many. Where a term is sharper than the mode's scale w, the rule works in s with
 *    t = psi(s) analytic and psi' = eps on a window that covers the sharp range of eta
 *    (zmap_* below), so one step size serves both scales. The map decides the cost
 *    only: whether a step is accepted is decided by the stopping rule, from the
 *    integrand's values as a function of t at the nodes the map places.
 *  - Stopping rule. Two successive levels can agree by coincidence while both are
 *    wrong, when both step over a feature (the error at one level is an oscillating
 *    function of where the feature falls between nodes, and can pass through zero).
 *    So agreement alone is not accepted as evidence. A level is accepted only when
 *    (a) the previous step resolved the integrand wherever it has mass, judged from the
 *    node values themselves, and (b) the two levels agree to a relative tol. For (a), on
 *    every three consecutive nodes the second divided difference gives the curvature
 *    kappa = -f''; with depth = f(mode) - f at the middle node:
 *     (G) kappa h^2 <= RESG (1 + depth): the step resolves the curvature;
 *     (E) the feature curvature c = kappa - 1 (the Gaussian factor contributes exactly
 *         1) changes between neighbouring triples by at most a factor
 *         exp(ARES / (LRES - depth)): the step resolves the rate at which the curvature
 *         changes. For these entire integrands that rate sets the scale l of their
 *         features, and with it how far into the complex plane they stay of moderate
 *         size, which is what sets the trapezoid error (the integrands have no
 *         singularities; the error model is a model of this behaviour, see below).
 *         A feature at relative height e^-depth of scale l, sampled at
 *         step h, contributes an error of order e^-depth exp(-pi^2 l / h); in that
 *         model (E) holds the accepted level (half the resolved step) to
 *         exp(-LRES) = 1e-11 relative, and ignores what lies more than LRES nats down.
 *         The rate is corrected for the feature's small amplitude in the tail of a term
 *         exp(-exp(eta)) (see hist_add).
 *    (a) does not use the agreement in (b), so a coincidental agreement of two levels is
 *    not enough to stop. (b) is needed for the accuracy too: (a) alone can accept a
 *    level that still differs from the previous one by more than tol, and at some
 *    parameter values it then leaves errors of about 5e-11 in log q and 1e-8 in the
 *    scores (the gradient), where with (b) the rule takes at least one more level and
 *    the errors are about 1e-13 and 1e-16. The scores depend on (b) the most.
 *    The error model behind (a) is a model of exp(-exp(eta)) shapes, not a proof. Its
 *    constants are set from that model, tightened where an independent quadrature
 *    shows it to be too lax, and validated against that quadrature on dense parameter
 *    grids, by a search over mu for the largest error, and along fine parameter paths
 *    (on which the largest jump in a returned log probability is about 4e-11, and in a
 *    log-likelihood of 100 to 300 species about 3e-12). The remaining
 *    risk: (a) estimates the scale of a feature from node values, so a feature much
 *    narrower than the step whose total effect on the slope is small could pass unseen,
 *    as could one whose shape departs from the model.
 *  - Extra level. Where a parameter change makes the rule accept one level later, the
 *    returned value moves by the error of the earlier level: a small jump. For pmf terms
 *    this part of the jumps above is below about 4e-14; for the truncation terms it is
 *    about 5e-12, which a likelihood of N species multiplies by N. So the truncation
 *    terms take one more halving after the accepted level: the value then moves by the
 *    error of the accepted level, which is far smaller (below 1e-14), except for the
 *    direct form of log P(N1 + N2 > 0) where q0 is small (1.6e-12, see above), which is
 *    therefore not used there. The inner pmf integrals of log1p(-q0) take the extra
 *    level too: without it the gradient of N log P(N1 + N2 > 0) jumps up to 15 times
 *    more across the switch at q0 = 1/2, and fits are no faster.
 *  - Posterior moments come from the same nodes: derivatives of an inner log integral
 *    with respect to its location (L' = E[g'], L'' = E[g''] + Var[g'], used by the
 *    outer Newton search), and the score E[d log integrand / d theta], which is the
 *    gradient of log q with respect to the parameters (used by the MLEs).
 * Function names: u... functions (ubi, upos, ubipos) evaluate a log integrand at u, with
 * its derivatives in u, for the mode search (type ufun); t... functions (tpmf, tpos, tbi,
 * tbipos) evaluate it at t = u - mode, relative to its value at the mode, with the
 * moments, for the trapezoid rule (type tfun).
 * sd = 0 gives the Poisson log pmf. The R functions allow sd <= 1000, the range over
 * which the kernel was validated; for much larger sd the sharp window becomes narrower
 * than the rounding of t. The .C entry points receive variances and take square roots;
 * the bivariate ones also receive c = sqrt(1 - ro^2), which the caller computes in a form
 * that does not round to 0 when ro is close to +-1, and they return the ro component of
 * the gradient multiplied by c (finite at ro = +-1, where d/dro itself is not). There is
 * no heap allocation, so error() can longjmp without leaking.
 */
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* nodes are summed outward until the log integrand is 36 nats (a factor 2e-16) below
 * its mode; with f'' <= -1 the omitted tail is smaller still, far below rounding */
#define DROP 36.0
/* at most 16 halvings of the first step (h down to 1.5e-5 w) to reach an accepted level,
 * then error() (the truncation terms take one more halving after that, see "extra
 * level"): this only guarantees that a rule that never accepts stops; no input in the
 * supported range needs more than a few halvings */
#define MAXLEV 16
/* nodes per integral, all levels together, then error(). With the change of variables
 * an integral takes a few hundred to a few thousand nodes. BUDGET below would also stop
 * an integral that runs away; this bound stops it much sooner (BUDGET is 800 times as
 * large), so that a failure is reported promptly */
#define MAXNODES 400000
/* nodes per probability, all nested integrals together, then error(): the same
 * guarantee for the nested bivariate integrals, whose inner integrals are each bounded
 * by MAXNODES but whose number is not. The largest number of nodes used by any
 * probability on dense validation grids over the supported range is 3.14e6 (the
 * bivariate P(N1 + N2 > 0) at sig = 1000); the bound is 100 times that, rounded up, so
 * that it stops only a computation that has run away. */
#define BUDGET 3.2e8
#define NMOM 7       /* g', g'', and up to 5 scores */

/* resolution tests (see the header). RESG: at the mode (depth 0), kappa h^2 <= 2 at the
 * resolved step is kappa h^2 <= 0.5 at the accepted one, where a Gaussian's trapezoid
 * error is e^-39; further down, (G) allows RESG (1 + depth), since an error there is
 * scaled by e^-depth. */
#define RESG 2.0
/* ARES = 2 pi^2 (19.739, rounded to 19.74): with the error model
 * e^-depth exp(-pi^2 l / h) at the accepted step h (half the resolved one), a rate of
 * change below ARES / (LRES - depth) per resolved step keeps the error below
 * exp(-LRES) */
#define ARES 19.74
#define LRES 25.0    /* nats: exp(-25) = 1.4e-11, the relative error aimed at */
/* curvature floor for (E): feature curvatures below 5% of the Gaussian's (which is 1)
 * are too weak to shape the integrand on the scale of a step and do not count */
#define EPSF 0.05
/* (b), the agreement of two levels, relative: TOL1 for univariate integrals (dpoilog,
 * P(N > 0)) and TOLI for the inner integrals of the bivariate. (b) is part of the
 * accuracy, not only a check (see the header): where (a) accepts but two levels still
 * differ by more than tol, the rule takes another level. With both tests the accepted
 * (finer) level is far more accurate than tol itself (the errors in the header), since
 * the rule converges geometrically once the step resolves the integrand. The values
 * were set together with the constants of (a) and validated with them; a change of tol
 * changes the results, the scores most. */
#define TOL1 1e-7
#define TOLI 1e-7
/* outer integrals of the bivariate: each outer node carries an inner integral, whose
 * agreement test bounds its error only to about TOLI, and these errors enter the outer
 * sum; so two outer levels cannot be required to agree much more closely than that.
 * (Changing TOLO to 1e-7 moves the bivariate scores by up to about 1e-10.) */
#define TOLO 1e-6
/* first step, in units of the local scale w: on a Gaussian of scale w the trapezoid
 * rule with step h has relative error about 2 exp(-2 pi^2 w^2 / h^2), 1e-8 at h = w and
 * far below rounding at h = w / 2, so for a near-Gaussian integrand one or two halvings
 * suffice */
#define H0 1.0
/* nodes used by the current probability, for BUDGET. A file-level counter: each .C
 * entry point resets it before each probability. (R calls the entry points from one
 * thread, so it is not shared between computations.) */
static double WORK = 0;

/* for the Poisson log pmf (lpois_eta) and the mode of the pmf integrand
 * (pmf_mode_guess): below eta = -700 (exp(-700) = 1e-304, just above the smallest normal
 * double, 2.2e-308) these are formed from the logarithm instead of from exp(), so that
 * they neither underflow nor lose precision. Elsewhere exp() of a smaller argument is
 * evaluated and allowed to underflow to 0 where 0 is the correct limit to double
 * precision (log1mexp_exp, s_pos, and the weights of the trapezoid nodes). */
#define LOGTINY -700.0

static void poilog1(int *x, double *my, double *sig, int *nrN, double *val);
static void poilog2(int *x, int *y, double *my1, double *my2, double *sig1, double *sig2,
                    double *ro, double *c, int *nrN, double *val);
static void poilogpos1(double *my, double *sig, double *val);
static void poilogpos2(double *my1, double *my2, double *sig1, double *sig2, double *ro,
                       double *c, double *val);
static void poilog1g(int *x, double *my, double *sig, int *nrN, double *val, double *gr);
static void poilog2g(int *x, int *y, double *my1, double *my2, double *sig1, double *sig2,
                     double *ro, double *c, int *nrN, double *val, double *gr);
static void poilogpos1g(double *my, double *sig, double *val, double *gr);
static void poilogpos2g(double *my1, double *my2, double *sig1, double *sig2, double *ro,
                        double *c, double *val, double *gr);

static const R_CMethodDef cMethods[] = {
  {"poilog2", (DL_FUNC) &poilog2, 10},
  {"poilog1", (DL_FUNC) &poilog1, 5},
  {"poilogpos1", (DL_FUNC) &poilogpos1, 3},
  {"poilogpos2", (DL_FUNC) &poilogpos2, 7},
  {"poilog2g", (DL_FUNC) &poilog2g, 11},
  {"poilog1g", (DL_FUNC) &poilog1g, 6},
  {"poilogpos1g", (DL_FUNC) &poilogpos1g, 4},
  {"poilogpos2g", (DL_FUNC) &poilogpos2g, 8},
  {NULL, NULL, 0}
};

void R_init_poilog(DllInfo *dll)
{
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

/* ---------------------------------------------------------------- small tools */

/* expm1(y) - y without cancellation. For |y| < 0.1 a Taylor series through y^10
 * (relative truncation error below 1e-16); above, expm1(y) - y loses at most a
 * factor 2/|y| <= 20 ulps to cancellation. */
static double em1x(double y)
{
  if (fabs(y) < 0.1)
    return y * y * (1.0 / 2 + y * (1.0 / 6 + y * (1.0 / 24 + y * (1.0 / 120 +
           y * (1.0 / 720 + y * (1.0 / 5040 + y * (1.0 / 40320 + y * (1.0 / 362880 +
           y / 3628800))))))));
  return expm1(y) - y;
}

/* log(1 - exp(-lambda)) for lambda = exp(eta), accurate for all eta. Below eta = -30
 * (lambda < 1e-13) the series eta - lambda/2 is exact to double precision (the next
 * term is lambda^2/24 < 1e-27) and stays finite where exp(eta) underflows. Above, the
 * standard switch at lambda = log 2 between log1p(-exp(-lambda)) and log(-expm1(-lambda))
 * keeps full relative accuracy on both sides. */
static double log1mexp_exp(double eta)
{
  if (eta < -30.0) return eta - 0.5 * exp(eta);
  double lam = exp(eta);
  return (lam > M_LN2) ? log1p(-exp(-lam)) : log(-expm1(-lam));
}

/* s(eta) = d/deta log(1 - exp(-lambda)) = lambda / expm1(lambda), lambda = exp(eta).
 * Below eta = -30 the series 1 - lambda/2 (next term lambda^2/12 < 1e-27) avoids 0/0.
 * Above lambda = 700, exp(-lambda) < 1e-304, so s = lambda exp(-lambda) = exp(eta -
 * lambda) to double precision; computed that way it is 0, not Inf * 0, when lambda
 * overflows. */
static double s_pos(double eta, double *lam_out)
{
  double lam = exp(eta);
  *lam_out = lam;
  if (eta < -30.0) return 1.0 - 0.5 * lam;
  if (lam > 700.0) return exp(eta - lam);
  return lam / expm1(lam);
}

/* ds/deta = s (1 - lambda - s); 0 when s = 0 (lambda infinite), where the product
 * would be 0 * Inf */
static double ds_pos(double s, double lam)
{
  return (s > 0.0) ? s * (1.0 - lam - s) : 0.0;
}

static double logaddexp(double a, double b)
{
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  return (a > b) ? a + log1p(exp(b - a)) : b + log1p(exp(a - b));
}

/* log Poisson pmf at lambda = exp(eta) */
static double lpois_eta(int x, double eta)
{
  return (eta > LOGTINY) ? dpois_raw((double) x, exp(eta), TRUE)
                         : x * eta - exp(eta) - lgammafn(x + 1.0);
}

/* Lambert W (principal branch) of z = exp(lz). Newton from a one-sided start converges
 * monotonically and quadratically, typically in 3 to 6 steps; it stops when the step
 * is below 1e-15 relative (a few ulps), and 100 steps is only a safety cap. The caller,
 * pmf_mode_guess, handles lz < LOGTINY itself, so lz >= LOGTINY here. */
static double lambertW_exp(double lz)
{
  int k;
  double w;
  if (lz <= 1.0) {
    /* w e^w = z from w0 = log1p(z) >= W(z): Newton on a convex increasing
     * function, monotone from the right */
    double z = exp(lz);
    w = log1p(z);
    for (k = 0; k < 100; k++) {
      double ew = exp(w), dw = (w * ew - z) / (ew * (1.0 + w));
      w -= dw;
      if (fabs(dw) <= 1e-15 * w) break;
    }
    return w;
  }
  /* w + log w = lz from w0 = lz - log lz <= W: concave increasing, monotone from
   * the left */
  w = lz - log(lz);
  for (k = 0; k < 100; k++) {
    double dw = (w + log(w) - lz) / (1.0 + 1.0 / w);
    w -= dw;
    if (fabs(dw) <= 1e-15 * w) break;
  }
  return w;
}

/* mode of the u-integrand of the pmf, x (my + sd u) - e^(my + sd u) - u^2/2:
 * s e^eta = W(s e^(s x + my)), s = sd^2, so u = sd x - W / sd. For a tiny argument
 * (below e^-700) W(z) = z to double precision. */
static double pmf_mode_guess(int x, double my, double sd)
{
  double s = sd * sd, c = s * x + my, lz = log(s) + c;
  return (lz < LOGTINY) ? sd * (x - exp(c)) : sd * x - lambertW_exp(lz) / sd;
}

/* ---------------------------------------------------------------- trapezoid core */

/* f(t, p, nm, mom): log integrand at t (t = u - mode) relative to its value at the
 * mode. nm = 0: value only. nm >= 2: mom[0] = g', mom[1] = g'' (derivatives with
 * respect to the location of the integrand, for an outer Newton search), and
 * mom[2 .. nm-1] = further quantities whose posterior means are wanted. */
typedef double (*tfun)(double t, void *p, int nm, double *mom);

/* Change of variables t = psi(s) that shrinks the step inside one window [a, b] of t
 * (where a sharp term of the integrand lives) by the factor eps:
 *   psi'(s) = 1 - (1 - eps) W(s),
 *   W(s) = (tanh((s - al) / tau) - tanh((s - be) / tau)) / 2,
 * a plateau of psi' = eps on [al, be] with smooth edges of width tau. psi is analytic
 * (its nearest singularities are at distance pi tau / 2 from the real axis), so the
 * trapezoid rule in s keeps its geometric convergence. The map decides the cost: with
 * it, one step size in s resolves both the sharp term and the rest of the integrand.
 * Whether a step is accepted is decided by the stopping rule, which tests the
 * transformed integrand and was validated empirically (see the header). */
struct zmap { int on; double al, be, eps, tau, C, w, s0; };

static double logcosh(double x)
{
  x = fabs(x);
  return x + log1p(exp(-2.0 * x)) - M_LN2;
}

static double zmap_raw(const struct zmap *Z, double s)
{
  return s - 0.5 * (1.0 - Z->eps) * Z->tau *
             (logcosh((s - Z->al) / Z->tau) - logcosh((s - Z->be) / Z->tau));
}

/* t = psi(s) and lj = log psi'(s) */
static double zmap_t(const struct zmap *Z, double s, double *lj)
{
  if (!Z->on) { *lj = 0.0; return s; }
  double W = 0.5 * (tanh((s - Z->al) / Z->tau) - tanh((s - Z->be) / Z->tau));
  *lj = log1p(-(1.0 - Z->eps) * W);
  return zmap_raw(Z, s) - Z->C;
}

/* window [a, b] in t (the mode at t = 0), step factor eps, coarse scale w. A map is set
 * only if it shrinks the step at least by half (eps < 0.5); otherwise the plain rule is
 * as cheap. The plateau is [0, 2 me + (b - a) / eps] in s, shifted so that psi(me) = a;
 * me = (tau / 2) log(2 / eps) is where an edge has brought psi' down to 2 eps, so
 * psi' <= 2 eps on all of [a, b] (psi(be - me) = b up to O(eps tau)). The grid of the
 * trapezoid rule is centred at s0 = psi^-1(0), found by bisection: psi increases, and
 * |psi(s) - s| is at most the plateau length plus |a| + tau, so [-hi, hi] with
 * hi = 2 (be + |a|) + 10 brackets s0 (the 10 exceeds tau = 2 w <= 2, which is all the
 * bracket needs); 200 bisection steps reach the stopping width 1e-14 relative from any
 * such bracket. */
static void zmap_set(struct zmap *Z, double a, double b, double eps, double w)
{
  int k;
  double lj, lo, hi;
  Z->on = 0;
  if (!(a < b) || !(eps < 0.5)) return;
  Z->on = 1; Z->eps = eps; Z->w = w;
  /* across an edge log psi' changes at the rate 2 / tau, and the nearest poles of psi'
   * are at distance pi tau / 2: tau = 2 w keeps the edges resolved at the coarse step */
  Z->tau = 2.0 * w;
  double me = 0.5 * Z->tau * log(2.0 / eps);
  Z->al = 0.0; Z->be = 2.0 * me + (b - a) / eps;
  Z->C = zmap_raw(Z, me) - a;   /* so that psi(me) = a */
  hi = 2.0 * (Z->be + fabs(a)) + 10.0; lo = -hi;
  for (k = 0; k < 200 && hi - lo > 1e-14 * (1.0 + fabs(lo)); k++) {
    double mid = 0.5 * (lo + hi);
    if (zmap_t(Z, mid, &lj) < 0.0) lo = mid; else hi = mid;
  }
  Z->s0 = 0.5 * (lo + hi);
}

/* window for a term of the log integrand in eta = e0 + b t (t from the mode) whose
 * sharp part lives in eta in [lo, hi] and varies on the scale 1 / beff in t (beff = |b|
 * for an exp(eta) term; smaller for a term smoothed by an inner integral), seen at the
 * coarse scale w; merged with any window already in Z. Nothing is done if the term is
 * smooth on the scale w (beff w <= 2: at the step w / 2 of the first halving eta moves
 * by at most 1 per node, and two more halvings give four nodes per unit of eta, so
 * uniform steps resolve it at little cost; the stopping rule decides how many halvings
 * are taken), or if the window misses the mass of the
 * integrand: |t| <= 9.5 holds all of it, since f(t) <= f(0) - t^2/2 and the nodes stop
 * at f(0) - DROP, i.e. within |t| <= sqrt(2 DROP) = 8.5. */
static void zmap_term(struct zmap *Z, double e0, double b, double beff, double lo,
                      double hi, double km, double w)
{
  if (!(beff * w > 2.0)) return;
  /* coarse scale: that of the integrand without this term (its curvature at the mode
   * is km), the step outside the window; the window brings it down to 1 / beff. In the
   * univariate integrals the term is the only one besides the Gaussian, and km is its
   * curvature at the same point as w, so 1 / w^2 - km = 1 and wr = 1, the Gaussian's
   * scale. For the inner-integral term of the bivariate outer integrand (bi_zoom), km
   * is passed as 0, so wr is the full scale w, which includes that term's own
   * curvature and is the smaller: the coarse step is then finer than needed, which is
   * conservative and costs nothing measurable (passing the curvature instead changes
   * results only by about 1e-10). */
  double wr = 1.0 / sqrt(fmax(1.0, 1.0 / (w * w) - km));
  if (Z->on) wr = fmin(wr, Z->w);
  double ta = (lo - e0) / b, tb = (hi - e0) / b, eps = 1.0 / (beff * wr);
  if (ta > tb) { double z = ta; ta = tb; tb = z; }
  ta = fmax(ta, -9.5); tb = fmin(tb, 9.5);
  if (!(ta < tb)) return;
  if (Z->on) {   /* merge with the window already set, in t */
    double lj, me = 0.5 * Z->tau * log(2.0 / Z->eps);
    double a0 = zmap_t(Z, Z->al + me, &lj), b0 = zmap_t(Z, Z->be - me, &lj);
    ta = fmin(ta, a0); tb = fmax(tb, b0); eps = fmin(eps, Z->eps);
  }
  zmap_set(Z, ta, tb, eps, wr);
}

/* resolution bookkeeping of one sweep: the last three nodes (v = f(t) relative to the
 * mode, v[0] the latest), how many are valid, the feature curvature at the previous
 * triple, the first three nodes of the sweep (to seed the sweep in the opposite
 * direction, so that the nodes across the centre are also checked), and the largest
 * violation r seen (resolved if r <= 1). The tests are on f as a function of t, so
 * they measure the integrand and not the change of variables. */
struct hist { double v[3], t[3], fv[3], ft[3], cp, r; int n, k, hasc; };

static void hist_start(struct hist *H, double r)
{
  H->n = H->k = H->hasc = 0; H->cp = 0.0; H->r = r;
}

static void hist_max(struct hist *H, double x)
{
  if (!(x <= H->r)) H->r = ISNAN(x) ? R_PosInf : x;
}

/* curvature -f'' from three nodes (second divided difference), and the half width ht
 * of the triple */
static double hist_curv(double v0, double v1, double v2, double t0, double t1, double t2,
                        double *ht)
{
  *ht = 0.5 * fabs(t0 - t2);
  return -2.0 * ((v0 - v1) / (t0 - t1) - (v1 - v2) / (t1 - t2)) / (t0 - t2);
}

/* continue on the other side of the centre: the first nodes of sweep R are the
 * neighbours of the first node of the next sweep */
static void hist_turn(struct hist *L, const struct hist *R)
{
  int j;
  double ht;
  hist_start(L, R->r);
  L->n = R->k;   /* at most 3: hist_add keeps only the first three nodes */
  for (j = 0; j < 3; j++) { L->v[j] = R->fv[j]; L->t[j] = R->ft[j]; }
  L->hasc = (L->n == 3);
  if (L->hasc)
    L->cp = fmax(hist_curv(R->fv[2], R->fv[1], R->fv[0], R->ft[2], R->ft[1], R->ft[0],
                           &ht) - 1.0, 0.0);
}

/* add node v = f(t) to the resolution check. Two tests on each triple of nodes, with
 * kappa = -f'' from the triple, ht its half width and depth = -v at its middle node:
 * (G) kappa ht^2 <= RESG (1 + depth): the step resolves the overall curvature;
 * (E) |log((c + EPSF) / (c' + EPSF))| <= ARES / (LRES - depth) for the feature
 *     curvatures c = kappa - 1 >= 0 of neighbouring triples (the Gaussian factor
 *     contributes exactly 1): the step resolves how fast the curvature changes. */
static void hist_add(struct hist *H, double v, double t)
{
  if (H->n >= 2) {
    double ht, kap = hist_curv(v, H->v[0], H->v[1], t, H->t[0], H->t[1], &ht);
    double c = fmax(kap - 1.0, 0.0), d1 = fmax(0.0, -H->v[0]);
    hist_max(H, kap * ht * ht / (RESG * (1.0 + d1)));
    /* depth and the amplitude correction below only relax (E), so a ratio within
     * exp(+-ARES / LRES) = exp(+-0.79) passes without them. This shortcut, which skips
     * the logarithms, uses a bound just inside that (0.999 of it, for rounding), so it
     * never passes a ratio that (E) would fail */
    double q = (c + EPSF) / (H->cp + EPSF), qmax = 0.999 * exp(ARES / LRES);
    if (H->hasc && (q > qmax || q < 1.0 / qmax)) {
      double d = fmin(d1, fmax(0.0, -H->v[1]));
      double rho = fabs(log(q));
      /* a feature of curvature c varying at the rate b = rho / ht adds about c / b^2
       * to f; a small amplitude counts as extra depth */
      double amp = fmax(c, H->cp) * ht * ht / (rho * rho);
      if (amp < 1.0) {
        /* in the tail of a term exp(-e^eta), the curvature changes at about -log(amp)
         * times the rate at which the term itself varies (its scale is 1 in eta), so
         * the rate of the term is rho / (-log(amp)). The divisor la - 1 instead of la
         * keeps a margin for the approximation, and it is used only where la > 2, so
         * that the divisor exceeds 1 and the correction can only relax the test.
         * These choices are empirical: they were set and checked against an
         * independent quadrature (see the header). */
        double la = -log(amp);
        d += la;
        if (la > 2.0) rho /= la - 1.0;
      }
      if (d < LRES) hist_max(H, rho * (LRES - d) / ARES);
    }
    H->cp = c; H->hasc = 1;
  }
  H->v[2] = H->v[1]; H->v[1] = H->v[0]; H->v[0] = v;
  H->t[2] = H->t[1]; H->t[1] = H->t[0]; H->t[0] = t;
  if (H->n < 3) H->n++;
  if (H->k < 3) { H->fv[H->k] = v; H->ft[H->k] = t; H->k++; }
}

/* sum the nodes s0, s0 + ds, ... outward until the log integrand is DROP below the
 * mode. The test is on f(t), which is unimodal in s because psi is increasing. The
 * counts for MAXNODES and BUDGET leave out the last node of each sweep (the one below
 * DROP). */
static void sweep(tfun f, void *p, const struct zmap *Z, double s0, double ds, int nm,
                  double *S, int *nodes, struct hist *H)
{
  int n, j;
  double g[NMOM];
  for (n = 0; ; n++) {
    double lj, t = zmap_t(Z, s0 + n * ds, &lj), v = f(t, p, nm, g);
    /* v is relative to the mode, so it is <= 0 up to the accuracy of the mode; a
     * value near the overflow threshold of exp (709.8) means the evaluation failed */
    if (ISNAN(v) || v > 700.0) error("poilog: integrand evaluation failed");
    hist_add(H, v, t);
    double e = exp(v + lj);
    S[0] += e;
    /* a node whose weight underflows contributes nothing; its moments may be
     * infinite (e.g. lambda = Inf), so skip them rather than form 0 * Inf */
    if (nm && e > 0.0) {
      S[1] += e * g[0]; S[2] += e * (g[1] + g[0] * g[0]);
      for (j = 2; j < nm; j++) S[j + 1] += e * g[j];
    }
    if (v < -DROP) break;
    if (++(*nodes) > MAXNODES || ++WORK > BUDGET)
      error("poilog: integration did not terminate");
  }
}

/* one level: the nodes halfway between those of step h, on both sides of the centre.
 * R and L receive the resolution check of the new nodes (spacing h). */
static void halve(tfun f, void *p, const struct zmap *Z, double s0, double h, int nm,
                  double *S, int *nodes, struct hist *R, struct hist *L)
{
  hist_start(R, 0.0); sweep(f, p, Z, s0 + 0.5 * h, h, nm, S, nodes, R);
  hist_turn(L, R);    sweep(f, p, Z, s0 - 0.5 * h, -h, nm, S, nodes, L);
}

/* log int exp(f(t)) dt, f with its mode at t = 0, w = 1 / sqrt(-f''(0)),
 * through the map Z (Z->on = 0: none). If nm > 0: M[0] = E[g'], M[1] = E[g''] + Var[g'],
 * M[j] = E[mom[j]] (j >= 2) under the normalised integrand.
 * Level lev (step h0 2^-lev in s) adds the nodes halfway between the previous ones;
 * these new nodes have spacing h0 2^-(lev-1), the previous step, so the resolution
 * test on them says whether the previous step resolved the integrand.
 * extra = 1: return the level after the accepted one (see "extra level" in the
 * header), for the truncation terms. */
static double trap_log(tfun f, void *p, double w, const struct zmap *Z,
                       double tol, int extra, int nm, double *M)
{
  double S[NMOM + 1], T, Tprev, h;
  int lev, j, nodes = 0;
  struct hist R, L;
  /* the first step, in s: w, or with a map its coarse scale (eps times that inside
   * the window) */
  h = H0 * (Z->on ? Z->w : w);
  for (j = 0; j <= NMOM; j++) S[j] = 0.0;
  double s0 = Z->on ? Z->s0 : 0.0;   /* the mode, in s */
  /* level 0. sweep also keeps the resolution check, but its result for level 0 is not
   * used: the first test is at level 1, and halve starts the check afresh */
  hist_start(&R, 0.0); sweep(f, p, Z, s0, h, nm, S, &nodes, &R);
  hist_turn(&L, &R);   sweep(f, p, Z, s0 - h, -h, nm, S, &nodes, &L);
  T = h * S[0];
  for (lev = 1; lev <= MAXLEV; lev++) {
    halve(f, p, Z, s0, h, nm, S, &nodes, &R, &L);
    h *= 0.5;
    Tprev = T; T = h * S[0];
    /* L.r is the worst resolution ratio of both sweeps (hist_turn starts the left
     * sweep from the right sweep's value), so this tests the whole level */
    if (L.r <= 1.0 && fabs(T - Tprev) <= tol * T) break;
  }
  if (lev > MAXLEV) error("poilog: integration did not converge");
  if (extra) {
    halve(f, p, Z, s0, h, nm, S, &nodes, &R, &L);
    h *= 0.5;
    T = h * S[0];
  }

  if (nm) {
    M[0] = S[1] / S[0]; M[1] = S[2] / S[0] - M[0] * M[0];
    for (j = 2; j < nm; j++) M[j] = S[j + 1] / S[0];
  }
  return log(T);
}

/* a log integrand in u: returns f(u), and if d is not NULL d[0] = f'(u), d[1] = f''(u) */
typedef double (*ufun)(double u, void *p, double *d);

/* a better first point for newton_mode than u (where f is fu), using uc, a point where f
 * is finite: move halfway towards uc as long as f is -Inf at the current point, or rises.
 * f is concave, so where it is finite is an interval, which contains uc, and on the
 * segment from u to uc f is first -Inf (if it is at u), then rises up to the mode (if
 * the mode lies on the segment) and then falls. So the moves reach a finite point, and
 * they stop near the mode, within a factor 2 of its distance from uc. This matters when
 * u lies where a term of f falls steeply, like -exp(sd u): there Newton moves only 1/sd
 * per step, and from a start far into that region (or where f overflows to -Inf) it
 * would not reach the mode. When u is near the mode, the first move already falls, and
 * u is kept, at the cost of one evaluation. At most 200 halvings (each halves the
 * distance to uc, so 200 reach uc to rounding from any double). Called only with a
 * finite uc (newton_mode skips it otherwise). */
static double approach_mode(ufun f, void *p, double u, double fu, double uc)
{
  int k;
  for (k = 0; k < 200; k++) {
    double un = 0.5 * (u + uc), fn;
    if (un == u) break;
    fn = f(un, p, NULL);
    if (ISNAN(fn) || !(fu == R_NegInf || fn > fu)) break;
    u = un; fu = fn;
  }
  return u;
}

/* safeguarded Newton for the mode of a concave f with f'' <= -1. Since f' falls at
 * least at unit rate, one evaluation brackets the mode within |f'(u)|. A Newton step
 * that leaves the bracket, or follows a point worse than the best so far (past the
 * mode, in a region where f falls like -exp(sd u) Newton moves only 1/sd per step), is
 * replaced by bisection, so 200 steps cover any bracket down to rounding. A point where
 * f overflows to -Inf is taken as lying past the mode on its side of the best point.
 * The search needs a first point where f is finite, and not far into a region where f
 * falls steeply; given a point uc where f is finite, approach_mode first looks for a
 * better start between u and uc. The search stops when the step, or the bracket, is
 * below 1e-6 of the local scale: the trapezoid rule is valid for any centre, and f
 * there is then within about 1e-12 of its maximum. It returns the best point evaluated
 * (with its f and scale), so the centre is consistent even if the derivatives (from
 * inner integrals) are noisy. */
static double newton_mode(ufun f, void *p, double u, double uc, double *fm, double *w)
{
  double lo = R_NegInf, hi = R_PosInf, d[2], ub = u, fb = R_NegInf, d2b = -1.0;
  int k;
  for (k = 0; k < 200; k++) {
    double fu = f(u, p, d);
    if (k == 0 && R_FINITE(uc) && !ISNAN(fu)) {
      double ua = approach_mode(f, p, u, fu, uc);
      if (ua != u) { u = ua; continue; }   /* start again from the better point */
    }
    if (fu == R_NegInf && fb > R_NegInf) {
      /* a step far past the mode, where exp() overflows: the mode lies between here
       * and the best point, so shrink the bracket and bisect */
      if (u > ub) hi = fmin(hi, u); else lo = fmax(lo, u);
      u = 0.5 * (lo + hi);
      continue;
    }
    if (!R_FINITE(fu) || ISNAN(d[0])) error("poilog: mode search failed");
    if (!(d[1] <= -1.0)) d[1] = -1.0;
    int worse = !(fu > fb);
    if (!worse) { ub = u; fb = fu; d2b = d[1]; }
    if (d[0] > 0) { lo = fmax(lo, u); hi = fmin(hi, u + d[0]); }
    else          { hi = fmin(hi, u); lo = fmax(lo, u + d[0]); }
    double un = u - d[0] / d[1], tol = 1e-6 / sqrt(-d[1]);
    if (worse || !(un >= lo && un <= hi)) un = 0.5 * (lo + hi);
    if (fabs(un - u) < tol || hi - lo < tol) break;
    u = un;
  }
  if (fb == R_NegInf) error("poilog: mode search failed");
  *fm = fb;
  *w = 1.0 / sqrt(-d2b);
  return ub;
}

/* ranges of eta where the exp(eta) term of a log integrand is sharp, at scale sd.
 * The low end is where its curvature in u, sd^2 exp(eta), is a tenth of the floor EPSF
 * (below that the resolution test ignores it). The high end is where it has removed
 * all mass, for the Poisson term x eta - exp(eta): there exp(eta) = 2 (x + DROP), and
 * the term has fallen more than DROP below its maximum at exp(eta) = x. For the term
 * log(1 - exp(-lambda)) of P(N > 0), lambda = exp(eta), the curvature in u is
 * sd^2 s (lambda - 1 + s), s = lambda / expm1(lambda), which for lambda >= 2 is below
 * sd^2 lambda^2 exp(-lambda); the high end is where that bound has fallen to the same
 * floor, lambda = L + 2 log(lambda) with L = log(sd^2 / (0.1 EPSF)). A window is only
 * used for sd > 2 (zmap_term), where L > 6.6 and the root exceeds 11; there the
 * fixed-point iteration from lambda = L + 3 (within 4 of the root) contracts by the
 * factor 2 / lambda <= 0.21 per step, so four steps leave an error below 0.01 in lambda.
 * Where L < 2 - 2 log 2 (sd < 0.096) the bound never reaches the floor and the iteration
 * has no fixed point; max(lambda, 2) (2 is where lambda^2 exp(-lambda) peaks) keeps
 * the result finite there, where it is not used. A fixed end at lambda = L + 3 would
 * leave the bound at lambda^2 exp(-3) times the floor, which exceeds EPSF itself once
 * lambda > 14 (sd above about 20): there the window would end while the term is still
 * sharp. */
static double sharp_lo(double sd) { return log(0.1 * EPSF / (sd * sd)); }
static double sharp_hi_pois(int x) { return log(2.0 * (x + DROP)); }
static double sharp_hi_pos(double sd)
{
  double L = log(sd * sd / (0.1 * EPSF)), lam = L + 3.0;
  int k;
  for (k = 0; k < 4; k++) lam = L + 2.0 * log(fmax(lam, 2.0));
  return log(fmax(lam, 2.0));
}

/* ---------------------------------------------------------------- univariate pmf */

/* integrand in t = u - m:  r t - A em1x(sd t) - t^2/2, A = exp(eta) at the mode (no
 * large terms, so full relative precision even for counts near 2^31). When A is below
 * 1e-300 or sd t is above 700 (where exp(sd t) is close to overflowing),
 * A em1x(sd t) is formed as exp(eta + sd t) - A (1 + sd t) instead, so that a tiny A
 * times an overflowing exponential is not 0 * Inf. Both forms are exact to rounding, so
 * the switch between them causes no jump. */
struct p1c { double A, eta, r, sd, x, xa, m; };   /* xa = x - A */

static double tpmf(double t, void *vp, int nm, double *mom)
{
  struct p1c *p = (struct p1c *) vp;
  double st = p->sd * t;
  int direct = (p->A > 1e-300 && st < 700.0);
  if (nm) {
    double lam;
    if (direct) {
      double e = expm1(st);
      mom[0] = p->xa - p->A * e;         /* g' = x - lambda */
      lam = p->A * (1.0 + e);
    } else {
      lam = exp(p->eta + st);
      mom[0] = p->x - lam;
    }
    mom[1] = -lam;                       /* g'' = -lambda */
    if (nm > 2) mom[2] = mom[0] * (p->m + t);   /* (x - lambda) u: the sd-score */
  }
  double c = direct ? p->A * em1x(st) : exp(p->eta + st) - p->A * (1.0 + st);
  return p->r * t - c - 0.5 * t * t;
}

/* log q(x; my, sd). nm = 2 or 3: M[0] = d/dmy log q, M[1] = d2/dmy2 log q,
 * M[2] = d/dsd log q = E[(x - lambda) u].
 * rel (if not NULL) receives log q - lpois(x, my), the log pmf relative to the Poisson
 * log pmf at the location my, formed without that (possibly huge) term: the bivariate
 * outer integrand needs differences of the inner log pmf between two locations, and
 * the Poisson part of such a difference can then be formed in t without that term,
 * where that is the more accurate form (see tbi). */
static double lpoilog_m(int x, double my, double sd, double tol, int extra, int nm,
                        double *M, double *rel)
{
  double s = sd * sd;
  if (s == 0.0) {
    /* Poisson: sd = 0 (sd is never negative: a square root, or sd2 c), or sd^2
     * underflows, where the mixing is below rounding */
    if (nm) {
      double lam = exp(my);
      M[0] = x - lam; M[1] = -lam; if (nm > 2) M[2] = 0.0;
    }
    if (rel) *rel = 0.0;
    return lpois_eta(x, my);
  }
  double m = pmf_mode_guess(x, my, sd);
  /* the closed form cancels for large s x; polish with Newton steps on
   * d(u) = x sd - sd e^(my + sd u) - u, whose rounding is relative to the slope.
   * Newton converges quadratically from the closed form, so at most 4 steps; the
   * mode is good enough when the step is below 1e-6 of the local scale sqrt(w2). */
  double eta = my + sd * m, A = exp(eta), r = x * sd - sd * A - m;
  double w2 = 1.0 / (s * A + 1.0);
  int k;
  for (k = 0; k < 4 && fabs(r) * w2 > 1e-6 * sqrt(w2); k++) {
    m += r * w2;
    eta = my + sd * m; A = exp(eta); r = x * sd - sd * A - m; w2 = 1.0 / (s * A + 1.0);
  }
  struct p1c p = { A, eta, r, sd, (double) x, x - A, m };
  struct zmap Z = { 0 };
  zmap_term(&Z, eta, sd, sd, sharp_lo(sd), sharp_hi_pois(x), s * A, sqrt(w2));
  double T = trap_log(tpmf, &p, sqrt(w2), &Z, tol, extra, nm, M);
  if (rel) {
    /* lpois(x, eta) - lpois(x, my) = x sd m - (A - lambda0), lambda0 = exp(my), with
     * A - lambda0 = lambda0 expm1(sd m) (or as a difference, under the conditions of
     * tpmf); +Inf where lambda0 overflows, so that lpois(x, my) is not finite */
    double st = sd * m, l0 = exp(my);
    double dA = (l0 > 1e-300 && st < 700.0) ? l0 * expm1(st) : A - l0;
    *rel = R_FINITE(l0) ? x * st - dA - 0.5 * m * m - M_LN_SQRT_2PI + T : R_PosInf;
  }
  return lpois_eta(x, eta) - 0.5 * m * m - M_LN_SQRT_2PI + T;
}

/* ------------------------------------------------------------ P(N > 0), univariate */

struct ppos { double my, sd, m, fm; };

static double upos(double u, void *vp, double *d)
{
  struct ppos *p = (struct ppos *) vp;
  double eta = p->my + p->sd * u, lam, s = s_pos(eta, &lam);
  if (d) { d[0] = p->sd * s - u; d[1] = p->sd * p->sd * ds_pos(s, lam) - 1.0; }
  return log1mexp_exp(eta) - 0.5 * u * u;
}

static double tpos(double t, void *vp, int nm, double *mom)
{
  struct ppos *p = (struct ppos *) vp;
  double u = p->m + t;
  if (nm) {
    double lam, s = s_pos(p->my + p->sd * u, &lam);
    mom[0] = s; mom[1] = ds_pos(s, lam);
    if (nm > 2) mom[2] = s * u;
  }
  return upos(u, p, NULL) - p->fm;
}

/* log P(N > 0); moments as for lpoilog_m. The mode search starts at the mode of the
 * pmf of a count of 1: for small lambda, log(1 - exp(-lambda)) = eta - lambda/2 + ...
 * is close to the Poisson log pmf of 1, eta - lambda, and for large lambda both terms
 * are flat or falling, so that start is near the mode. P(N > 0) enters only the
 * truncation terms, so it always takes the extra level (see the header). */
static double lpoilog_pos_m(double my, double sd, double tol, int nm, double *M)
{
  if (sd * sd == 0.0) {   /* as in lpoilog_m */
    double lam, s = s_pos(my, &lam);
    if (nm) { M[0] = s; M[1] = ds_pos(s, lam); if (nm > 2) M[2] = 0.0; }
    return log1mexp_exp(my);
  }
  struct ppos p = { my, sd, 0.0, 0.0 };
  double w;
  struct zmap Z = { 0 };
  p.m = newton_mode(upos, &p, pmf_mode_guess(1, my, sd), R_NaN, &p.fm, &w);
  double lam, sm = s_pos(my + sd * p.m, &lam);
  zmap_term(&Z, my + sd * p.m, sd, sd, sharp_lo(sd), sharp_hi_pos(sd),
            -sd * sd * ds_pos(sm, lam), w);
  return p.fm + trap_log(tpos, &p, w, &Z, tol, 1, nm, M) - M_LN_SQRT_2PI;
}

/* ---------------------------------------------------------------- bivariate */

/* u: standard normal for the first log intensity eta1 = my1 + sd1 u; given u the
 * second is eta2 = my2 + ro sd2 u + sdc v, sdc = sd2 c, c = sqrt(1 - ro^2), v ~ N(0, 1).
 * Scores of the log integrand: a1 = x - lambda1, a2 = y - lambda2,
 *   d/dmy1 = a1, d/dsd1 = a1 u, d/dmy2 = a2, d/dsd2 = a2 (ro u + c v),
 *   c d/dro = a2 sd2 (c u - ro v)   (d/dro itself is a2 sd2 (u - ro v / c)).
 * The inner integral over v returns E[a2 | u] and E[a2 v | u].
 * extra: the truncation path (see the header), passed on to the inner integrals (the
 * inner P(N2 > 0 | u) of the direct form always takes the extra level).
 * At the mode m (bi_centre): e1 = eta1, A1 = lambda1, xa1 = x - A1; e2 = my2 + ro sd2 m,
 * the location of the inner integral, and A2 = exp(e2); L, the inner log integral, and
 * for the pmf R = L - lpois(y, e2) (see tbi), for the trapezoid nodes in t = u - m (tbi,
 * tbipos). */
struct p2 { int x, y, extra; double my1, my2, sd1, sd2, ro, sdc, m, fm, e1, A1, xa1, e2,
            A2, L, R; };

static double ubi(double u, void *vp, double *d)
{
  struct p2 *p = (struct p2 *) vp;
  double eta = p->my1 + p->sd1 * u, M[3];
  double L = lpoilog_m(p->y, p->my2 + p->ro * p->sd2 * u, p->sdc, TOLI, p->extra,
                       d ? 2 : 0, M, NULL);
  if (d) {
    double lam = exp(eta), b = p->ro * p->sd2;
    d[0] = p->sd1 * (p->x - lam) - u + b * M[0];
    d[1] = -p->sd1 * p->sd1 * lam - 1.0 + b * b * M[1];
  }
  return lpois_eta(p->x, eta) - 0.5 * u * u + L;
}

/* The outer integrand at the trapezoid nodes, in t = u - m and relative to its value
 * at the mode, formed term by term from the differences, as tpmf does for the univariate
 * pmf: the Gaussian -(u^2 - m^2) / 2 = -t (m + t / 2); for the first intensity
 * lambda1(u) - A1 = A1 expm1(sd1 t) (formed as exp(e1 + sd1 t) - A1 under the conditions
 * explained at tpmf), so x - lambda1 = xa1 - A1 expm1(sd1 t); and the inner integral L(u)
 * - L(m). Computing f(u) - f(m) from the two values instead would lose the digits of
 * the large terms, which are of the order of m^2 / 2, A1 and L: rounding errors of these
 * would then appear in every node, as noise in the scores and, when they are large
 * (|m| of 1e6 or more, for large my1 and small sd1), as noise that the stopping rule
 * cannot resolve. */
static double bi_lam1(const struct p2 *p, double st, double *dlam)
{
  if (p->A1 > 1e-300 && st < 700.0) {
    double e = expm1(st);
    *dlam = p->A1 * e;
    return p->A1 * (1.0 + e);
  }
  double lam = exp(p->e1 + st);
  *dlam = lam - p->A1;
  return lam;
}

static double tbi(double t, void *vp, int nm, double *mom)
{
  struct p2 *p = (struct p2 *) vp;
  double u = p->m + t, st = p->sd1 * t, dlam, R, dL, M[3];
  /* the inner log pmf L(loc), loc = my2 + ro sd2 u, relative to its value at the mode,
   * in one of two forms that are equal up to rounding: (a) L(u) - L(m); (b) with
   * L = lpois(y, loc) + R(loc), the Poisson part of the difference, y s2t - (lambda2 -
   * A2), s2t = ro sd2 t, formed as for the first intensity, plus R(u) - R(m). Where the
   * Poisson part carries most of L (the mixing over v changes little: sdc^2 lambda2
   * small, as for a large my2 with a small sd2), L itself can be huge and (a) would lose
   * its digits, while R is small; where the mixing dominates, R and the Poisson part
   * are the huge, cancelling terms instead. Each node uses the form whose terms are
   * smaller in magnitude (so whose rounding error is smaller). */
  double L = lpoilog_m(p->y, p->my2 + p->ro * p->sd2 * u, p->sdc, TOLI, p->extra,
                       nm ? 3 : 0, M, &R);
  double s2t = p->ro * p->sd2 * t, dlam2 = (p->A2 > 1e-300 && s2t < 700.0) ?
    p->A2 * expm1(s2t) : exp(p->e2 + s2t) - p->A2;
  double ea = fabs(L) + fabs(p->L);
  double eb = fabs(p->y * s2t) + fabs(dlam2) + fabs(R) + fabs(p->R);
  dL = (eb < ea) ? (p->y * s2t - dlam2) + (R - p->R) : L - p->L;
  double lam = bi_lam1(p, st, &dlam);
  if (nm) {
    double a1 = (p->A1 > 1e-300 && st < 700.0) ? p->xa1 - dlam : p->x - lam;
    /* mom[0], mom[1] (derivatives in the location, for an outer search) are not used
     * for an outer integral; the scores follow */
    mom[0] = mom[1] = 0.0;
    mom[2] = a1; mom[3] = a1 * u; mom[4] = M[0]; mom[5] = u * M[0]; mom[6] = M[2];
  }
  return p->x * st - dlam - t * (p->m + 0.5 * t) + dL;
}

/* windows for the outer integrand: the term in eta1 = my1 + sd1 u (Poisson for count
 * x, or exp(-lambda1) for P(N1 + N2 > 0)), and the inner integral as a function of
 * its location my2 + ro sd2 u. The inner integral smooths the exp(eta2) term over
 * v ~ N(0, sdc^2), so it varies on the scale sqrt(1 + sdc^2) / |ro sd2| in u, and
 * its sharp range is widened by 3 sdc on each side (3 standard deviations of the
 * smoothing, which hold 99.7% of its mass). The windows only set the cost; the
 * stopping rule checks the result. Uses e1, A1 and e2 from bi_centre. */
static void bi_zoom(struct zmap *Z, const struct p2 *p, int pos, double w)
{
  double b2 = p->ro * p->sd2, sw = 3.0 * p->sdc;
  zmap_term(Z, p->e1, p->sd1, p->sd1, sharp_lo(p->sd1), sharp_hi_pois(p->x),
            p->sd1 * p->sd1 * p->A1, w);
  if (b2 != 0.0)
    zmap_term(Z, p->e2, b2, fabs(b2) / sqrt(1.0 + p->sdc * p->sdc),
              sharp_lo(fabs(b2)) - sw,
              (pos ? sharp_hi_pos(fabs(b2)) : sharp_hi_pois(p->y)) + sw, 0.0, w);
}

/* the start of the outer mode search: the mode of the first factor alone (Poisson
 * term for count x and the Gaussian), and the point uc where the second factor's
 * intensity is about y + 1, the centre of its Poisson term (for P(N1 + N2 > 0), y = 0:
 * intensity 1, the centre of the kink of 1 - exp(-lambda2)). The mode of the sum of the
 * two concave terms lies between the modes of the two terms, and the integrand is
 * finite at uc. At the first point the second intensity can be huge (or overflow), so
 * that the integrand falls steeply there; approach_mode then moves towards uc, as long
 * as the integrand rises. Without a dependence on u (ro sd2 = 0) the second factor is
 * constant in u, so no uc is needed. */
static double bi_start(const struct p2 *p, double *uc)
{
  double b2 = p->ro * p->sd2;
  *uc = (b2 != 0.0) ? (log(p->y + 1.0) - p->my2) / b2 : R_NaN;
  return pmf_mode_guess(p->x, p->my1, p->sd1);
}

/* the quantities at the mode for tbi and tbipos; L is computed afresh (a value-only
 * inner integral at m), not as a difference of fm and the other terms, which would lose
 * the digits that the t-form keeps */
static void bi_centre(struct p2 *p, int pos)
{
  p->e1 = p->my1 + p->sd1 * p->m;
  p->A1 = exp(p->e1);
  p->xa1 = p->x - p->A1;
  p->e2 = p->my2 + p->ro * p->sd2 * p->m;
  p->A2 = exp(p->e2);
  if (pos) p->L = lpoilog_pos_m(p->e2, p->sdc, TOLI, 0, NULL);
  else p->L = lpoilog_m(p->y, p->e2, p->sdc, TOLI, p->extra, 0, NULL, &p->R);
}

/* map posterior means (E a1, E a1 u, E a2, E a2 u, E a2 v) to the gradient in
 * (my1, my2, sd1, sd2, ro), the last one multiplied by c (finite at c = 0) */
static void bi_grad(const double *E, double sd2, double ro, double c, double *g)
{
  g[0] = E[0];
  g[1] = E[2];
  g[2] = E[1];
  g[3] = ro * E[3] + c * E[4];
  g[4] = sd2 * (c * E[3] - ro * E[4]);
}

static double lbipoilog_g(int x, int y, double my1, double my2, double sd1, double sd2,
                          double ro, double c, int extra, double *grad)
{
  struct p2 p = { .x = x, .y = y, .extra = extra, .my1 = my1, .my2 = my2, .sd1 = sd1,
                  .sd2 = sd2, .ro = ro, .sdc = sd2 * c };   /* the rest: bi_centre */
  double w, uc, M[NMOM];
  struct zmap Z = { 0 };
  double u0 = bi_start(&p, &uc);
  p.m = newton_mode(ubi, &p, u0, uc, &p.fm, &w);
  bi_centre(&p, 0);
  bi_zoom(&Z, &p, 0, w);
  double l = p.fm + trap_log(tbi, &p, w, &Z, TOLO, extra, grad ? 7 : 0, M) -
             M_LN_SQRT_2PI;
  if (grad) bi_grad(M + 2, sd2, ro, c, grad);
  return l;
}

/* the direct form of log P(N1 + N2 > 0) (see the header):
 *   log [ P(N1 > 0) + int phi(u) exp(-lambda1) P(N2 > 0 | u) du ];
 * in the second term the scores are -lambda1 (for eta1) and s(eta2) (for eta2) */
static double ubipos(double u, void *vp, double *d)
{
  struct p2 *p = (struct p2 *) vp;
  double eta = p->my1 + p->sd1 * u, M[3], lam = exp(eta);
  double L = lpoilog_pos_m(p->my2 + p->ro * p->sd2 * u, p->sdc, TOLI, d ? 2 : 0, M);
  if (d) {
    double b = p->ro * p->sd2;
    d[0] = -p->sd1 * lam - u + b * M[0];
    d[1] = -p->sd1 * p->sd1 * lam - 1.0 + b * b * M[1];
  }
  return -lam - 0.5 * u * u + L;
}

/* as tbi (in t, relative to the mode, term by term). The inner log P(N2 > 0 | u) is
 * differenced directly: it lies between log(1 - exp(-lambda2)) for no mixing and 0, so
 * it is never a huge number built from cancelling terms as the inner log pmf can be. */
static double tbipos(double t, void *vp, int nm, double *mom)
{
  struct p2 *p = (struct p2 *) vp;
  double u = p->m + t, st = p->sd1 * t, dlam, M[3];
  double L = lpoilog_pos_m(p->my2 + p->ro * p->sd2 * u, p->sdc, TOLI, nm ? 3 : 0, M);
  double lam = bi_lam1(p, st, &dlam);
  if (nm) {
    mom[0] = mom[1] = 0.0;   /* not used, as in tbi */
    mom[2] = -lam; mom[3] = -lam * u; mom[4] = M[0]; mom[5] = u * M[0]; mom[6] = M[2];
  }
  return -dlam - t * (p->m + 0.5 * t) + (L - p->L);
}

static double lbipoilog_pos_direct(double my1, double my2, double sd1, double sd2,
                                   double ro, double c, double *grad)
{
  struct p2 p = { .x = 0, .y = 0, .extra = 1, .my1 = my1, .my2 = my2, .sd1 = sd1,
                  .sd2 = sd2, .ro = ro, .sdc = sd2 * c };   /* the rest: bi_centre */
  double w, uc, M[NMOM], M1[3];
  int nm = grad ? 7 : 0;
  struct zmap Z = { 0 };
  double u0 = bi_start(&p, &uc);
  p.m = newton_mode(ubipos, &p, u0, uc, &p.fm, &w);
  bi_centre(&p, 1);
  bi_zoom(&Z, &p, 1, w);
  double l01 = p.fm + trap_log(tbipos, &p, w, &Z, TOLO, 1, nm, M) - M_LN_SQRT_2PI;
  double l1 = lpoilog_pos_m(my1, sd1, TOL1, grad ? 3 : 0, M1);
  double l = logaddexp(l1, l01);
  if (grad) {
    double g01[5], w1 = exp(l1 - l), w01 = exp(l01 - l);
    int j;
    bi_grad(M + 2, sd2, ro, c, g01);
    for (j = 0; j < 5; j++) grad[j] = w01 * g01[j];
    grad[0] += w1 * M1[0];
    grad[2] += w1 * M1[2];
  }
  return l;
}

/* log P(N1 + N2 > 0): log1p(-q0) where q0 = q(0, 0) <= 1/2, the direct form where
 * q0 > 1/2 (see the header). q0 and its gradient come from the pmf at the cell (0, 0),
 * on the truncation path (extra level). The gradient of log1p(-q0) is
 * -q0 / (1 - q0) times that of log q0. */
static double lbipoilog_pos_g(double my1, double my2, double sd1, double sd2, double ro,
                              double c, double *grad)
{
  double l0 = lbipoilog_g(0, 0, my1, my2, sd1, sd2, ro, c, 1, grad);
  if (l0 <= -M_LN2) {
    double l = log1p(-exp(l0));
    if (grad) {
      double r = exp(l0 - l);
      int j;
      for (j = 0; j < 5; j++) grad[j] *= -r;
    }
    return l;
  }
  return lbipoilog_pos_direct(my1, my2, sd1, sd2, ro, c, grad);
}

/* ---------------------------------------------------------------- .C entry points
 * sig, sig1, sig2 are variances (R passes sig^2), c = sqrt(1 - rho^2). val receives log
 * probabilities; gr the gradients with respect to (mu, sig) or (mu1, mu2, sig1, sig2,
 * rho), with sig the standard deviation and the rho component multiplied by c, one row
 * of the column-major nrN x k matrix per count. */

static void poilog2(int *x, int *y, double *my1, double *my2, double *sig1, double *sig2,
                    double *ro, double *c, int *nrN, double *val)
{
  int i;
  for (i = 0; i < *nrN; i++) {
    WORK = 0;
    val[i] = lbipoilog_g(x[i], y[i], *my1, *my2, sqrt(*sig1), sqrt(*sig2), *ro, *c, 0,
                         NULL);
    R_CheckUserInterrupt();
  }
}

static void poilog2g(int *x, int *y, double *my1, double *my2, double *sig1, double *sig2,
                     double *ro, double *c, int *nrN, double *val, double *gr)
{
  int i, j, n = *nrN;
  double g[5];
  for (i = 0; i < n; i++) {
    WORK = 0;
    val[i] = lbipoilog_g(x[i], y[i], *my1, *my2, sqrt(*sig1), sqrt(*sig2), *ro, *c, 0, g);
    for (j = 0; j < 5; j++) gr[i + j * n] = g[j];
    R_CheckUserInterrupt();
  }
}

static void poilog1(int *x, double *my, double *sig, int *nrN, double *val)
{
  int i;
  for (i = 0; i < *nrN; i++) {
    WORK = 0; val[i] = lpoilog_m(x[i], *my, sqrt(*sig), TOL1, 0, 0, NULL, NULL);
    R_CheckUserInterrupt();
  }
}

static void poilog1g(int *x, double *my, double *sig, int *nrN, double *val, double *gr)
{
  int i, n = *nrN;
  double M[3];
  for (i = 0; i < n; i++) {
    WORK = 0; val[i] = lpoilog_m(x[i], *my, sqrt(*sig), TOL1, 0, 3, M, NULL);
    gr[i] = M[0]; gr[i + n] = M[2];
    R_CheckUserInterrupt();
  }
}

static void poilogpos1(double *my, double *sig, double *val)
{
  WORK = 0; val[0] = lpoilog_pos_m(*my, sqrt(*sig), TOL1, 0, NULL);
}

static void poilogpos1g(double *my, double *sig, double *val, double *gr)
{
  double M[3];
  WORK = 0; val[0] = lpoilog_pos_m(*my, sqrt(*sig), TOL1, 3, M);
  gr[0] = M[0]; gr[1] = M[2];
}

static void poilogpos2(double *my1, double *my2, double *sig1, double *sig2, double *ro,
                       double *c, double *val)
{
  WORK = 0;
  val[0] = lbipoilog_pos_g(*my1, *my2, sqrt(*sig1), sqrt(*sig2), *ro, *c, NULL);
}

static void poilogpos2g(double *my1, double *my2, double *sig1, double *sig2, double *ro,
                        double *c, double *val, double *gr)
{
  WORK = 0;
  val[0] = lbipoilog_pos_g(*my1, *my2, sqrt(*sig1), sqrt(*sig2), *ro, *c, gr);
}
