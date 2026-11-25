\name{MLE for Poisson lognormal distribution}
\alias{poilogMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum Likelihood Estimation for Poisson Lognormal Distribution }
\description{
  \code{poilogMLE} fits the Poisson lognormal distribution to data and estimates 
  parameters mean \code{mu} and standard deviation \code{sig} in the lognormal distribution
}
\usage{
poilogMLE(n, startVals = c(mu=1, sig=2), 
          nboot = 0, zTrunc = TRUE,
          method = "BFGS", control = list(maxit=1000))
}
\arguments{
  \item{n}{ A vector of counts }
  \item{startVals}{ Starting values of parameters, see details }
  \item{nboot}{ Number of parametric bootstraps, defaults to zero }
  \item{zTrunc}{ Zero-truncated distribution, defaults to TRUE }
  \item{method}{ Method to use during optimization, see details }
  \item{control}{ A list of control parameters for the optimization routine, see details }
  
}
\details{
 
  The function estimates parameters mean \code{mu} and standard deviation \code{sig}. 
  In cases of incomplete sampling the estimate of \code{mu} will be confounded with the sampling
  intensity (see \code{\link{rpoilog}}). Assuming sampling intensity \eqn{\nu}, 
  the estimates of the mean is \eqn{\code{mu}+\ln(\nu)}. Parameter \code{sig} can be estimated without 
  any knowledge of sampling intensity.\cr
  The parameters must be given starting values for the optimization procedure (default starting values are 
  used if starting values are not specified in the function call). The start value of \code{sig} must be
  larger than 0 and at most 1000, the range supported by \code{\link{dpoilog}}, and the search stays in that
  range: beyond it the likelihood is treated as not evaluable, and the optimizer moves away.
  
  The function uses the optimization procedures in \code{\link{optim}} to obtain the maximum likelihood estimate.
  The \code{method} and \code{control} arguments are passed to \code{\link{optim}}, see the help page for this
  function for its control parameters. The methods \code{"BFGS"} (the default) and \code{"CG"}
  use the analytic gradient of the log-likelihood; \code{"Nelder-Mead"} and \code{"SANN"} are
  also accepted. \code{"L-BFGS-B"} is not supported, since the function takes no bounds on the
  parameters (without them it would stop at the first point where the likelihood cannot be
  evaluated, where the other methods move away). During the search, a point where the numerical
  integration fails, which happens only at extreme parameter values, is treated as a point where
  the likelihood cannot be evaluated.

  The function returns an estimate only if the fit ends at a maximum of the likelihood. It stops
  with an error, instead, if
  \itemize{
    \item the data have no maximum likelihood estimate: every count is 1 with
      \code{zTrunc = TRUE}, or every count is 0 with \code{zTrunc = FALSE}. The likelihood
      then increases towards a limit as \code{mu} goes to \eqn{-\infty}{-Inf};
    \item the likelihood cannot be evaluated at the start values;
    \item the optimization does not converge (for example when \code{maxit} is reached);
    \item the fit does not end at a maximum of the likelihood. After the optimization has
      converged, the Hessian of the negative log-likelihood in (\code{mu}, \code{sig}) is
      computed from its analytic gradient. The end point counts as a maximum if every eigenvalue
      of the Hessian exceeds \eqn{10^{-6}}{1e-6} (a smaller one means a standard error above
      1000 in some direction, so the data do not determine the parameters) and if a Newton
      step would increase the log-likelihood by less than 0.05. Likely reasons for a failure are
      start values from which the optimizer stopped short of the maximum or at a point that is
      not one, data for which no maximum exists, an estimate
      that approaches \code{sig = 1000}, and the method \code{"SANN"} or loose \code{control}
      settings (a large \code{reltol} or a small \code{maxit}), which can stop short of the
      maximum. Other start values, a larger \code{maxit} or a smaller \code{reltol} may help;
    \item the numerical integration fails at or near the estimate.
  }

  A zero-truncated distribution (see \code{\link{dpoilog}}) is assumed by default (\code{zTrunc = TRUE}).
  In cases where the number of zeros is known the \code{zTrunc} argument should be set to \code{FALSE}.
  The counts must be non-negative integers; with \code{zTrunc = TRUE} they must also be positive.

  The approximate fraction of species revealed by the sample is \eqn{1-q(0;\code{mu},\code{sig})}.

  Parametric bootstrapping is done by simulating new sets of observations using the estimated parameters,
  conditional on the observed number of species (see \code{\link{rpoilog}}), and refitting the model to
  each, starting from the estimates. A refit that does not converge or does not end at a maximum of the
  likelihood (including a simulated data set that has no maximum) is discarded and replaced by a new simulated data set. The bootstrap distribution is
  therefore conditioned on simulated data sets for which the maximum likelihood estimate exists; for
  sparse data, with a few tens of species, about a fifth of the refits can be discarded. Any other error
  in a refit stops the bootstrap. The bootstrap stops with an error if \code{nboot} successful refits
  are not obtained from \code{10*nboot + 50} simulated data sets. A simulated data set with a count
  above \code{.Machine$integer.max}, which cannot be fitted, is replaced by a new one as well, which
  conditions the bootstrap on data sets without such counts; the number of data sets replaced for this
  reason is reported in a message.
  
  
}
\value{
  \item{par }{Maximum likelihood estimates of the parameters}
  \item{p }{Approximate fraction of species revealed by the sample}
  \item{logLval }{Log likelihood of the data given the estimated parameters}
  \item{gof }{Goodness of fit measure: the proportion of bootstrap log-likelihoods that are smaller than
                logLval. (gof<0.05) or (gof>0.95) indicates lack of fit. NULL if \code{nboot = 0}}
  \item{boot }{A data frame with the bootstrap replicates of the parameters (columns \code{mu} and
                \code{sig}) and their log-likelihoods (\code{logLval}). NULL if \code{nboot = 0}}
}
\references{
Bulmer, M. G. 1974. On fitting the Poisson lognormal distribution to species abundance data. Biometrics \bold{30}, 651-660.\cr
Engen, S., R. Lande, T. Walla and P. J. DeVries. 2002. Analyzing spatial structure of communities using the two-dimensional 
            Poisson lognormal species abundance model. American Naturalist \bold{160}, 60-73.}
\author{ Vidar Grotan \email{vidar.grotan@ntnu.no}, Steinar Engen}

\seealso{ \code{\link{optim}}, \code{\link{dpoilog}}, \code{\link{rpoilog}}  }
\examples{
### simulate observations
\dontshow{set.seed(200)}
n <- rpoilog(S=80,mu=1,sig=2)

### obtain estimates of parameters
est <- poilogMLE(n)

### similar, but now with bootstrapping ###
est <- poilogMLE(n,nboot=10)

### change start values and request tracing information 
### from optimization procedure
est <- poilogMLE(n,startVals=c(2,3),
                 control=list(maxit=1000,trace=1, REPORT=1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
