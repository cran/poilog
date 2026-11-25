\name{MLE for bivariate Poisson lognormal distribution}
\alias{bipoilogMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum Likelihood Estimation for Bivariate Poisson Lognormal Distribution }
\description{
  \code{bipoilogMLE} fits the bivariate Poisson lognormal distribution to data
}
\usage{
bipoilogMLE(n1, n2 = NULL, 
            startVals = c(mu1=1, mu2=1, sig1=2, sig2=2, rho=0.5), 
            nboot = 0, zTrunc = TRUE, file = NULL, 
            method = "BFGS", control = list(maxit=1000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n1}{ a vector or a matrix with two columns of pairwise counts of observed individuals for each species }
  \item{n2}{ if n1 is not given as a matrix, a vector of counts with same ordering of species as in argument n1}
  \item{startVals}{ starting values of parameters }
  \item{nboot}{ number of parametric bootstraps, defaults to zero }
  \item{zTrunc}{ logical; if TRUE (default) the zero-truncated distribution is fitted }
  \item{file}{ text file to hold copies of bootstrap estimates }
  \item{method}{ method to use during optimization, see details }
  \item{control}{ a list of control parameters for the optimization routine, see details }
}
\details{
  
  The function estimates the parameters \code{mu1}, \code{sig1}, \code{mu2}, \code{sig2} and \code{rho}. 
  In cases of incomplete sampling the estimates of \code{mu1} and \code{mu2} will be confounded with the sampling
  intensities (see \code{\link{rbipoilog}}). Assuming sampling intensities \eqn{\nu_1}{\nu1} and \eqn{\nu_2}{\nu2}, 
  the estimates of the means are \eqn{\code{mu1}+\ln \nu_1}{\code{mu1}+ \ln nu1} and \eqn{\code{mu2}+\ln\nu_2}{\code{mu2}+ \ln nu2}. Parameters
  \code{sig1}, \code{sig2} and \code{rho} can be estimated without any knowledge of sampling intensities.
  The parameters must be given starting values for the optimization procedure (default starting values are 
  used if starting values are not specified in the function call). The start values of \code{sig1} and
  \code{sig2} must be larger than 0 and at most 1000, the range supported by \code{\link{dbipoilog}}, and
  the search stays in that range: beyond it the likelihood is treated as not evaluable, and the optimizer
  moves away. The start value of \code{rho} must lie strictly between -1 and 1.
  
  A zero-truncated distribution (see \code{\link{dbipoilog}}) is assumed by default (\code{zTrunc = TRUE}).
  In cases where the number of zeros is known the \code{zTrunc} argument should be set to \code{FALSE}.
  
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
    \item the data have no maximum likelihood estimate: \code{n1 + n2} is 1 for every
      species with \code{zTrunc = TRUE}, or both counts of every species are 0 with
      \code{zTrunc = FALSE}. The likelihood then increases towards a limit as \code{mu1}
      and \code{mu2} go to \eqn{-\infty}{-Inf};
    \item the likelihood cannot be evaluated at the start values;
    \item the optimization does not converge (for example when \code{maxit} is reached);
    \item the fit does not end at a maximum of the likelihood. After the optimization has
      converged, the Hessian of the negative log-likelihood in (\code{mu1}, \code{mu2},
      \code{sig1}, \code{sig2}, \code{rho}) is computed from its analytic gradient. The end
      point counts as a maximum if every eigenvalue of the Hessian exceeds
      \eqn{10^{-6}}{1e-6} (a smaller one means a standard error above 1000 in some direction,
      so the data do not determine the parameters) and if a Newton step would increase the
      log-likelihood by less than 0.05. Likely reasons for a failure are start values from which
      the optimizer stopped short of the maximum or at a point that is not one, data for which no
      maximum exists, an estimate that approaches \code{sig1 = 1000} or \code{sig2 = 1000},
      and the method \code{"SANN"} or loose \code{control} settings (a large \code{reltol} or a
      small \code{maxit}), which can stop short of the maximum. Other start values, a larger
      \code{maxit} or a smaller \code{reltol} may help. When \code{sig1} or \code{sig2} is
      close to 0, \code{rho} has little effect on the likelihood; if the data then do not
      determine \code{rho}, the fit stops with this error;
    \item the numerical integration fails at or near the estimate.
  }
  The counts must be non-negative integers; with \code{zTrunc = TRUE}, no species may have
  zero counts in both samples.
  
  The approximate fraction of species revealed by each sample is estimated as 1 minus the zero term of the 
  univariate Poisson lognormal distribution: \eqn{1-q(0;\code{mu1},\code{sig1})} and \eqn{1-q(0;\code{mu2},\code{sig2})}.

  Parametric bootstrapping could be time consuming for large data sets. If argument \code{file} is specified, e.g. 
  \code{file =}\sQuote{\code{C:\\\\myboots.txt}}, the matrix with bootstrap estimates are copied into a tab-separated text-file 
  providing extra backup. Bootstrapping is done by simulating new sets of observations conditioned on the observed number of
  species (see \code{\link{rbipoilog}}), and refitting the model to each, starting from the estimates. A
  refit that does not converge or does not end at a maximum of the likelihood (including a simulated
  data set that has no maximum) is discarded and replaced
  by a new simulated data set, so the bootstrap distribution is conditioned on simulated data sets for
  which the maximum likelihood estimate exists. Any other error in a refit stops the bootstrap. The
  bootstrap stops with an error if \code{nboot} successful refits are not obtained from
  \code{10*nboot + 50} simulated data sets. A simulated data set with a count above
  \code{.Machine$integer.max}, which cannot be fitted, is replaced by a new one as well, which
  conditions the bootstrap on data sets without such counts; the number of data sets replaced for this
  reason is reported in a message.
  
}
\value{
  \item{par }{Maximum likelihood estimates of the parameters}
  \item{p }{Approximate fraction of species revealed by the samples for sample 1 and 2 respectively}
  \item{logLval }{Log likelihood of the data given the estimated parameters}
  \item{gof }{Goodness of fit measure: the proportion of bootstrap log-likelihoods that are smaller than
              logLval. (\code{gof} < 0.05) or (\code{gof} > 0.95) indicates lack of fit. NULL if \code{nboot = 0}}
  \item{boot }{A data frame with the bootstrap replicates of the parameters (columns \code{mu1},
              \code{mu2}, \code{sig1}, \code{sig2} and \code{rho}) and their log-likelihoods
              (\code{logLval}). NULL if \code{nboot = 0}}
  
}
\references{ 
Engen, S., R. Lande, T. Walla and P. J. DeVries. 2002. Analyzing spatial structure of communities using the two-dimensional 
            Poisson lognormal species abundance model. American Naturalist \bold{160}, 60-73.}
\author{ Vidar Grotan \email{vidar.grotan@ntnu.no}, Steinar Engen   }

\seealso{ \code{\link{optim}}, \code{\link{dbipoilog}}, \code{\link{rbipoilog}}  }
\examples{
\dontshow{set.seed(200)}
## simulate observations
xy  <- rbipoilog(S=30,mu1=1,mu2=1,sig1=2,sig2=2,rho=0.5)

## obtain estimates of parameters
est <- bipoilogMLE(xy)

## change start values and request tracing information 
## from optimization procedure
\dontrun{est <- bipoilogMLE(xy,startVals=c(2,2,4,4,0.3),
       control=list(maxit=1000,trace=1, REPORT=1))}

## effect of sampling intensity 
xy  <- rbipoilog(S=100,mu1=1,mu2=1,sig1=2,sig2=2,rho=0.5,nu1=0.5,nu2=0.5)
\dontrun{est <- bipoilogMLE(xy)}
## the expected estimates of mu1 and mu2 are now 1+log(0.5) = 0.31 (approximately)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
