\name{Poisson lognormal}
\alias{dpoilog}
\alias{rpoilog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Poisson lognormal distribution }
\description{
  Density and random generation for the Poisson lognormal distribution with parameters \code{mu} and \code{sig}.
}
\usage{
dpoilog(n, mu, sig, log = FALSE)
rpoilog(S, mu, sig, nu=1, condS=FALSE, keep0=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ vector of observed individuals for each species }
  \item{S}{ number of species in the community }
  \item{mu}{ mean of the log abundances (of the normal distribution of ln \emph{y}, see
    details) }
  \item{sig}{ standard deviation of the log abundances; for \code{dpoilog} larger than 0
    and at most 1000 }
  \item{log}{ logical; if TRUE the log probability is returned }
  \item{nu}{ sampling intensity, defaults to 1 }
  \item{condS}{ logical; if TRUE random deviates are conditional on S}
  \item{keep0}{ logical; if TRUE species with count 0 are included in the random deviates }
  
}
\details{

  The following is written from the perspective of using the Poisson lognormal distribution 
  to describe community structure (the distribution of species when sampling individuals 
  from a community of several species). 
  
  Under the assumption of random sampling, the number of individuals sampled from a given 
  species with abundance \emph{y}, say \emph{N}, is Poisson distributed with mean \eqn{\code{nu}\,y}{\code{nu}\emph{y}} 
  where the parameter \code{nu} expresses the sampling intensity. If ln \emph{y} is normally distributed 
  with mean \code{mu} and standard deviation \code{sig} among species, then the vector of individuals sampled 
  from all \code{S} species then constitutes a sample from the Poisson lognormal distribution
  with parameters (\code{mu} + ln \code{nu}, \code{sig}), where \code{mu} and \code{sig} 
  are the mean and standard deviation of the log abundances. For \code{nu} = 1, this is the Poisson 
  lognormal distribution with parameters (\code{mu},\code{sig}) which may be written in the form
    
  \deqn{P(N=\code{n};\code{mu},\code{sig}) = q(\code{n};\code{mu},\code{sig}) = \int\limits_{-\infty}^{\infty} g_\code{n}(\code{mu},\code{sig},u)\phi(u)\;du,}{
        P(N=\code{n}; \code{mu},\code{sig}) = q(\code{n}; \code{mu},\code{sig}) = \int_-infty^infty g_\code{n}(\code{mu},\code{sig},u) phi(u) du,}
  
  where \eqn{\phi(u)} is the standard normal distribution and
    
  \deqn{g_\code{n}(\code{mu},\code{sig},u) = \frac{\exp(u\,\code{sig}\,\code{n} + \code{mu}\,\code{n} - \exp(u\,\code{sig} + \code{mu}))}{\code{n}!}}{
        g_n(\code{mu},\code{sig},u) = exp(u \code{sig} \code{n} + \code{mu} \code{n} - exp(u \code{sig} + \code{mu})) / \code{n}!}
       
   
  Since \code{S} is usually unknown, we only consider the observed number of individuals for the observed species.
  With a general sampling intensity \code{nu}, the distribution of the number of individuals then follows the
  zero-truncated Poisson lognormal distribution
    
   \deqn{\frac{q(\code{n};\code{mu}+\ln\code{nu},\code{sig})}{1-q(0;\code{mu}+\ln\code{nu},\code{sig})}}{
          q(\code{n}; \code{mu} + ln \code{nu},\code{sig})/(1 - q(0; \code{mu} + ln \code{nu},\code{sig}))}

  The integral is computed numerically on the log scale. The error of the log
  probability is about \eqn{10^{-8}}{1e-8} or less relative to its size, and about
  \eqn{10^{-16}}{1e-16} or less in absolute terms where the log probability is
  between -1 and 0 (a probability close to 1). Probabilities below about
  \eqn{10^{-308}}{1e-308}, the smallest normal number in double precision, lose
  precision, and those below about \eqn{5 \times 10^{-324}}{5e-324} are returned as 0;
  use \code{log = TRUE} to obtain their logarithm, which is accurate also there.

  Counts must be non-negative integers no larger than \code{.Machine$integer.max};
  \code{mu} and \code{sig} must be single finite numbers. \code{dpoilog} requires
  \code{sig > 0} and supports \code{sig} up to 1000, the range over which the numerical
  integration has been validated; larger values give an error. (\code{sig = 1000}
  spreads the abundances over a factor of \eqn{e^{1000}}{exp(1000)}, far beyond any
  data.) \code{rpoilog} does not use the integration and accepts any \code{sig >= 0}
  for which the Poisson intensities it draws stay below the largest double (a log
  intensity \code{mu} + ln \code{nu} + \code{sig} \emph{u} below 709.78); otherwise it
  stops with an error.
  With \code{condS = TRUE} and \code{keep0 = FALSE}, \code{rpoilog} draws until
  \code{S} positive counts are found, and stops with an error when this is impossible
  (\code{nu = 0}) or when positive counts are too rare to be found: when the
  simulation would draw more than \code{max(2e7, 100*S)} values. With
  \code{keep0 = TRUE} it returns \code{S} counts, zeros included (all zeros when
  \code{nu = 0}).
}
\value{
  \command{dpoilog} returns the probabilities \eqn{q(\code{n};\code{mu},\code{sig})},
  or their logarithms if \code{log = TRUE}\cr
  \command{rpoilog} returns random deviates
}
\references{ Engen, S., R. Lande, T. Walla & P. J. DeVries. 2002. Analyzing spatial structure of communities using the two-dimensional 
            Poisson lognormal species abundance model. American Naturalist 160: 60-73.
 }
\author{ Vidar Grotan \email{vidar.grotan@ntnu.no} and Steinar Engen   }

\seealso{\code{\link{poilogMLE}} for ML estimation}
\examples{

### plot density for given parameters 
barplot(dpoilog(n=0:20,mu=2,sig=1),names.arg=0:20)

### draw random deviates from a community of 50 species 
rpoilog(S=50,mu=2,sig=1)

### draw random deviates including zeros 
rpoilog(S=50,mu=2,sig=1,keep0=TRUE)

### draw random deviates with sampling intensity = 0.5 
rpoilog(S=50,mu=2,sig=1,nu=0.5)

### how many species are likely to be observed 
### (given S,mu,sig2 and nu)? 
hist(replicate(1000,length(rpoilog(S=30,mu=0,sig=3,nu=0.7))))

### how many individuals are likely to be observed
### (given S,mu,sig2 and nu)? 
hist(replicate(1000,sum(rpoilog(S=30,mu=0,sig=3,nu=0.7))))


}
\keyword{ distribution }
