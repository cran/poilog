   # ---------------------------------------------------------------- input checks
   # Shared by all exported functions, so that invalid input fails with the same
   # message everywhere, and always before anything reaches the C code.

   check_counts <- function(x, name, allow_empty = FALSE){
      if (!is.numeric(x)) stop(paste(name,'must be numeric'))
      if (length(x)==0 && !allow_empty) stop(paste(name,'is empty'))
      if (any(!is.finite(x))) stop(paste('one or several values of',name,'are NA or not finite'))
      if (any(x<0)) stop(paste('one or several values of',name,'are negative'))
      if (any(x!=round(x))) stop(paste('all values of',name,'must be integers'))
      if (any(x>.Machine$integer.max))
        stop(paste('values of',name,'larger than .Machine$integer.max are not supported'))
   }

   # one finite number
   check_param <- function(x, name){
      # a single NA (logical NA is not numeric) is a missing value, not a vector
      if (length(x)==1 && is.na(x)) stop(paste(name,'must be finite'))
      if (!is.numeric(x) || length(x)!=1)
        stop(paste(name,'must be a single number (vectorization of parameters is not implemented)'))
      if (!is.finite(x)) stop(paste(name,'must be finite'))
   }

   check_flag <- function(x, name){
      if (!is.logical(x) || length(x)!=1 || is.na(x)) stop(paste(name,'should be TRUE or FALSE'))
   }

   check_size <- function(S){
      check_param(S, 'S')
      if (S<1) stop('S is not positive')
      if (S!=round(S)) stop('S is not an integer')
   }

   # The supported range of sig (and sig1, sig2) in the probabilities and the fits.
   # The integration in the C code was validated against an independent reference up
   # to this value, on grids of counts up to .Machine$integer.max; it fails from about
   # sig = 1e4. sig = 1000 spreads the abundances over a factor of exp(1000), far beyond
   # any data. The random generators do not use the C code; they need only that every
   # Poisson intensity they draw is finite (poisson_draws).
   SIG_MAX <- 1000
   check_sig_max <- function(x, name){
      if (x>SIG_MAX) stop(paste(name,'larger than',SIG_MAX,'is not supported'))
   }

   # sqrt(1 - rho^2) for the C code, computed as sqrt((1 - rho) (1 + rho)): both factors
   # are exact for rho in [-1, 1], so the result has full relative accuracy also when
   # rho is close to +-1 (1 - rho^2 itself loses digits there)
   rho_c <- function(rho) sqrt((1-rho)*(1+rho))

   # Conditional simulation (condS = TRUE) draws batches until S positive counts are
   # found. The first batch has 2 S draws. Each further batch has 2 / p times n draws,
   # where p is the share of positive counts in the previous batch (1000 times n if it
   # had none) and n is S minus the number of positive counts in the previous batch
   # alone, not in all batches so far, so it can draw more than needed. In total at most
   # MAX_DRAWS(S) values are drawn; a request that needs more stops with an error
   # before the batch that would exceed it is drawn. This bounds the memory a batch can
   # take: with rare positive counts one batch can otherwise hold billions of draws.
   # The S term keeps every request with P(N > 0) above about 2% working for any S;
   # 2e7 draws take about 0.5 GB (univariate) or 1.3 GB (bivariate) at the peak.
   MAX_DRAWS <- function(S) max(2e7, 100*S)

   # Poisson counts with intensities lam. An intensity above the largest double
   # (about 1.8e308, a log intensity above 709.78) is Inf, for which rpois returns NA
   # with only a warning; that NA would be returned as a count, and with condS = TRUE
   # counted as a positive one. Checking the intensities before the draws uses no random numbers, so seeded
   # draws are unchanged wherever they are finite.
   poisson_draws <- function(nr, lam){
      if (any(lam==Inf))
        stop(paste('a Poisson intensity exp(mu + log(nu) + sig u) of the simulation exceeds',
                   'the largest double (log intensity above 709.78): mu, sig or nu are too large'))
      rpois(nr, lam)
   }


   # ---------------------------------------------------------------- random deviates

    rbipoilog <- function(S,mu1,mu2,sig1,sig2,rho,nu1=1,nu2=1,condS=FALSE,keep0=FALSE){

       sim <- function(nr){
         lamx <- rnorm(nr)
         lamy <- rho*lamx+sqrt(1-rho^2)*rnorm(nr)
         x <- poisson_draws(nr,exp(sig1*lamx+mu1+log(nu1)))
         y <- poisson_draws(nr,exp(sig2*lamy+mu2+log(nu2)))
         sel <- as.double(x)+y   # in double: integer x + y overflows near 2^31
         xy <- cbind(x,y)
         if (keep0) xy <- cbind(x,y) else xy <- cbind(x[sel>0],y[sel>0])
         return(xy)
       }

       check_size(S)
       for (p in c('mu1','mu2','sig1','sig2','rho','nu1','nu2')) check_param(get(p), p)
       if (sig1<0) stop('sig1 is not positive')
       if (sig2<0) stop('sig2 is not positive')
       if (nu1<0)    stop('nu1 is not positive')
       if (nu2<0)    stop('nu2 is not positive')
       if (-1>rho | rho>1) stop('rho must be in the range -1 to 1')
       check_flag(condS, 'condS')
       check_flag(keep0, 'keep0')
       if (condS && !keep0 && nu1==0 && nu2==0) stop('with nu1 = nu2 = 0 all counts are 0, so condS = TRUE with keep0 = FALSE is impossible')

       if (condS){
         simMat <- matrix(NA,0,2)
         fac <- 2
         nr  <- S
         drawn <- 0
         while (nrow(simMat)<S){
           drawn <- drawn+nr*fac
           if (drawn>MAX_DRAWS(S))
             stop(paste('P(N1 + N2 > 0) is too small to simulate S positive counts (condS = TRUE):',
                        'more than', MAX_DRAWS(S), 'draws would be needed'))
           simvals <- sim(nr*fac)
           simMat  <- rbind(simMat,simvals)
           fac     <- (1/(nrow(simvals)/(nr*fac)))*2
           fac     <- ifelse(is.finite(fac),fac,1000)
           nr      <- S-nrow(simvals)
         }
         simMat <- simMat[1:S,,drop=FALSE]
       }

       else simMat <- sim(S)
       return(simMat)
    }


    rpoilog <- function(S,mu,sig,nu=1,condS=FALSE,keep0=FALSE){

       sim <- function(nr){
         lamx <- rnorm(nr)
         x <- poisson_draws(nr,exp(sig*lamx+mu+log(nu)))
         if (!keep0) x <- x[x>0]
         return(x)
       }

       check_size(S)
       check_param(mu, 'mu')
       check_param(sig, 'sig')
       check_param(nu, 'nu')
       if (sig<0) stop('sig is not positive')
       if (nu<0)    stop('nu is not positive')
       check_flag(condS, 'condS')
       check_flag(keep0, 'keep0')
       if (condS && !keep0 && nu==0) stop('with nu = 0 all counts are 0, so condS = TRUE with keep0 = FALSE is impossible')

       if (condS){
         simVec <- vector('numeric',0)
         fac <- 2
         nr  <- S
         drawn <- 0
         while (length(simVec)<S){
           drawn <- drawn+nr*fac
           if (drawn>MAX_DRAWS(S))
             stop(paste('P(N > 0) is too small to simulate S positive counts (condS = TRUE):',
                        'more than', MAX_DRAWS(S), 'draws would be needed'))
           simvals <- sim(nr*fac)
           simVec  <- c(simVec,simvals)
           fac     <- (1/(length(simvals)/(nr*fac)))*2
           fac     <- ifelse(is.finite(fac),fac,1000)
           nr      <- S-length(simvals)
         }
         simVec <- simVec[1:S]
       }

       else simVec <- sim(S)
       return(simVec)
    }


   # ---------------------------------------------------------------- probabilities

   dbipoilog <- function(n1,n2,mu1,mu2,sig1,sig2,rho,log=FALSE){

      check_counts(n1, 'n1', allow_empty = TRUE)
      check_counts(n2, 'n2', allow_empty = TRUE)
      if (length(n1)!=length(n2)) stop('n1 and n2 must be of equal length')
      for (p in c('mu1','mu2','sig1','sig2','rho')) check_param(get(p), p)
      if (-1>rho | rho>1) stop('rho must be in the range -1 to 1')
      if (sig1<=0 | sig2<=0) stop('sig1 and/or sig2 is not larger than 0')
      check_sig_max(sig1, 'sig1')
      check_sig_max(sig2, 'sig2')
      check_flag(log, 'log')

     lv <- .C('poilog2',n1=as.integer(n1),n2=as.integer(n2),mu1=as.double(mu1),mu2=as.double(mu2),
        sig1=as.double(sig1^2),sig2=as.double(sig2^2),rho=as.double(rho),c=as.double(rho_c(rho)),
        nrN=as.integer(length(n1)),val=double(length(n1)), PACKAGE = "poilog")$val
     if (log) lv else exp(lv)
    }

   dpoilog <- function(n,mu,sig,log=FALSE){
     check_counts(n, 'n', allow_empty = TRUE)
     check_param(mu, 'mu')
     check_param(sig, 'sig')
     if (sig<=0) stop('sig is not larger than 0')
     check_sig_max(sig, 'sig')
     check_flag(log, 'log')
     lv <- .C("poilog1", n = as.integer(n), mu = as.double(mu),
              sig2 = as.double(sig^2), nrN = as.integer(length(n)),
              val = double(length(n)), PACKAGE = "poilog")$val
     if (log) lv else exp(lv)
   }

   # log P(N > 0) and log P(N1 + N2 > 0), the zero-truncation terms. The C code computes
   # them so that they stay accurate when P(N = 0) is close to 1 (see its header).
   # lpoilogpos is used by the MLEs for p; lbipoilogpos is the value-only
   # counterpart of the bivariate truncation term, kept for testing.
   lpoilogpos <- function(mu, sig){
     check_sig_max(sig, 'sig')
     .C("poilogpos1", mu = as.double(mu), sig2 = as.double(sig^2), val = double(1),
        PACKAGE = "poilog")$val
   }
   lbipoilogpos <- function(mu1, mu2, sig1, sig2, rho){
     check_sig_max(sig1, 'sig1')
     check_sig_max(sig2, 'sig2')
     .C("poilogpos2", mu1 = as.double(mu1), mu2 = as.double(mu2), sig1 = as.double(sig1^2),
        sig2 = as.double(sig2^2), rho = as.double(rho), c = as.double(rho_c(rho)), val = double(1),
        PACKAGE = "poilog")$val
   }


   # ---------------------------------------------------------------- likelihoods

   # optim works with z = (mu, log sig) or (mu1, mu2, log sig1, log sig2, w), with
   # w = log((1 + rho) / (1 - rho)), i.e. rho = tanh(w / 2). sqrt(1 - rho^2) is passed to
   # the C code as 1 / cosh(w / 2), computed from w: from the rounded rho it would lose
   # its digits, and round to 0, when rho is close to +-1.
   rhoTrans    <- function(rho) log((1+rho)/(1-rho))
   invRhoTrans <- function(w) tanh(w/2)
   rhoCTrans   <- function(w) 1/cosh(w/2)

   # Negative (zero-truncated) log-likelihood and its analytic gradient in z. The
   # kernel returns the gradient of log q from the same quadrature nodes as q itself.
   # The objective is defined for every z: the true value where every sig is at most
   # SIG_MAX, and +Inf beyond, which optim treats as a point it cannot evaluate and backs
   # off from. (optim asks for the gradient only at points with a finite value; there
   # is none beyond SIG_MAX, so NA is returned.) For very small sig, sig^2 = exp(2 z)
   # underflows to 0 and the C code returns the Poisson limit, which is the correct value
   # there. Inside the range the C code can fail only at extreme combinations, such as
   # a large mu with a very small sig; nll1 and nll2 pass such an error on, and fit_ml
   # decides what it means during the search.
   nll1 <- function(z, un, nr, zTrunc) {
     sig <- exp(z[2]); k <- length(un)
     if (sig>SIG_MAX) return(list(value = Inf, gradient = rep(NA_real_, 2)))
     r <- .C("poilog1g", n = as.integer(un), mu = as.double(z[1]), sig2 = as.double(sig^2),
             nrN = as.integer(k), val = double(k), gr = double(2 * k), PACKAGE = "poilog")
     v <- sum(r$val * nr); g <- colSums(matrix(r$gr, k, 2) * nr)
     if (zTrunc) {
       b <- .C("poilogpos1g", mu = as.double(z[1]), sig2 = as.double(sig^2), val = double(1),
               gr = double(2), PACKAGE = "poilog")
       v <- v - sum(nr) * b$val; g <- g - sum(nr) * b$gr
     }
     list(value = -v, gradient = -g * c(1, sig))
   }
   nll2 <- function(z, un, nr, zTrunc) {
     s <- exp(z[3:4]); k <- nrow(un)
     if (any(s>SIG_MAX)) return(list(value = Inf, gradient = rep(NA_real_, 5)))
     rho <- invRhoTrans(z[5]); cr <- rhoCTrans(z[5])
     r <- .C("poilog2g", n1 = as.integer(un[, 1]), n2 = as.integer(un[, 2]), mu1 = as.double(z[1]),
             mu2 = as.double(z[2]), sig1 = as.double(s[1]^2), sig2 = as.double(s[2]^2),
             rho = as.double(rho), c = as.double(cr), nrN = as.integer(k), val = double(k),
             gr = double(5 * k), PACKAGE = "poilog")
     v <- sum(r$val * nr); g <- colSums(matrix(r$gr, k, 5) * nr)
     if (zTrunc) {
       b <- .C("poilogpos2g", mu1 = as.double(z[1]), mu2 = as.double(z[2]), sig1 = as.double(s[1]^2),
               sig2 = as.double(s[2]^2), rho = as.double(rho), c = as.double(cr), val = double(1),
               gr = double(5), PACKAGE = "poilog")
       v <- v - sum(nr) * b$val; g <- g - sum(nr) * b$gr
     }
     # the C code returns c d/drho (c = sqrt(1 - rho^2)), and d rho / d w = c^2 / 2, so
     # d/dw = (c / 2) times that; finite also where c rounds to 0
     jac <- c(1, 1, s, cr / 2)
     list(value = -v, gradient = -g * jac)
   }
   # fn and gr for optim sharing one kernel call per point (optim asks for the value
   # and the gradient at the same point separately)
   nll_obj <- function(f) {
     key <- NULL; res <- NULL
     at <- function(z, ...) {
       k <- list(z, ...)
       if (!identical(k, key)) { res <<- f(z, ...); key <<- k }
       res
     }
     list(fn = function(z, ...) at(z, ...)$value, gr = function(z, ...) at(z, ...)$gradient)
   }


   # ---------------------------------------------------------------- fitting

   # optim methods this package supports. The analytic gradient is passed only to
   # the methods that use one: for 'SANN', optim would take it as the generator of
   # candidate points. 'Brent' is not supported: it is for one parameter only, and the
   # models have two or five. 'L-BFGS-B' is the bounded
   # method, and the fitting functions take no bounds; without them it stops with an
   # error at any point where the likelihood cannot be evaluated, where the other
   # methods back off, so it is refused with its own message.
   OPTIM_METHODS   <- c('BFGS','CG','Nelder-Mead','SANN')
   GRADIENT_METHODS <- c('BFGS','CG')

   check_fit_args <- function(startVals, k, nboot, zTrunc, method){
      if (!is.numeric(startVals)) stop('startVals must be numeric')
      if (length(startVals)!=k) stop(paste('length of startVals is not',k))
      if (any(!is.finite(startVals))) stop('all startVals must be finite')
      check_param(nboot, 'nboot')
      if (nboot<0 || nboot!=round(nboot)) stop('nboot must be a non-negative integer')
      check_flag(zTrunc, 'zTrunc')
      if (identical(method, 'L-BFGS-B'))
        stop(paste("method 'L-BFGS-B' is not supported: it is the method for bounded",
                   "parameters, and this function takes no bounds. Use one of", paste(OPTIM_METHODS, collapse = ', ')))
      if (!is.character(method) || length(method)!=1 || !(method %in% OPTIM_METHODS))
        stop(paste('method must be one of', paste(OPTIM_METHODS, collapse = ', ')))
   }

   # The two models in the form the fitting code needs. theta: the natural parameters
   # (mu, sig) or (mu1, mu2, sig1, sig2, rho) at z; z: back from theta, through log|sig|
   # (the likelihood depends on sig only through sig^2, so it is an even, smooth function
   # of sig, and a step may cross sig = 0); dz: d z / d theta, from z where rho is
   # involved (d w / d rho = 2 / (1 - rho^2) = 2 cosh(w / 2)^2); ndeps: the steps of the
   # numerical Hessian in theta (see at_maximum); no_max: for data whose likelihood has
   # no maximum, the reason, else NULL (see fit_ml).
   #  - Without zero truncation, data in which every count is 0 (in the bivariate model,
   #    both counts of every species): the likelihood increases towards 1 as mu (mu1
   #    and mu2) goes to -Inf.
   #  - With zero truncation, data in which every count is 1 (in the bivariate model,
   #    n1 + n2 = 1 for every species): the conditional probability of a count of 1
   #    given a positive count increases towards 1 as mu goes to -Inf (with mu1 - mu2
   #    matching the split between the two samples).
   # In both cases the likelihood comes as close to its limit as one likes but never
   # reaches it, so any finite point is not a maximum.
   MODEL1 <- list(
      theta = function(z) c(z[1], exp(z[2])),
      z     = function(th) c(th[1], log(abs(th[2]))),
      dz    = function(th, z) c(1, 1/th[2]),
      ndeps = function(th) 1e-3*pmax(abs(th), 0.1),
      no_max = function(un, zTrunc){
         if (!zTrunc && all(un==0)) return('every count is 0')
         if (zTrunc && all(un==1)) return('every count is 1')
         NULL
      })
   MODEL2 <- list(
      theta = function(z) c(z[1:2], exp(z[3:4]), invRhoTrans(z[5])),
      z     = function(th) c(th[1:2], log(abs(th[3:4])), rhoTrans(th[5])),
      dz    = function(th, z) c(1, 1, 1/th[3:4], 2*cosh(z[5]/2)^2),
      ndeps = function(th) c(1e-3*pmax(abs(th[1:4]), 0.1), 1e-3*(1-abs(th[5]))),
      no_max = function(un, zTrunc){
         tot <- as.double(un[,1])+as.double(un[,2])
         if (!zTrunc && all(tot==0)) return('both counts of every species are 0')
         if (zTrunc && all(tot==1)) return('n1 + n2 is 1 for every species')
         NULL
      })

   # Does the fit end at a maximum of the likelihood? optim can report convergence at a
   # point that is not one: a saddle, a ridge along which the likelihood still rises, a
   # drift towards a boundary, or data for which no maximum exists. The test is made in
   # the natural parameters theta, where a real maximum is well conditioned (also on the
   # boundary sig = 0, the Poisson limit, where in log sig the likelihood flattens). H is
   # the Hessian of -logL from differences of the analytic gradient (optimHess), with
   # steps of 1e-3 of each parameter (at least 1e-4), and for rho 1e-3 of its distance
   # from +-1 so that no step crosses +-1. Such steps are large against the numerical
   # noise of the gradient (below 1e-7 relative), so that its differences resolve the
   # curvature, and small against the scale on which the curvature itself changes near a
   # maximum. Two conditions:
   #  1. every eigenvalue of H exceeds EIGEN_MIN = 1e-6. A smaller one means a standard
   #     error above 1000 in some direction of theta: the data do not determine the
   #     parameters there, so there is no identifiable maximum; a negative one means the
   #     point is not a maximum at all.
   #  2. the Newton decrement g' H^-1 g / 2 (g the gradient), the log-likelihood that a
   #     Newton step would still gain, is below DECREMENT_MAX = 0.05: a gain of 0.05 (a
   #     likelihood-ratio statistic of 0.1) is negligible against sampling variation, so
   #     a fit that could gain more stopped short of the maximum.
   # Both thresholds are statistical, not numerical: at real maxima the smallest
   # eigenvalue is far above 1e-6 and the decrement far below 0.05, and the numerical
   # error of H (from the gradient) is far below both. If a step of H leaves the supported
   # range of sig, or H is not finite, the test fails.
   EIGEN_MIN <- 1e-6
   DECREMENT_MAX <- 0.05
   at_maximum <- function(z, ob, tab, model){
      th <- model$theta(z)
      fn <- function(t) ob$fn(model$z(t), nr = tab$nr, un = tab$un)
      gr <- function(t) { zt <- model$z(t); ob$gr(zt, nr = tab$nr, un = tab$un) * model$dz(t, zt) }
      H <- optimHess(th, fn, gr, control = list(ndeps = model$ndeps(th)))
      g <- gr(th)
      if (!all(is.finite(H)) || !all(is.finite(g))) return(FALSE)
      # the decrement from the eigen decomposition H = V diag(e) V', g' H^-1 g =
      # sum((V' g)^2 / e), which exists whenever condition 1 holds (no solve() that
      # could fail on a badly conditioned H)
      E <- eigen(H, symmetric = TRUE)
      if (min(E$values)<=EIGEN_MIN) return(FALSE)
      0.5*sum(crossprod(E$vectors, g)^2/E$values)<DECREMENT_MAX
   }

   # A fit that failed: optim did not converge, or the end point is not a maximum. It
   # has its own condition class, so that a bootstrap refit can reject exactly these
   # failures and let any other error through.
   fit_failure <- function(msg)
      structure(class = c('poilog_fit_failure', 'error', 'condition'),
                list(message = msg, call = NULL))

   # Fit from the start values z (in the optimiser's coordinates); returns optim's
   # result, or stops. The one routine used by poilogMLE, bipoilogMLE and their
   # bootstrap refits.
   # Trial points the C code cannot evaluate. The first step of BFGS (and CG) is the
   # negative gradient itself, and the gradient of the log-likelihood grows with the
   # number of species, so with a few hundred species the first trial point can lie at
   # absurd values, such as mu above 100 with sig below 1e-20, from the default start on
   # ordinary data. There the integration can fail (an error whose message starts with
   # 'poilog:'). During the search such a point gets the value +Inf, which optim treats
   # as a point it cannot evaluate: BFGS, CG and Nelder-Mead shorten the step or move
   # away and never accept it; SANN never moves to it from a point it can evaluate.
   # Only these errors are treated so; any other error stops the fit. The start values
   # are evaluated first without this rule, so that a failure there stops the fit with
   # the C code's own message. The end point is checked without the rule too: its value
   # must be finite, and at_maximum evaluates the likelihood and gradient around it
   # directly, so a failure of the C code at or near the estimate still stops the fit
   # with its error.
   # Data whose likelihood has no maximum (model$no_max) stop before the search, with a
   # message that names the cause. Without this, a fit on such data runs towards the
   # limit; with the default control the end-point test then stops it, but only with its
   # general message, and a fit that stops within about 0.1 of the limit (as with a
   # reltol of 1e-2 or 1e-3) passes that test, since there the Newton decrement is about
   # half the log-likelihood still to be gained. In a bootstrap refit this is a fit
   # failure like the others, so such a simulated data set is replaced.
   fit_ml <- function(z, ob, tab, control, method, model, zTrunc){
      why <- model$no_max(tab$un, zTrunc)
      if (!is.null(why))
        stop(fit_failure(paste0('the likelihood has no maximum for these data (', why,
          '): it increases towards a limit as mu goes to -Inf, so there is no estimate')))
      gr <- if (method %in% GRADIENT_METHODS) ob$gr else NULL
      ob$fn(z, nr=tab$nr, un=tab$un)
      fn <- function(z, ...) tryCatch(ob$fn(z, ...), error = function(e)
        if (startsWith(conditionMessage(e), 'poilog:')) Inf else stop(e))
      fit <- optim(z,fn,gr,nr=tab$nr,un=tab$un,control=control,method=method)
      if (fit$convergence!=0){
        if (fit$convergence==1) stop(fit_failure('the iteration limit has been reached!   try different startVals or increase maxit'))
        if (fit$convergence==10) stop(fit_failure('degeneracy of the Nelder Mead simplex ....'))
        else stop(fit_failure(paste('unknown error in optimization', fit$message)))
      }
      # with the start values finite, optim ends only at a point it evaluated as finite,
      # so this would be a defect
      if (!is.finite(fit$value) || any(!is.finite(fit$par)))
        stop('the fit ended at a point where the likelihood is not finite')
      if (!at_maximum(fit$par, ob, tab, model))
        stop(fit_failure(paste('the fit did not end at a maximum of the likelihood.',
          'Likely reasons: start values from which the optimiser stopped short of the maximum,',
          'or at a point that is not one; data for which no maximum exists (for example',
          'bivariate data that do not determine rho); or an estimate that approaches the',
          'edge of the supported range',
          paste0('of sig (', SIG_MAX, ').'), 'Try other startVals, a larger maxit, or a smaller reltol',
          'in control.')))
      fit
   }

   # Parametric bootstrap shared by both MLEs. simulate() draws a data set from the
   # fitted model; refit(sim) returns c(estimates, logLval), or NULL when the refit
   # failed. A failed refit is replaced by a new draw, up to boot_draws(nboot) draws in
   # total, so the bootstrap distribution is conditioned on data sets that can be fitted.
   # The cap, 10 draws per requested refit plus 50 (so that a small nboot still allows
   # some failures), stops a bootstrap in which most refits fail instead of letting it
   # run on; refits that fail are typically far fewer (none to a quarter of the draws).
   # The C code takes counts as C integers, so a data set with a count above
   # .Machine$integer.max cannot be refitted either; it is redrawn as well, which
   # conditions the bootstrap on data sets without such counts. How often this happens
   # depends on the fitted model (at mu = 0, sig = 6 and 300 species, at least one such
   # count occurs in about 5% of the data sets), so the number redrawn is reported.
   boot_draws <- function(nboot) 10*nboot+50
   run_boot <- function(nboot, simulate, refit, cols, report){
      bMat <- matrix(NA,nboot,length(cols))
      colnames(bMat) <- cols
      count <- 0; draws <- 0; redrawn <- 0
      while (count<nboot){
        draws <- draws+1
        if (draws>boot_draws(nboot))
          stop(paste('the bootstrap stopped after', boot_draws(nboot), 'simulated data sets with only',
                     count, 'successful refits out of', nboot))
        sim <- simulate()
        if (any(sim>.Machine$integer.max)) { redrawn <- redrawn+1; next }
        row <- refit(sim)
        if (!is.null(row)){
          count <- count+1
          bMat[count,] <- row
          report(count, bMat)
        }
      }
      if (redrawn>0)
        message(paste(redrawn, 'simulated data set(s) had a count above .Machine$integer.max',
                      'and were replaced by new draws'))
      bMat
   }

   # a bootstrap refit from the fitted values zhat: optim's result, or NULL if the refit
   # failed (did not converge, or did not end at a maximum)
   boot_refit <- function(zhat, ob, sim, control, method, model, zTrunc)
      tryCatch(fit_ml(zhat, ob, tabulate_counts(sim), control, method, model, zTrunc),
               poilog_fit_failure = function(e) NULL)

   # goodness of fit: the share of bootstrap log-likelihoods below the observed one
   boot_gof <- function(logLval, bootLogL) mean(bootLogL<logLval)

   # number of occurrences of each unique count, or each unique row of a 2-column matrix
   tabulate_counts <- function(x){
      if (is.matrix(x)){
        un <- unique(x)
        nr <- vapply(seq_len(nrow(un)), function(i) sum(x[,1]==un[i,1] & x[,2]==un[i,2]), 0)
      } else {
        un <- unique(x)
        nr <- vapply(un, function(u) sum(x==u), 0)
      }
      list(un = un, nr = nr)
   }

   poilogMLE <- function(n,startVals=c(mu=1,sig=2),nboot=0,zTrunc=TRUE,
                         method='BFGS',control=list(maxit=1000)){

      if (is.matrix(n) | (is.data.frame(n))) stop(paste('n has',ncol(n),'colums, supply a vector or use function bipoilogMLE',sep=' '))
      check_counts(n, 'n')
      check_fit_args(startVals, 2, nboot, zTrunc, method)
      if (zTrunc && any(n==0))
        stop('n contains zeros, which a zero-truncated fit cannot use: remove them or use zTrunc = FALSE')
      if (startVals[2]<=0) stop('start value of sig is not larger than 0')
      check_sig_max(startVals[2], 'start value of sig')
      startVals[2] <- log(startVals[2])
      tab <- tabulate_counts(n)

      ob <- nll_obj(function(z, nr, un) nll1(z, un, nr, zTrunc))
      fit <- fit_ml(startVals, ob, tab, control, method, MODEL1, zTrunc)

      fit$par <- c(as.numeric(fit$par),exp(lpoilogpos(fit$par[1],exp(fit$par[2]))))
      est <- list('par'=c('mu'=fit$par[1],'sig'=exp(fit$par[2])),'p'=fit$par[3],'logLval'=-fit$value,'gof'=NULL,boot=NULL)

      if (nboot>0){
        cat(paste('estimates: mu: ',round(fit$par[1],3),' sig ',round(exp(fit$par[2]),3),sep=''),'\n')
        cat('********     bootstrapping    ********\n')
        bStartVals <- fit$par[1:2]
        # simulations are conditional on the number of species in the observed data set
        simulate <- function() rpoilog(length(n),fit$par[1],exp(fit$par[2]),condS=TRUE,keep0=!zTrunc)
        refit <- function(sim){
          bfit <- boot_refit(bStartVals, ob, sim, control, method, MODEL1, zTrunc)
          if (is.null(bfit)) return(NULL)
          c(bfit$par[1],exp(bfit$par[2]),-bfit$value)
        }
        report <- function(count, bMat) if (count%%100==0) cat('   boot',count,'of',nboot,'\n')
        bMat <- run_boot(nboot, simulate, refit, c('mu','sig','logLval'), report)
        est$boot <- data.frame(bMat)
        est$gof  <- boot_gof(est$logLval, bMat[,3])
      }
      return(est)
   }



   bipoilogMLE <- function(n1,n2=NULL,startVals=c(mu1=1,mu2=1,sig1=2,sig2=2,rho=0.5),
                           nboot=0,zTrunc=TRUE,file=NULL,
                           method='BFGS',control=list(maxit=1000)){


      if (is.null(n2)){
        if (!is.matrix(n1) & !is.data.frame(n1)) {stop('n2 is missing and n1 is not a matrix or data frame') }
        if (ncol(n1)!=2) stop('n1 should have 2 columns when n2 is not given')
        n2 <- n1[,2]
        n1 <- n1[,1]
      }
      if (!is.null(n2) & (is.matrix(n1) | is.data.frame(n1))) stop('n1 is a matrix and n2 is specified')
      if (is.matrix(n2) | is.data.frame(n2)) stop('n2 must be a vector')
      if (length(n1)!=length(n2)) stop('n1 and n2 must be of equal length')
      check_counts(n1, 'n1')
      check_counts(n2, 'n2')
      check_fit_args(startVals, 5, nboot, zTrunc, method)
      # in double: integer n1 + n2 overflows near 2^31
      if (zTrunc && any(as.double(n1)+as.double(n2)==0))
        stop('some species have zero counts in both samples, which a zero-truncated fit cannot use: remove them or use zTrunc = FALSE')

      z <- startVals
      if (z[3]<=0) stop('start value of sig1 is not larger than 0')
      if (z[4]<=0) stop('start value of sig2 is not larger than 0')
      check_sig_max(z[3], 'start value of sig1')
      check_sig_max(z[4], 'start value of sig2')
      if (-1>=z[5] | z[5]>=1) stop('start value of rho must be in the range -1 to 1 (exclusive)')

      z[3:4] <- log(z[3:4])
      z[5]   <- rhoTrans(z[5])
      tab <- tabulate_counts(cbind(n1,n2))

      ob <- nll_obj(function(z, nr, un) nll2(z, un, nr, zTrunc))
      fit <- fit_ml(z, ob, tab, control, method, MODEL2, zTrunc)

      fit$par <- c(as.numeric(fit$par),exp(lpoilogpos(fit$par[1],exp(fit$par[3]))),exp(lpoilogpos(fit$par[2],exp(fit$par[4]))))
      est <- list('par'=c('mu1'=fit$par[1],'mu2'=fit$par[2],'sig1'=exp(fit$par[3]),
             'sig2'=exp(fit$par[4]),'rho'=invRhoTrans(fit$par[5])),'p'=c(fit$par[6],fit$par[7]),'logLval'=-fit$value,'gof'=NULL,boot=NULL)

      if (nboot>0){
        cat(paste('estimates: mu1:',round(fit$par[1],3), ' mu2:',round(fit$par[2],3),
                    ' sig1:',round(exp(fit$par[3]),3), ' sig2:',round(exp(fit$par[4]),3),
                    ' rho:',  round(invRhoTrans(fit$par[5]),3), sep=''),'\n')
        cat('********     bootstrapping    ********\n')
        bStartVals <- fit$par[1:5]
        # simulations are conditional on the number of species in the observed data set
        simulate <- function() rbipoilog(length(n1),fit$par[1],fit$par[2],exp(fit$par[3]),exp(fit$par[4]),
                                         invRhoTrans(fit$par[5]),condS=TRUE,keep0=!zTrunc)
        refit <- function(sim){
          bfit <- boot_refit(bStartVals, ob, sim, control, method, MODEL2, zTrunc)
          if (is.null(bfit)) return(NULL)
          c(bfit$par[1],bfit$par[2],exp(bfit$par[3]),exp(bfit$par[4]),invRhoTrans(bfit$par[5]),-bfit$value)
        }
        report <- function(count, bMat){
          cat('   boot',count,'of',nboot,': ',bMat[count,1:5],'\n')
          if (!is.null(file)) write.table(bMat[1:count,,drop=FALSE],file,sep='\t')
        }
        bMat <- run_boot(nboot, simulate, refit, c('mu1','mu2','sig1','sig2','rho','logLval'), report)
        est$boot <- data.frame(bMat)
        est$gof  <- boot_gof(est$logLval, bMat[,6])
      }

      return(est)
   }
