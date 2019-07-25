
# (c) Pablo Vena, 2019

library(Rcpp)                    # C++
library(RcppArmadillo)           # C++
sourceCpp('huber.cpp')  # Huber without scale estimation



## This function is used internally by FPLMBsplines(), not meant
## for "public" use.

minimize <- function (yy, xx_coef, uu, spl_kn, freq, fLoss, norder, pars,
                      rob.control=lmrob.control(trace.level = 0,         # 0
                                                nResample   = 5000,      # 500 default
                                                tuning.psi  = 3.443689,      # para 85% eff usar 3.443689
                                                subsampling = 'simple',  #
                                                rel.tol     = 1e-5,      # 1e-7
                                                refine.tol  = 1e-5,      # 1e-7
                                                k.max       = 2e3,       # 200
                                                maxit.scale = 2e3,       # 200
                                                max.it      = 2e3)) {

  ## Initialize
  cv  <- tr <- aicM <- NA

  ## B-spline basis
  kns    <- seq(min(uu), max(uu), length = spl_kn - norder + 2)
  base   <- create.bspline.basis(rangeval = range(uu),
                                 norder = norder,
                                 breaks = kns)
  spl_uu <- getbasismatrix(uu, base)

  ## Design matrix

  X <- cbind(xx_coef, spl_uu)

  if( !(fLoss %in% c('ls','huang', 'lmrob') ) )
    stop("Invalid fLoss. Should be one of \'ls\' or \'huang\' or \'lmrob\' ")

  ## Minimization menu
  switch(fLoss,
         ls = {
           fit <- lm(yy ~ X - 1)
           cf  <- fit$coef
           vv  <- sum(fit$res^2)  # 1
           ss  <- sqrt(mean(fit$res^2))
         },
         huang = {
           init <- lm(yy ~ X - 1)$coef
           # X1=cbind(rep(1,length(yy)), X)
           fit  <- huber_estimates(X, yy, init, 1.345, 1e-8)
           cf   <- as.vector(fit$param)
           ss   <- 1
           vv   <- fit$value
           cv   <- 2
         },
         lmrob = {
           fit  <- lmrob(yy ~ X - 1, control = rob.control)
           if (fit$init.S$converged) {
             cf <- fit$coef
             ss <- fit$scale
             vv <- sum(Mpsi(fit$res / ss,
                            cc  = rob.control$tuning.psi,
                            psi = rob.control$psi, deriv = -1))
             cv <- fit$converged
           } else {
             stop('S-estimator did not converge.')
           }
         },
         stop('Invalid fLoss.')
  )

  ## Estimated parameters
  # intercept <- cf[1]
  # slope_par <- cf[2:(freq + 1)]
  # spl_par <- cf[-(1:(freq + 1))]
  spl_par   <- cf[-(1:freq)]
  slope_par <- cf[  1:freq ]

  return(list(spl = spl_par, slope = slope_par, #intercept = intercept,
              value = vv, scale = ss, conv = cv, traza = tr, matias = aicM))
}


FPLMBsplines_fit <- function (y, x, u, t, freq, spl,
                              norder, fLoss) {

  ## Integration step
  dt <- min(diff(t))  # width of grid
  xcenter <- x

  ## Covariance decomposition
  grilla_spl <- t #seq(0, 1, length = length(t))
  nodos_spl <- seq(min(grilla_spl), max(grilla_spl),
                   length = freq - norder + 2)
  base_spl <- create.bspline.basis(rangeval = range(t),
                                   norder = norder,
                                   breaks = nodos_spl)
  beta_spl <- getbasismatrix(grilla_spl, base_spl)
  cov_dec <- beta_spl[, 1:freq]

  ## Estimated Fourier coefficients (by row)
  xx_coef <- xcenter %*% cov_dec * dt

  ## Parameter estimation
  est <- minimize(y, xx_coef, u, spl, freq, fLoss,
                  norder, pars)

  est$slope_fun <- cov_dec %*% est$slope

  return (est)
}

goodness <- function(nn, scl, val, spl,freq, criterion){
  switch(criterion,
         hic     = log(scl^2 * val/ nn) + spl * log(nn) / (2 * nn) + 2 * freq / nn,
         hic2    = log(scl^2 * val/ nn) + freq * log(nn) / (2 * nn) + 2 * spl / nn,,
         hicGR2  = log(scl^2) + val / nn + freq * log(nn) / (2 * nn) + 2 * spl / nn,
         hicGR   = log(scl^2) + val / nn + spl * log(nn) / (2 * nn) + 2 * freq / nn,
         bic     = log(scl^2 * val/ nn) + (spl + freq) * log(nn) / (2 * nn),
         bic1     = log(scl^2 * val/ nn) + (spl + freq) * log(nn) / ( nn),
         bicGR   = log(scl^2) + val / nn  + (spl + freq) * log(nn) / (2 * nn),
         aicCL   = log(scl^2 * val/ nn) +  2 * (spl + freq) / nn,
         aicCLGR = log(scl^2) + val / nn  +  2 * (spl + freq) / nn)
}


## Model: y = <X, beta> + g(u) + e
## y: vector of n responses
## x: matrix of functional covariates (n x p_1)
## u: vector of n observations for g(u_i)
## t: vector of "times" where X was observed
## "x[i, j] = X_i( t_j )"
## freq:
## spl:
## norder:
## fLoss: loss function / method to use
##        valid options: 'ls', 'lad', 'lmrob', 'huang'


FPLMBsplines <- function (y, x, u, t, range_freq = (floor(max(n^(1/5), norder)):floor(2*(norder+n^(1/5)))),
                          range_spl= (floor(max(n^(1/5), norder)):floor(2*(norder+n^(1/5)))),
                          norder=4, fLoss='lmrob', criterion = 'bic1',
                          trace=FALSE) {

  opt <- spl_opt <- freq_opt <- fit_opt <- Inf

  n <- length(y)

  for(spl in range_spl) {
    for(freq in range_freq) {
      fit <- FPLMBsplines_fit(y, x, u, t, freq, spl, norder, fLoss)
      val <- fit$value
      scl <- fit$scale
      # nn  <- length(y)
      crt <- goodness(n, scl, val, spl, freq, criterion)
      if (crt < opt) {
        opt      <- crt
        spl_opt  <- spl
        freq_opt <- freq
        fit_opt  <- fit
      }
      if(trace) print(c('spl' = spl, 'freq' = freq, 'crit' = crt))
    }
  }

  kns    <- seq(min(u), max(u), length = spl_opt - norder + 2)
  base   <- create.bspline.basis(rangeval = range(u),
                                 norder = norder,
                                 breaks = kns)
  spl_uu <- getbasismatrix(u, base)
  fit_opt$eta_est <- spl_uu %*% fit_opt$spl
  dt <- min(diff(t))
  fit_opt$fitted <- as.vector( x %*% fit_opt$slope_fun * dt + fit_opt$eta_est )
  return(list(fit = fit_opt, spl = spl_opt, freq = freq_opt,
              u = u, t = t))
}

