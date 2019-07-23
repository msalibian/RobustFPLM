## ---------------------------------------------------------------------
## Functions for the analysis of the TECATOR data
## ---------------------------------------------------------------------

## Tuning for Tukey's bisquare is 3.443689
##
## The function "submuestra" is TECATOR-specific: it builds the
## response vector and design matrix for the specific varying coefficients model
## fat = < beta, absorp2> + water * eta(protein)

## fLoss is the type of estimator
## fLoss = 'ls' least squares
## fLoss = 'lmrob' MM-estimator with tuning constant cterho in the M-step

minimizar <- function (yy, xx_coef, uu, interac, spl_kn, freq, fLoss, norder, pars,cterho=3.443689, nresamp=5000) {
  
  ## Initialize
  cv  <- tr <- aicM <- NA
  
  ## B-spline basis
  kns    <- seq(min(uu), max(uu), length = spl_kn - norder + 2)
  base   <- create.bspline.basis(rangeval = range(uu),
                                 norder = norder,
                                 breaks = kns)
  spl_uu <- getbasismatrix(uu, base)
  
  ## Design matrix
  X <- cbind(xx_coef, spl_uu * interac)
  
  if( !(fLoss %in% c('ls', 'lmrob') ) )
    stop("Invalid fLoss. Should be one of \'ls\' or \'lmrob\' ")
  
  ## Minimization menu
  if (fLoss == 'ls') {
    fit <- lm(yy ~ X)
    cf  <- fit$coef
    vv  <- sum(fit$res^2)  # 1
    ss  <- sqrt(mean(fit$res^2))
  }
  
  if (fLoss == 'lmrob') {
    control <- lmrob.control(trace.level = 0,         # 0
                             nResample   =  nresamp,  # 500 default
                             tuning.psi = cterho,     # for 85% eff set to 3.443689
                             subsampling = 'simple',  #
                             rel.tol     = 1e-5,      # 1e-7
                             refine.tol  = 1e-5,      # 1e-7
                             k.max       = 2e3,       # 200
                             maxit.scale = 2e3,       # 200
                             max.it      = 2e3)       # 50
    fit  <- lmrob(yy ~ X, control = control)
    
    if (fit$init.S$converged) {
      cf <- fit$coef
      ss <- fit$scale
      vv <- sum(Mpsi(fit$res / ss,
                     cc  = control$tuning.psi,
                     psi = control$psi, deriv = -1))
      cv <- fit$converged
    } else {
      stop('S-estimator did not converge.')
    }
  }
  
  ## Estimated parameters
  intercept <- cf[1]
  slope_par <- cf[2:(freq + 1)]
  spl_par <- cf[-(1:(freq + 1))]
  
  return(list(spl = spl_par, slope = slope_par, intercept = intercept,
              value = vv, scale = ss, conv = cv, traza = tr, matias = aicM))
}



ajustar_uno <- function (y, x, u, t, interac, freq, spl,
                         norder, fLoss) {

    ## Integration step
    dt <- min(diff(t))  # practicamente uniforme
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
    est <- minimizar(y, xx_coef, u, interac, spl, freq, fLoss,
                     norder, pars)

    est$slope_fun <- cov_dec %*% est$slope

    return (est)
}

medida_bondad <- function(nn, scl, val, spl,freq, criterio){
    switch(criterio, # todos verificados
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

ajustar_todos <- function (y, x, u, t, interac, rango_freq,
                           rango_spl, norder, fLoss, criterio = 'bic') {

    opt <- spl_opt <- freq_opt <- fit_opt <- Inf

    for(spl in rango_spl) {
        for(freq in rango_freq) {
            fit <- ajustar_uno(y, x, u, t, interac, freq, spl, norder, fLoss)
            val <- fit$value
            scl <- fit$scale
            nn  <- length(y)
            crt <- medida_bondad(nn, scl, val, spl, freq, criterio)
            if (crt < opt) {
                opt      <- crt
                spl_opt  <- spl
                freq_opt <- freq
                fit_opt  <- fit
            }
            print(c('spl' = spl, 'freq' = freq, 'crit' = crt))
        }
    }
    return(list(fit = fit_opt, spl = spl_opt, freq = freq_opt,
                interac = interac, u = u, t = t))
}

predecir_una <- function (obs, t, intercept, beta, u, eta_coef) {
    ## no param
    uu      <- u
    norder  <- 4
    spl_kn  <- length(eta_coef)
    kns     <- seq(min(uu), max(uu), length = spl_kn - norder + 2)
    base    <- create.bspline.basis(rangeval = range(uu),
                                    norder = norder, breaks = kns)
    spl_uu  <- getbasismatrix(obs$u, base)
    noparam <- (spl_uu * obs$interac) %*% eta_coef
    ## funcional
    dt <- min(diff(t)) # practicamente uniforme
    funcional <- sum(obs$x * beta * dt)
    return(list(verdad = obs$y, predicho = intercept + funcional +
                                    noparam))
}

submuestra <- function (cuales, covariable, interac) {

    ## Elijo/construyo la submuestra (para train/test)

    y  <- tecator$y$Fat[cuales]
    u  <- tecator$y$Protein[cuales]
    ab <- tecator$absorp.fdata
    t  <- ab$argvals
    nn <- length(y)

    if (interac) {
        interac = tecator$y$Water[cuales]
    } else {
        interac = rep(1, nn)
    }
    if (covariable == 'd1') {
        ab1      <- fdata.deriv(ab, nderiv = 1) # primera derivada
        ab1$data <- ab1$data[cuales, ]
        x        <- ab1$data
    }
    if (covariable == 'd2'){
        ab2      <- fdata.deriv(ab, nderiv = 2) # segunda derivada
        ab2$data <- ab2$data[cuales, ]
        x        <- ab2$data
    }
    return(list(y = y, x = x, u = u, t = t, interac = interac, n = nn))
}

residuos <- function(modelo, datos) {

    ## Residuos
    res <- t(sapply(1:datos$n, function (cual) {
        obs_new <- list(y       = datos$y[cual],
                        x       = datos$x[cual, ],
                        u       = datos$u[cual],
                        interac = datos$interac[cual])
        unlist(predecir_una(obs_new,
                            t         = modelo$t,
                            intercept = modelo$fit$intercept,
                            beta      = modelo$fit$slope_fun,
                            u         = modelo$u,
                            eta_coef  = modelo$fit$spl))
    }))

    return(res)
}

performance <- function (test, res, out_test) {
    mm2 <- mad(test$y)^2
    res2 <- res^2
    if (!missing(out_test)) {
        res2_wo <- res[-out_test]^2
        ret <- list(MSPE = mean(res2) / mm2,
                    MedSPE = median(res2) / mm2,
                    MSPE_clean = mean(res2_wo) / mm2)
    } else {
        ret <- list(MSPE = mean(res2) / mm2,
                    MedSPE = median(res2) / mm2)

    }
    return(ret)
}

performance_trim <- function (test, res, out_test) {
    mm2 <- mad(test$y)^2
    res2 <- res^2
    if (!missing(out_test)) {
        res2_wo <- res[-out_test]^2
        ret <- list(MSPE = mean(res2) / mm2,
                    MSPE_trim = media_trim(res2) / mm2,
                    MSPE_clean = mean(res2_wo) / mm2)
    } else {
        ret <- list(MSPE = mean(res2) / mm2,
                    MSPE_trim = media_trim(res2) / mm2)

    }
    return(ret)
}

media_trim <- function(a, alpha=0.1) {
	n <- length(a)
	n0 <- n - floor( alpha * n )
	return( mean( (sort(a))[1:n0] ) )
}




