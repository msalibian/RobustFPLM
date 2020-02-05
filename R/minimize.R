#' Minimization rutine
#'
#' @description
#' This function is used internally by FPLMBsplines(), not meant
#' for 'public' use.
#'
#' @param yy description
#' @param xx_coef description
#' @param uu description
#' @param spl_kn description
#' @param freq description
#' @param fLoss description
#' @param norder description
#' @param pars description
#' @param rob.control description
#'
#' @return Minimization stuff
#' @export
minimize <- function(yy, xx_coef, uu, spl_kn, freq, fLoss, norder, pars,
                     rob.control = lmrob.control(
                         trace.level = 0, # 0
                         nResample = 5000, # 500 default
                         tuning.psi = 3.443689, # 85%
                         subsampling = "simple", #
                         rel.tol = 1e-5, # 1e-7
                         refine.tol = 1e-5, # 1e-7
                         k.max = 2e3, # 200
                         maxit.scale = 2e3, # 200
                         max.it = 2e3
                     )) {

    ## B-spline basis
    kns <- seq(min(uu), max(uu), length = spl_kn - norder + 2)
    base <- create.bspline.basis(
        rangeval = range(uu),
        norder = norder,
        breaks = kns
    )
    spl_uu <- getbasismatrix(uu, base)

    ## Design matrix

    X <- cbind(xx_coef, spl_uu)

    if (!(fLoss %in% c("ls", "huang", "lmrob"))) {
        stop("Invalid fLoss. Should be one of \'ls\' or \'huang\' or \'lmrob\' ")
    }

    ## Minimization menu
    switch(fLoss,
           ls = {
               fit <- lm(yy ~ X - 1)
               cf <- fit$coef
               vv <- sum(fit$res^2) # 1
               ss <- sqrt(mean(fit$res^2))
           },
           huang = {
               init <- lm(yy ~ X - 1)$coef
               fit <- huber_estimates(X, yy, init, 1.345, 1e-8)
               cf <- as.vector(fit$param)
               ss <- 1
               vv <- fit$value
           },
           lmrob = {
               fit <- lmrob(yy ~ X - 1, control = rob.control)
               if (fit$init.S$converged) {
                   cf <- fit$coef
                   ss <- fit$scale
                   vv <- sum(Mpsi(fit$res / ss,
                                  cc = rob.control$tuning.psi,
                                  psi = rob.control$psi, deriv = -1
                                  ))
               } else {
                   stop("S-estimator did not converge.")
               }
           },
           stop("Invalid fLoss.")
           )

    spl_par <- cf[-(1:freq)]
    slope_par <- cf[  1:freq ]

    return(list(
        spl = spl_par,
        slope = slope_par,
        value = vv,
        scale = ss
    ))
}
