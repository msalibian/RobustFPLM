#' Minimization rutine
#'
#' @description
#' This function is used internally by FPLMBsplines_fit(), not meant
#' for 'public' use.
#'
#' @param y the vector of scalar responses.
#' @param x_coef covariate coefficients.
#' @param u the values of the explanatory variable that enters the model
#'     non-parametrically.
#' @param spl_kn number of splines for the non-parametric term.
#' @param freq number of splines for the functional term.
#' @param fLoss string specifying the loss function. 'ls' for least squares,
#'     'huang' for Huber, 'lmbrob' for MM-estimator.
#' @param norder spline order (n=4 for cubic splines)
#' @param rob.control optional parameters if \code{fLoss == "lmbrob"}.
#'
#' @return A vector with components:
#' \itemize{
#'   \item{spl}{estimate for the non parametric term},
#'   \item{slope}{estimate for the functional slope},
#'   \item{value}{the minimum value of the loss},
#'   \item{scale}{ the scale estimate.} 
#' }
#' 
#' @import fda robustbase
#' @importFrom stats lm
minimize <- function(y, x_coef, u, spl_kn, freq, fLoss, norder,
                     rob.control = lmrob.control(
                         trace.level = 0,
                         nResample = 5000, 
                         tuning.psi = 3.443689, # 85% eff
                         subsampling = "simple",
                         rel.tol = 1e-5, 
                         refine.tol = 1e-5,
                         k.max = 2e3,
                         maxit.scale = 2e3,
                         max.it = 2e3
                     )) {

    ## B-spline basis
    kns <- seq(min(u), max(u), length = spl_kn - norder + 2)
    base <- create.bspline.basis(
        rangeval = range(u),
        norder = norder,
        breaks = kns
    )
    spl_u <- getbasismatrix(u, base)

    ## Design matrix
    X <- cbind(x_coef, spl_u)

    if (!(fLoss %in% c("ls", "huang", "lmrob"))) {
        stop("Invalid fLoss. Should be one of \'ls\' or \'huang\' or \'lmrob\' ")
    }

    ## Minimization menu
    switch(fLoss,
           ls = {
               fit <- lm(y ~ X - 1)
               cf <- fit$coef
               vv <- sum(fit$res^2) # 1
               ss <- sqrt(mean(fit$res^2))
           },
           huang = {
               init <- lm(y ~ X - 1)$coef
               fit <- huber_estimates(X, y, init, 1.345, 1e-8)
               cf <- as.vector(fit$param)
               ss <- 1
               vv <- fit$value
           },
           lmrob = {
               fit <- lmrob(y ~ X - 1, control = rob.control)
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
