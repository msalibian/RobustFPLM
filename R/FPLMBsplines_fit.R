#' FPLM fit with splines
#'
#' @description
#' FPLM fit with splines
#'
#' 
#' @param y vector of n responses
#' @param x matrix of functional covariates (n x p_1)
#' @param u vector of n observations for g(u_i)
#' @param t vector of "times" where X was observed "x[i, j] = X_i( t_j )"
#' @param freq description
#' @param spl description
#' @param norder description
#' @param fLoss loss function to be minimize (see Details)
#' @return A list including several.
#'
#' @details
#' Valid options for the loss function are 'ls', 'lad', 'lmrob', 'huang'.
#' 
#' @examples
#' Example.
#' 
#' @references Paper
#' @export
FPLMBsplines_fit <- function(y, x, u, t, freq, spl, norder, fLoss) {

    ## Integration step
    dt <- min(diff(t)) # width of grid
    xcenter <- x

    ## Covariance decomposition
    grilla_spl <- t # seq(0, 1, length = length(t))
    nodos_spl <- seq(min(grilla_spl), max(grilla_spl),
                     length = freq - norder + 2)
    base_spl <- create.bspline.basis(
        rangeval = range(t),
        norder = norder,
        breaks = nodos_spl
    )
    beta_spl <- getbasismatrix(grilla_spl, base_spl)
    cov_dec <- beta_spl[, 1:freq]

    ## Estimated Fourier coefficients (by row)
    xx_coef <- xcenter %*% cov_dec * dt

    ## Parameter estimation
    est <- minimize(
        y, xx_coef, u, spl, freq, fLoss,
        norder, pars
    )

    est$slope_fun <- cov_dec %*% est$slope

    return(est)
}
