#' FPLM fit with splines
#'
#' @description
#' FPLM fit with splines
#'
#' @param y the vector of scalar responses.
#' @param x a matrix of the functional covariates, where each row contains the
#'     functions evaluated on a (common) grid.
#' @param u the values of the explanatory variable that enters the model
#'     non-parametrically.
#' @param t the grid over which the functional covariates were evaluated.
#' @param freq basis size for the functional regression coefficient.
#' @param spl basis size for the non-parametric component.
#' @param norder the order of the B-Splines.
#' @param fLoss string specifying the loss function. 'ls' for least squares,
#'     'huang' for Huber, 'lmbrob' for MM-estimator.
#'
#' @return A list with components:
#' \itemize{
#'   \item{spl}{estimated coefficients for the non parametric term},
#'   \item{slope}{estimated coefficients for the functional slope},
#'   \item{value}{the minimum value of the loss},
#'   \item{scale}{the scale estimate},
#'   \item{slope_fun}{functional slope estimate}.
#' }
#'
#' @import fda robustbase
#' 
#' @references Pending.
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
        norder
    )

    est$slope_fun <- cov_dec %*% est$slope

    return(est)
}
