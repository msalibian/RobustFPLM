#' Best FPLMBsplines_fit given by a model selection criterion.
#'
#' @description Fit a FPLM model for different spline basis sizes and picks the
#'     best one according to a specified model selection criterion.
#' @param y the vector of scalar responses.
#' @param x a matrix of the functional covariates, where each row contains the
#'     functions evaluated on a (common) grid.
#' @param u the values of the explanatory variable that enters the model
#'     non-parametrically.
#' @param t the grid over which the functional covariates were evaluated.
#' @param range_freq a vector of B-spline basis sizes to try for the functional
#'     regression coefficient.
#' @param range_spl a vector of B-spline basis sizes to try for the
#'     non-parametric component.
#' @param norder the order of the B-Splines.
#' @param fLoss loss function to be minimized.
#' @param criterion criterion for model selection.
#' @param trace a logical argument indicating whether partial results are
#'     printed.
#'
#' @details Add some details regarding the criteria, valid ranges (asympotics)
#'     and floss options
#'
#' @return A list including the following components:
#' \itemize{
#' \item{fit}{fitted parameters}
#' \item{spl}{chosen number of splines for the non-parametric component}
#' \item{freq}{chosen number of splines for the funcitonal regression coefficient}
#' }
#' @examples
#'
#' n <- 100
#' m <- 50
#' u <- runif(n)
#' t <- runif(m)
#' b <- function(x) x^3
#' g <- function(x) sin(x)
#' x <- matrix(rnorm(n * m), nrow = n)
#' y <- x %*% b(t) * min(diff(t)) + g(u) + rnorm(n, sd = 0.1)
#'
#' FPLM_fit <- FPLMBsplines(y, x, u, t,
#'   range_freq = 4:13, range_spl = 4:13,
#'   norder = 4, fLoss = "ls", criterion = "bic1", trace = FALSE
#' )
#'
#' par(mfrow = c(2, 1))
#' plot(t, FPLM_fit$fit$slope_fun, pch = 16)
#' plot(u, FPLM_fit$fit$eta_est, pch = 16)
#' @export
FPLMBsplines <- function(y, x, u, t, range_freq = range_default,
                         range_spl = range_default, norder = 4,
                         fLoss = "lmrob", criterion = "bic1",
                         trace = FALSE) {

    ## Some Setup
    opt <- spl_opt <- freq_opt <- fit_opt <- Inf
    n <- length(y)
    range_freq <- range_spl <- floor(max(n^(1 / 5), norder)):
        floor(2 * (norder + n^(1 / 5)))

    ## Double loop
    for (spl in range_spl) {
        for (freq in range_freq) {
            fit <- FPLMBsplines_fit(y, x, u, t, freq, spl, norder, fLoss)
            val <- fit$value
            scl <- fit$scale
            crt <- goodness(n, scl, val, spl, freq, criterion)
            if (crt < opt) {
                opt <- crt
                spl_opt <- spl
                freq_opt <- freq
                fit_opt <- fit
            }
            if (trace) print(c("spl" = spl, "freq" = freq, "crit" = crt))
        }
    }

    ## Best fit
    kns <- seq(min(u), max(u), length = spl_opt - norder + 2)
    base <- create.bspline.basis(
        rangeval = range(u),
        norder = norder,
        breaks = kns
    )
    spl_uu <- getbasismatrix(u, base)
    fit_opt$eta_est <- spl_uu %*% fit_opt$spl
    dt <- min(diff(t))
    fit_opt$fitted <- as.vector(x %*% fit_opt$slope_fun * dt + fit_opt$eta_est)

    return(list(fit = fit_opt, spl = spl_opt, freq = freq_opt))
}
