#' Best FPLMBsplines fit by criterion.
#'
#' @param y response
#' @param x description
#' @param u description
#' @param t description 
#' @param range_freq description
#' @param range_spl description
#' @param norder Order of the splines.
#' @param fLoss Loss function
#' @param criterion Criterion for model selection
#' @param trace Logical.
#'
#' @return list(fit = fit_opt, spl = spl_opt, freq = freq_opt, u, t)
#'
#' @examples List
#'
#' @export
FPLMBsplines <- function (y, x, u, t, range_freq = range_default,
                          range_spl = range_default, norder = 4,
                          fLoss = 'lmrob', criterion = 'bic1',
                          trace = FALSE) {

    ## Some Setup
    opt <- spl_opt <- freq_opt <- fit_opt <- Inf
    n <- length(y)
    range_freq <- range_spl <- floor(max(n^(1/5), norder)):
        floor(2 * (norder + n^(1/5)))

    ## Double loop
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

    ## Best fit
    kns    <- seq(min(u), max(u), length = spl_opt - norder + 2)
    base   <- create.bspline.basis(rangeval = range(u),
                                   norder = norder,
                                   breaks = kns)
    spl_uu <- getbasismatrix(u, base)
    fit_opt$eta_est <- spl_uu %*% fit_opt$spl
    dt <- min(diff(t))
    fit_opt$fitted <- as.vector(x %*% fit_opt$slope_fun * dt + fit_opt$eta_est)

    return(list(fit = fit_opt, spl = spl_opt, freq = freq_opt,
                u = u, t = t))
}
