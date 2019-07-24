## ---------------------------------------------------------------------
## TECATOR - http://lib.stat.cmu.edu/datasets/tecator
##

library('robustbase') # lmrob
source('FPLM-Bsplines-functions.R')

## Load data
library(fda.usc)
data(tecator)
## Prepare training set
## splines approximation to the second derivative
absorp2 <- fdata.deriv(tecator$absorp.fdata, nderiv = 2)
y  <- tecator$y$Fat[1:155]
u  <- tecator$y$Protein[1:155]
t  <- absorp2$argvals
x  <- absorp2$data[1:155,]

## Model: fat = < beta, absorp2> + eta(u)
criterion <- 'bic1'
## range for possible dimension for splines for the regression component
range_beta <- 4:13
## range for possible dimension for splines for the non-parametric component
range_eta <- 4:13
norder <- 4
# test  <- submuestra(156:215, covariable = 'd2', interac = TRUE)

# MM estimator
set.seed(124) # for lmrob()
a <- FPLMBsplines(y = y, x = x, u = u, t = t, range_freq = range_beta,
                  range_spl = range_eta, norder = norder, fLoss='lmrob', trace=TRUE)
# M- (Huang)
d <- FPLMBsplines(y = y, x = x, u = u, t = t, range_freq = range_beta,
                  range_spl = range_eta, norder = norder, fLoss='huang', trace=TRUE)
# LS
b <- FPLMBsplines(y = y, x = x, u = u, t = t, range_freq = range_beta,
                  range_spl = range_eta, norder = norder, fLoss='ls', trace=TRUE)

## Plot hat(eta)
u_order <- order(u)
ma <- max(rra <- c(a$fit$eta_est, b$fit$eta_est, d$fit$eta_est) )
mi <- min(rra)
plot(u[u_order], a$fit$eta_est[u_order], type = "l", xlab = "Protein", pch = 1,
     lwd = 6, col = "black", ylab = expression(hat(eta)), lty=1, ylim=c(mi, ma))
lines(u[u_order], d$fit$eta_est[u_order], lwd = 6, col = "tomato3")
lines(u[u_order], b$fit$eta_est[u_order], lwd = 6, col = "skyblue3")
legend('bottomright', lty=1, lwd=3, col=c('black', 'tomato3', 'skyblue3'),
       legend=c('MM', 'M', 'LS'))
## Plot robust hat(beta)
ma <- max(rrb <- c(a$fit$slope_fun, b$fit$slope_fun, d$fit$slope_fun) )
mi <- min(rrb)
plot(t, a$fit$slope_fun, col = "black", lwd = 6,
     type = "l", xlab = "Wavelength", ylab = expression(hat(beta)), lty=1, ylim=c(mi, ma))
lines(t, d$fit$slope_fun, col = "tomato3", lwd = 6, lty=1)
lines(t, b$fit$slope_fun, col = "skyblue3", lwd = 6, lty=1)
legend('bottom', lty=1, lwd=3, col=c('black', 'tomato3', 'skyblue3'),
       legend=c('MM', 'M', 'LS'))

