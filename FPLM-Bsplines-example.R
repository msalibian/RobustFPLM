## ---------------------------------------------------------------------
## TECATOR - http://lib.stat.cmu.edu/datasets/tecator
##

library('robustbase') # lmrob

## Load data
library(fda.usc)
datos <- data(tecator)
absorp <- tecator$absorp.fdata
## splines approximation to the second derivative
absorp2 <- fdata.deriv(absorp,nderiv = 2)

## Load the necessary functions
source('tecator-specific-functions.R')

source('FPLM-Bsplines-functions.R')

## Model: fat = < beta, absorp2> + water * eta(protein)
## splines dimension selection criterion
criterio <- 'bic1'
## range for possible dimension for splines for the regression component
rango_beta <- 4:13

## range for possible dimension for splines for the non-parametric component

rango_eta <- 4:13

## Create the training and test data sets
## "interact" indicates whether the 2nd term in the model is
## "water * g(protein)" (interac == TRUE) or
## "g(protein)" (interac == FALSE)
## covariable = 'd1' or 'd2' indicates which derivativate will
## be returned
##
## This function is designed to work with the TECATOR data
## as available in the package fda.usc

train <- submuestra(1:155, covariable = 'd2', interac = TRUE)
test  <- submuestra(156:215, covariable = 'd2', interac = TRUE)

set.seed(124) # for lmrob()
a <- FPLMBsplines(y = train$y, x = train$x, u = train$u,
                  t = train$t, range_freq = 4:13,
                  range_spl = 4:13, norder = 4, fLoss='lmrob', trace=TRUE)


d <- FPLMBsplines(y = train$y, x = train$x, u = train$u,
                  t = train$t, range_freq = 4:13,
                  range_spl = 4:13, norder = 4, fLoss='huang', trace=TRUE)



norder <- 4
spl_kn <- a$spl
kns    <- seq(min(train$u), max(train$u), length = spl_kn - norder + 2)
base   <- create.bspline.basis(rangeval = range(train$u),
                               norder = norder,
                               breaks = kns)
spl_uu <- getbasismatrix(train$u, base)
eta_est_rob <- spl_uu %*% a$fit$spl

## Plot robust hat(eta)

u_order <- order(train$u)
plot(train$u[u_order], eta_est_rob[u_order], type = "l", xlab = "Protein", pch = 1,
     lwd = 6,
     col = "black", ylab = expression(hat(eta)),lty=1)

norder <- 4
spl_kn <- d$spl
kns    <- seq(min(train$u), max(train$u), length = spl_kn - norder + 2)
base   <- create.bspline.basis(rangeval = range(train$u),
                               norder = norder,
                               breaks = kns)
spl_uu <- getbasismatrix(train$u, base)
eta_est_huber <- spl_uu %*% d$fit$spl
lines(train$u[u_order], eta_est_huber[u_order], lwd = 6, col = "tomato3")




## Plot robust hat(beta)

plot(train$t, a$fit$slope_fun, col = "black", lwd = 6,
     type = "l", xlab = "Wavelength", ylab = expression(hat(beta)),
     ylim = c(-6000, 6000),lty=1)


b <- FPLMBsplines(y = train$y, x = train$x, u = train$u,
                  t = train$t, range_freq = 4:13,
                  range_spl = 4:13, norder = 4, fLoss='ls')


norder <- 4
spl_kn <- b$spl
kns    <- seq(min(train$u), max(train$u), length = spl_kn - norder + 2)
base   <- create.bspline.basis(rangeval = range(train$u),
                               norder = norder,
                               breaks = kns)
spl_uu <- getbasismatrix(train$u, base)
eta_est_ls <- spl_uu %*% b$fit$spl



u_order <- order(train$u)
plot(train$u[u_order], eta_est_rob[u_order], type = "l", xlab = "Protein", pch = 1,
     lwd = 6,
     col = "black", ylab = expression(hat(eta)),lty=1)
lines(train$u[u_order], eta_est_ls[u_order], type='l', lwd=6, col='tomato3')
lines(train$u[u_order], eta_est_huber[u_order], lwd = 6, col = "skyblue3")

plot(train$t, b$fit$slope_fun, col = "black", lwd = 6,
     type = "l", xlab = "Wavelength", ylab = expression(hat(beta)),
     ylim = c(-6000, 6000),lty=1)

