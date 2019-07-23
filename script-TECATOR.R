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

## Plot the second derivatives of the spectra
plot(absorp2, lty=1, col='gray')  # derivada segunda

## or

matplot(absorp2$argvals, t(absorp2$data), type = "l", lty = 1,
        xlab = "Wavelength", ylab = "X''(t)", col='gray50')

## Functional boxplot of the 2nd derivatives
fd2 <- fbplot(t(absorp2$data), x = absorp2$argvals, method = "MBD",
       xlab = "Wavelength", ylab = "d(Absorbance, 2) ",
       xlim = c(850, 1050), ylim = c(-0.004, 0.004))

## Functional boxplot of the original spectra
fd0 <- fbplot(t(absorp$data), x = absorp$argvals, method = "MBD",
       xlab = "Wavelength", ylab = "Absorbance",
       xlim = c(850, 1050), ylim = c(2, 5.5))

## Analysis

## Load the necessary functions
source('tecator-specific-functions.R')

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

## Compute the optimal robust estimator, using 'criterio' to
## select best combination of splines bases

set.seed(124) # for lmrob()
best_rb <- ajustar_todos(y = train$y, x = train$x, u = train$u,
                         t = train$t, interac = train$interac,
                         rango_freq = rango_beta,
                         rango_spl = rango_eta, norder = 4, 'lmrob',
                         criterio = criterio)

## Residuals on the training set from the robust fit

res_best_rb_train <- residuos(best_rb, train)
residuos_train <- res_best_rb_train[, 2] - res_best_rb_train[, 1]
names(residuos_train) <- 1:length(residuos_train)

## identify outliers

( atipicos_train <- as.numeric(names(boxplot(residuos_train)$out)) )
# 7  10  11  28  29  31  34  86  89 122 132 140 141

## Build the robust hat(eta)

norder <- 4
spl_kn <- best_rb$spl
kns    <- seq(min(train$u), max(train$u), length = spl_kn - norder + 2)
base   <- create.bspline.basis(rangeval = range(train$u),
                               norder = norder,
                               breaks = kns)
spl_uu <- getbasismatrix(train$u, base)
eta_est_rob <- spl_uu %*% best_rb$fit$spl

## Plot robust hat(eta)

u_order <- order(train$u)
plot(train$u[u_order], eta_est_rob[u_order], type = "l", xlab = "Protein", pch = 1,
     ylim=c(-1.05,-0.85),
     lwd = 6,
     col = "black", ylab = expression(hat(eta)),lty=1)

## Plot robust hat(beta)

plot(train$t, best_rb$fit$slope_fun, col = "black", lwd = 6,
     type = "l", xlab = "Wavelength", ylab = expression(hat(beta)),
     ylim = c(-6000, 6000),lty=1)


## Residuals & performance of the robust fit on the test set

res_RB <- residuos(best_rb, test)
residuos_test <- res_RB[, 2] - res_RB[, 1]
names(residuos_test) <- 1:length(residuos_test) #156:215
 ( atipicos_test <- as.numeric(names(boxplot(residuos_test)$out)) ) + 155

# 173 174 177 180 181 185

( performance_RB <- performance(test, residuos_test, atipicos_test) )


## Compute the least squares estimator

best_cl <- ajustar_todos(y = train$y, x = train$x,
                         u = train$u, t = train$t,
                         interac = train$interac,
                         rango_freq = rango_beta, rango_spl = rango_eta,
                         norder = 4, 'ls', criterio = criterio)

## Build  the estimate of eta

norder <- 4
spl_kn <- best_cl$spl
kns    <- seq(min(train$u), max(train$u), length = spl_kn - norder + 2)
base   <- create.bspline.basis(rangeval = range(train$u),
                               norder = norder,
                               breaks = kns)
spl_uu <- getbasismatrix(train$u, base)
eta_CL <- spl_uu %*% best_cl$fit$spl

## Plot robust and LS hat(eta) on the same plot

plot(train$u[u_order], eta_est_rob[u_order], type = "l", xlab = "Protein", pch = 1,
     ylim=c(-1.05,-0.85),
     lwd = 6,
     col = "black", ylab = expression(hat(eta)),lty=1)
lines(train$u[u_order], eta_CL[u_order], lwd = 6,
      col = "gray70", ylab = expression(hat(eta)),lty=1)
legend('topright', lwd=6, col=c('black', 'gray70'), legend=c('MM', 'LS'))

## Plot robust and LS hat(beta) on the same plot


plot(train$t, best_rb$fit$slope_fun, col = "black", lwd = 6,
     type = "l", xlab = "Wavelength", ylab = expression(hat(beta)),
     ylim = c(-6000, 6000),lty=1)
lines(train$t, best_cl$fit$slope_fun,col = "gray70",
      lwd = 6, lty = 1)
legend('top', lwd=6, col=c('black', 'gray70'), legend=c('MM', 'LS'))

## Performance of the LS estimator on the test set

res_CL <- residuos(best_cl, test)
residuos_CL <- res_CL[, 2] - res_CL[, 1]
( performance_CL <- performance(test, residuos_CL, atipicos_test) )

## Performance table

rbind('CL' = performance_CL, 'RB' = performance_RB)

## Remove outliers from the training set

train_sin_out <- submuestra((1:155)[-atipicos_train], covariable='d2', interac=TRUE)

## Best LS fit on "clean" training set

best_cl_sin_out <- ajustar_todos(y = train_sin_out$y, x = train_sin_out$x,
                                 u = train_sin_out$u, t = train_sin_out$t,
                                 interac = train_sin_out$interac,
                                 rango_freq = rango_beta, rango_spl = rango_eta,
                                 norder = 4, 'ls', criterio = criterio)

## Build  hat(eta) based on  "clean" training set

norder <- 4
spl_kn <- best_cl_sin_out$spl
kns    <- seq(min(train_sin_out$u), max(train_sin_out$u), length = spl_kn - norder + 2)
base   <- create.bspline.basis(rangeval = range(train_sin_out$u),
                               norder = norder,
                               breaks = kns)
spl_uu <- getbasismatrix(train_sin_out$u, base)
eta_CL_sin_out <- spl_uu %*% best_cl_sin_out$fit$spl

## Plot robust and LS hat(eta) on clean data

uu_order <- order(train_sin_out$u)
plot(train$u[u_order], eta_est_rob[u_order], type = "l", xlab = "Protein", pch = 1,
     ylim=c(-1.05,-0.85),
     lwd = 6,
     col = "black", ylab = expression(hat(eta)),lty=1)
lines(train_sin_out$u[uu_order], eta_CL_sin_out[uu_order], lwd = 6,
      col = "gray70", ylab = expression(hat(eta)),lty=1)
legend('topright', lwd=6, col=c('black', 'gray70'), legend=c('MM', 'LS w/o outliers'))


## Plot robust and LS hat(beta) on clean data

plot(train$t, best_rb$fit$slope_fun, col = "black", lwd = 6,
     type = "l", xlab = "Wavelength", ylab = expression(hat(beta)),
     ylim = c(-6000, 6000),lty=1)
lines(train_sin_out$t, best_cl_sin_out$fit$slope_fun,col = "gray70",
      lwd = 6, lty = 1)
legend('topleft', lwd=6, col=c('black', 'gray70'), legend=c('MM', 'LS w/o outliers'))

