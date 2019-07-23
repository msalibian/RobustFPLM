# RobustFPLM
## Robust estimators for Functional Partial Linear Models 

This repository contains code to compute the robust MM-estimators for Functional Partial Linear
models in Boente, Salibian-Barrera, Vena (2019).

Because the regression model of the Tecator data includes an interaction between the non-parametric 
term and a scalar covariate (`y = <X, beta> + w g(z) + e`), setting the appropriate design matrices
(including splitting the data into a training and a testing set, and computing the 2nd derivatives
of the absorbance spectra) requires special care. 
Thus, we prepared a suite of specific 
functions and scripts to reproduce our analysis in the paper. See here.

Code to fit the usual Functional Partial linear Model (`y = <X, beta> + g(z) + e`) is here. A 
simple step-by-step example can be found here. 
