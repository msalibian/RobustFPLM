# RobustFPLM
## Robust estimators for Functional Partial Linear Models 

This repository contains code to compute the robust MM-estimators for Functional Partial Linear
models in Boente, Salibian-Barrera, Vena (2019).

Because the regression model of the Tecator data includes an interaction between the non-parametric 
term and a scalar covariate (`y = <X, beta> + w g(z) + e`), we prepared a suite of specific functions and scripts to reproduce 
our analysis in the paper. See here.

Code to fit the usual Functional Partial linear Model (`y = <X, beta> + g(z) + e`), is here. A 
simple step-by-step example can be found here. 
