
// (c) Pablo Vena (https://github.com/anevolbap)

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec huber_weight(const arma::vec &x, const double &cw) {
  int n = x.size();
  arma::vec out = zeros(n);
  for (int iter = 0; iter < n; iter++) {
    if (fabs(x[iter]) <= cw)
      out[iter] = 1;
    else 
      out[iter] = cw / fabs(x[iter]);
  }
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec huber_rho(const arma::vec &x, const double &cw) {
  int n = x.size();
  arma::vec out = zeros(n);
  for (int iter = 0; iter < n; iter++) {
    if (fabs(x[iter]) <= cw)
      out[iter] = pow(x[iter], 2)/2;
    else 
      out[iter] = cw * fabs(x[iter]) - pow(cw,2)/2;
  }
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List huber_estimates(const arma::mat &X, const arma::vec &y, arma::vec &beta,
		     const double &cw, const double &tol) {
  int n = y.size();
  arma::vec beta_old = beta;
  arma::vec residuals = zeros(n);
  arma::mat weights = zeros(n, n);
  double diff = tol;
  while (diff >= tol) {		
    residuals = y - X * beta;
    weights = diagmat(huber_weight(residuals, cw));
    arma::mat tX = trans(X);
    diff = norm(beta - beta_old, 2)/norm(beta_old, 2);
    beta = solve(tX * weights * X, tX * weights * y);
  }
  double value = accu(huber_rho(residuals, cw));
  List out(2);
  out("value") = value;
  out("param") = beta;
  return out;
}
