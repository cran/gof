/*!
  @file utils.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief Various utility functions and constants

*/
#include "utils.hpp"

namespace cumres {
  

  arma::mat softmax(arma::mat &lp, bool ref=true, bool log=false) {
    if (ref) lp.insert_cols(0, arma::zeros(lp.n_rows));
    arma::colvec lpmax = arma::max(lp, 1);
    lp.each_col() -= lpmax;
    arma::colvec denom = sum(exp(lp), 1);
    lp.each_col() -= arma::log(denom);
    if (log) return(lp);
    return(exp(lp));
  }


  arma::mat deriv(cx_func f, arma::vec theta) {
    arma::cx_vec thetac = arma::conv_to<arma::cx_vec>::from(theta);
    arma::cx_mat val0 = f(thetac);
    unsigned n = val0.n_elem;
    unsigned p = theta.n_elem;
    arma::mat res(n,p);
    double h = DBL_MIN;
    cx_dbl h0 = cx_dbl(0, h);
    for (unsigned i=0; i<p; i++) {
      arma::cx_vec theta0 = thetac;
      theta0[i] += h0;
      arma::mat val = imag(f(theta0))/h;
      for (unsigned j=0; j<n; j++)
	res(j,i) = val[j];      
    }
    return(res);
  }  

  // template arma::mat expit<double>(const arma::mat&);
  // template arma::cx_mat expit<Complex>(const arma::cx_mat&);

  arma::mat expit(arma::mat x) {
    for (unsigned i=0; i<x.n_elem; i++) {
      double z = x(i);
      if (z>=0) {
	x(i) = 1/(1+exp(-z));
      } else {
	z = exp(z);
	x(i) = z/(1+z);
      }
    }
    return(x);
  }

  arma::cx_mat expit(arma::cx_mat x) {
    return 1.0/(1+exp(-x));
  }

  arma::vec softmax(arma::vec u) {
    double umax = u.max();
    u -= umax;
    double denom = sum(exp(u));  
    return u - log(denom);
  }

  double KolmogorovSmirnov(const arma::vec &x) {
    return arma::max(arma::abs(x));
  };

  double CramerVonMises(const arma::vec &x, const arma::vec &t) {
    arma::vec delta(t.n_elem);
    for (unsigned i=0; i<t.n_elem-1; i++) delta(i) = t[i+1]-t[i];
    delta(delta.n_elem-1) = 0;
    return std::sqrt(sum(delta % x % x));
  }

  arma::mat const EmptyMat = arma::mat();
  arma::vec const EmptyVec = arma::vec();  

  
  // Foreground colors are in form of 3x, bacground are 4x
  const char* COL_RESET = "\x1b[0m";
  const char* COL_DEF   = "\x1b[39m";
  const char* BLACK     = "\x1b[30m";
  const char* RED       = "\x1b[31m";
  const char* MAGENTA   = "\x1b[35m";
  const char* YELLOW    = "\x1b[33m";
  const char* GREEN     = "\x1b[32m";
  const char* BLUE      = "\x1b[34m";
  const char* CYAN      = "\x1b[36m";
  const char* WHITE     = "\x1b[37m";
  const char* GRAY      = "\x1b[90m";
  const char* LRED      = "\x1b[91m";
  const char* LGREEN    = "\x1b[92m";
  const char* LYELLOW   = "\x1b[93m";
  const char* LBLUE     = "\x1b[94m";
  const char* LMAGENTA  = "\x1b[95m";
  const char* LCYAN     = "\x1b[96m";
  const char* LWHITE    = "\x1b[97m";

} // namespace cumres
