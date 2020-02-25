/*!
  @file rinterface.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief R interface for the cumres class. 

  The relevant bindings are created in \c RcppExports.cpp, \c RcppExports.h
  after running \c Rcpp::compileAttributes()
*/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

#include "cumres.hpp"
#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::export]]
double KolmogorovSmirnov(const arma::vec &x) {
  return cumres::KolmogorovSmirnov(x);
}

// [[Rcpp::export]]
double CramerVonMises(const arma::vec &x, const arma::vec &t) {
  return cumres::CramerVonMises(x, t);
}

RCPP_MODULE(gofmod) {
    using namespace Rcpp;
    class_<cumres::cumres>("CumRes")
    // expose the constructor
      .constructor<arma::vec, 
		   arma::mat,
		   arma::mat>("Constructor")
      .field( "t", &cumres::cumres::t )
      .field( "ord", &cumres::cumres::ord )
      .field( "r", &cumres::cumres::r )
      .field( "qt", &cumres::cumres::qt )

      .method("samplestat", (arma::mat (cumres::cumres::*)(unsigned, arma::uvec, bool ) )( &cumres::cumres::sample),
      	       "sample process and return KS and CvM statistic")
      .method("sample1", (arma::vec (cumres::cumres::*)(arma::uvec) )( &cumres::cumres::sample),
	       "sample process")      
      .method("obs",    &cumres::cumres::obs,   "Return observed process")      
      .method("rnorm",  &cumres::cumres::rnorm,   "Sample from Gaussian")
      .method("reorder",  &cumres::cumres::reorder,  "Order observations after input variable")
      ;
}

