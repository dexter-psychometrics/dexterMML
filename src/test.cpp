
#include <RcppArmadillo.h>
#include "item_ll.h"
#include "minimize.h"

using namespace arma;
using Rcpp::Named;

// [[Rcpp::export]]
double test_ll(arma::ivec& a, arma::vec& b, arma::vec& theta, arma::mat& r)
{
	ll_nrm f(a, theta, r);
	return f(b);
}

// [[Rcpp::export]]
arma::vec test_df(arma::ivec& a, arma::vec& b, arma::vec& theta, arma::mat& r)
{
	ll_nrm f(a, theta, r);
	vec g(b.n_elem);
	f.df(b,g);
	return g;
}

// [[Rcpp::export]]
Rcpp::List test_minimize(arma::ivec& a, arma::vec& b, arma::vec& theta, arma::mat& r)
{
	ll_nrm f(a, theta, r);
	vec pars = b;
	
	const double tol = 1e-8;
	int itr=0;
	double ll_itm=0;

	dfpmin(pars, tol, itr, ll_itm, f);

	return Rcpp::List::create(Named("ll")=ll_itm,Named("iter")=itr,Named("b")=pars);
}