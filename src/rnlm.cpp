
#include <RcppArmadillo.h>
#include "minimize.h"
#include "pl2_item.h"
#include "nrm_item.h"

using namespace arma;

// [[Rcpp::export]]
void test_nlm(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par_in)
{
	ll_pl2 f(a.memptr(), theta.memptr(), r); 
	const double gtol=1e-10;
	int iter=0,err=0;
	double fret;
	vec par = par_in;

	nlm(par, gtol, iter, fret, f, err);
	
	printf("par: %f, %f\nvalue at min: %f\niter: %i\nerror code: %i\n",par[0],par[1],fret, iter,err);
	fflush(stdout);
}

// [[Rcpp::export]]
void test_D1(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par_in)
{
	const int max_a = max(a), nt = theta.n_elem;
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	ll_nrm f(a.memptr(), exp_at, r);
	
	const double tol=1e-10;
	int iter=0,err=0;
	double fret;
	vec par = par_in;
	D1min(par, tol, iter, fret, f, err); 
	
	printf("par: %f\nvalue at min: %f\niter: %i\nerror code: %i\n",par[0],fret, iter,err);
	fflush(stdout);
}
