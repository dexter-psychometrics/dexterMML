#include <RcppArmadillo.h>
#include "nrm_item.h"
using namespace arma;

mat nrm_trace(const vec& theta, const ivec& a, const vec& b, const int ncat, const mat& exp_at)
{
	const int nt = theta.n_elem;
	mat out(nt,ncat);		
	
	vec exp_b = exp(b);
	
	for(int j=0; j< nt; j++)
	{
		double sm=0;
		for(int k=0;k<ncat;k++)
			sm += exp_b[k]*exp_at.at(a[k],j);
		for(int k=0;k<ncat;k++)
			out.at(j,k) = exp_b[k]*exp_at.at(a[k],j)/sm;
	}
		
	return out;
}


// [[Rcpp::export]]
double test_ll_nrm(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par)
{
	const int max_a = max(a), nt = theta.n_elem;
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	ll_nrm f(a.memptr(), exp_at, r);
	
	return f(par);
} 

// [[Rcpp::export]]
arma::vec test_gradient_nrm(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par)
{
	const int max_a = max(a), nt = theta.n_elem;
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	ll_nrm f(a.memptr(), exp_at, r);
	vec g(par.n_elem);
	f.df(par,g);
	return g;
}

// [[Rcpp::export]]
arma::mat test_hess_nrm(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par)
{
	const int max_a = max(a), nt = theta.n_elem;
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	ll_nrm f(a.memptr(), exp_at, r);
	mat h(par.n_elem,par.n_elem);
	f.hess(par,h);
	return h;
}