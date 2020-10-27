#include <RcppArmadillo.h>
#include "poly2pl_item.h"
using namespace arma;

mat poly2_trace(const vec& theta, const ivec& a, const double A, const vec& b, const int ncat)
{
	const int nt = theta.n_elem;
	mat out(nt,ncat);	
	vec p(ncat);
	p[0] = 1;
	
	for(int t=0; t<nt; t++)
	{
		double s=1;		
		for(int k=1; k<ncat; k++)
		{
			p[k] = std::exp(A*a[k]*(theta[t]-b[k])); 
			s += p[k];
		}

		for(int k=0; k<ncat; k++)
			out.at(t,k) = p[k]/s;
	}

	return out;
}


// [[Rcpp::export]]
double test_ll_p2(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par)
{
	ll_poly2 f(a.memptr(), theta.memptr(), r); 
	
	return f(par);
} 

// [[Rcpp::export]]
arma::vec test_gradient_p2(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par)
{
	ll_poly2 f(a.memptr(), theta.memptr(), r); 
	vec g(par.n_elem);
	f.df(par,g);
	return g;
}

// [[Rcpp::export]]
arma::mat test_hess_p2(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par)
{
	ll_poly2 f(a.memptr(), theta.memptr(), r); 
	mat h(par.n_elem,par.n_elem);
	f.hess(par,h);
	return h;
}