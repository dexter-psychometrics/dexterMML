#include <RcppArmadillo.h>
#include "pl2_item.h"
using namespace arma;

//in place
void pl2_trace(const vec& theta, const ivec& a, const double A, const vec& b, const int ncat, mat& out)
{
	const int nt = theta.n_elem;
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
}

mat pl2_trace(const vec& theta, const ivec& a, const double A, const vec& b, const int ncat)
{
	mat out(theta.n_elem,ncat);	
	pl2_trace(theta, a, A, b, ncat, out);

	return out;
}

cube pl2_trace_GH(const mat& theta, const ivec& a, const double A, const vec& b, const int ncat)
{
	const int ng = theta.n_cols;
	cube out(theta.n_rows, ncat, ng);
	for(int g=0;g<ng;g++)
		pl2_trace(theta.col(g), a, A, b, ncat, out.slice(g));
		
	return out;
}


// all constants related to icc
void pl2_icc(const vec& theta, const ivec& a, const double A, const vec& b, const int ncat, 
				mat& itrace, double* nc_ptr, double* nca_ptr, double* ncab_ptr)
{
	const int nt = theta.n_elem;
	vec p(ncat);
	p[0] = 1;
	
	vec norm_const(nc_ptr,ncat,false,true), norm_const_a(nca_ptr,ncat,false,true), norm_const_ab(ncab_ptr,ncat,false,true);
	
	for(int t=0; t<nt; t++)
	{
		double s=1,sa=0,sab=0,sabt=0;		
		for(int k=1; k<ncat; k++)
		{
			p[k] = std::exp(A*a[k]*(theta[t]-b[k])); 
			s += p[k];
			sa += p[k]*a[k];
			sab += p[k]*a[k]*b[k];
			sabt += p[k]*a[k]*(b[k]-theta[t]);
		}

		for(int k=0; k<ncat; k++)
			itrace.at(t,k) = p[k]/s;			
		norm_const[t] = s;
		norm_const_a[t] = sa;
		norm_const_ab[t] = sab;
	}
}


// [[Rcpp::export]]
double test_ll_p2(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par, const int prior=0)
{
	ll_pl2 f(a.memptr(), theta.memptr(), r,prior); 
	
	return f(par);
} 

// [[Rcpp::export]]
arma::vec test_gradient_p2(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par, const int prior=0)
{
	ll_pl2 f(a.memptr(), theta.memptr(), r,prior); 
	vec g(par.n_elem);
	f.df(par,g);
	return g;
}

// [[Rcpp::export]]
arma::mat test_hess_p2(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par, const int prior=0)
{
	ll_pl2 f(a.memptr(), theta.memptr(), r,prior); 
	mat h(par.n_elem,par.n_elem);
	f.hess(par,h);
	return h;
}


// [[Rcpp::export]]
double test_ll_v2(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta, 
			const arma::ivec& ip, const arma::ivec& pi, const arma::ivec& pcni, const arma::ivec& px, 
			const arma::ivec& pgroup, const arma::ivec& inp, const arma::ivec& icnp,
			const arma::vec& mu, const arma::vec& sigma, const int item, const arma::vec& pars)
{
	const int nit = ncat.n_elem;
	
	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);

	ll_pl2_v2 f(itrace, theta, ip, pi, pcni, px, 
			pgroup, inp, icnp, mu, sigma, item, a.col(item));
	
	
	return f(pars);
}			

// [[Rcpp::export]]
arma::vec test_gr_v2(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta, 
			const arma::ivec& ip, const arma::ivec& pi, const arma::ivec& pcni, const arma::ivec& px, 
			const arma::ivec& pgroup, const arma::ivec& inp, const arma::ivec& icnp,
			const arma::vec& mu, const arma::vec& sigma, const int item, const arma::vec& pars)
{
	const int nit = ncat.n_elem;
	
	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);

	ll_pl2_v2 f(itrace, theta, ip, pi, pcni, px, 
			pgroup, inp, icnp, mu, sigma, item, a.col(item));
	
	vec g(pars.n_elem);
	f.df(pars,g);
	return g;
}			


// [[Rcpp::export]]
arma::mat test_hess_v2(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta, 
			const arma::ivec& ip, const arma::ivec& pi, const arma::ivec& pcni, const arma::ivec& px, 
			const arma::ivec& pgroup, const arma::ivec& inp, const arma::ivec& icnp,
			const arma::vec& mu, const arma::vec& sigma, const int item, const arma::vec& pars)
{
	const int nit = ncat.n_elem;
	
	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);

	ll_pl2_v2 f(itrace, theta, ip, pi, pcni, px, 
			pgroup, inp, icnp, mu, sigma, item, a.col(item));
	
	mat h(pars.n_elem,pars.n_elem);
	f.hess(pars,h);
	return h;
}			
