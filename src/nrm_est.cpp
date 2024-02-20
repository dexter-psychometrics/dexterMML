#include <limits>
#include <RcppArmadillo.h>
#include "minimize.h"
#include "nrm_item.h"
#include "shared.h"
#include "posterior.h"

using namespace arma;
using Rcpp::Named;


// [[Rcpp::export]]
long double loglikelihood_nrm(const arma::imat& a, const arma::mat& b, const arma::ivec& ncat, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup)
{
	const int nit = ncat.n_elem, nt = theta.n_elem;
	const int max_a = a.max();
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = nrm_trace(theta, a.col(i), b.col(i), ncat[i], exp_at);
		
	return loglikelihood(itrace, pcni, pi, px, theta, mu, sigma, pgroup);	

}


// if fixed parameters, set ref_group to a negative nbr

void identify_1pl(vec& mu, const int ref_group, mat& b)
{
	const double mm = mu[ref_group];
	
	mu -= mm;
	b += mm;
	b.row(0).zeros();
}


// [[Rcpp::export]]
Rcpp::List estimate_nrm(arma::imat& a, const arma::mat& b_start, const arma::ivec& ncat,
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta_start, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const arma::ivec& item_fixed,
						const int ref_group=0, const int max_iter = 200, const int pgw=80)
{
	const int nit = a.n_cols, nt = theta_start.n_elem, np = pni.n_elem, ng=gn.n_elem;	
	const int max_a = a.max();
	
	vec sigma = sigma_start, mu=mu_start, theta=theta_start;
	
	progress_est prog(max_iter, pgw);
	
	//lookup table
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	field<mat> itrace(nit);
	
	mat b = b_start, old_b=b_start;

	field<mat> r(nit);
	for(int i=0; i<nit; i++)
		r(i) = mat(nt,ncat[i]);
	
	vec thetabar(np,fill::zeros);
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const double tol = 1e-10;
	int iter = 0,min_error=0,stop=0;
	long double ll, old_ll=-std::numeric_limits<long double>::max();
	double maxdif_b=0;
	
	bool adapt_theta = ng > 1 || ref_group < 0;

	for(; iter<max_iter; iter++)
	{
		for(int i=0; i<nit; i++)
			itrace(i) = nrm_trace(theta, a.col(i), b.col(i), ncat[i], exp_at);
		
		estep(itrace, pni, pcni, pi, px, theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);

		if(ll < old_ll)
		{
			if(adapt_theta) adapt_theta = false;
			else
			{
				stop += 2;
				break;
			}	
		}
		
		old_b=b;
#pragma omp parallel for reduction(+:min_error)
		for(int i=0; i<nit; i++)
		{	
			if(item_fixed[i] == 1)
				continue;
			ll_nrm f(a.colptr(i), exp_at, r(i));
			vec pars(b.colptr(i)+1,ncat[i]-1);
			int itr=0,err=0;
			double ll_itm=0;
			if(ncat[i] == 2)
				D1min(pars, tol, itr, ll_itm, f, err); // 1 dimensional minimization
			else
				nlm(pars, tol, itr, ll_itm, f, err);	

			min_error += err;
			for(int k=1;k<ncat[i];k++)
			{
				b.at(k,i) = pars[k-1];
			}
		}
		
		for(int g=0;g<ng;g++)
		{			
			mu[g] = sum_theta[g]/gn[g];		
			sigma[g] = std::sqrt(sum_sigma2[g]/gn[g] - mu[g] * mu[g]);
		}
			
		if(ref_group >= 0) identify_1pl(mu, ref_group, b);
		
		if(!mu.is_finite() || !sigma.is_finite() || sigma.min() < 1e-8)
		{
			sigma.elem(find_nonfinite(sigma)).zeros();
			stop += 16;
			break;
		}

		if(min_error>0)
		{
			stop += 1;
			break;
		}

		maxdif_b = abs(b-old_b).max();
		prog.update(maxdif_b, iter);		

		if(maxdif_b < 0.0001)
			break;
		
		if(adapt_theta)
		{
			scale_theta(mu, sigma, gn, theta_start, theta);
			for(int t=0; t< nt; t++)
				for(int k=1; k<=max_a;k++)
					exp_at.at(k,t) = std::exp(k*theta[t]);		
		}		
		
		old_ll = ll;
		
	}
	if(iter>=max_iter-1)
		stop += 4;
	
	prog.close();
	
	return Rcpp::List::create(Named("a")=a, Named("b")=b, Named("thetabar") = thetabar, Named("mu") = mu, Named("sigma") = sigma, 
								Named("niter")=iter, Named("err")=stop, Named("theta")=theta,
		Named("debug")=Rcpp::List::create( Named("error")=stop, Named("maxdif_b")=maxdif_b,Named("ll")=ll)); 
}



