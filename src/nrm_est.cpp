#include <RcppArmadillo.h>
#include "minimize.h"
#include "nrm_item.h"
#include "shared.h"

using namespace arma;
using Rcpp::Named;


void estep_nrm(imat& a, mat& b, const mat& exp_at, const ivec& ncat, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, field<mat>& r, vec& thetabar, vec& sumtheta, vec& sumsig2, const vec& mu, const vec& sigma, const ivec& pgroup, long double& ll)
{
	const int nit = ncat.n_elem, nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;
	
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
	{
		itrace(i) = nrm_trace(theta, a.col(i), b.col(i), ncat[i], exp_at);
		r(i).zeros();
	}
	
	mat sigma2(nt, ng, fill::zeros);
	sumtheta.zeros();
	
	ll=0;

#pragma omp parallel
	{
		vec posterior(nt);
#pragma omp for reduction(+: r, sigma2, sumtheta, ll)
		for(int p=0; p<np;p++)
		{
			int g = pgroup[p];
			posterior = posterior0.col(g);
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				posterior %= itrace(pi[indx]).col(px[indx]);

			double sp = accu(posterior);
			// LL according to Bock/Aitkin 1981 eq (5) and (6)
			ll += std::log(sp); 
			posterior = posterior / sp;
			sumtheta[g] += thetabar[p] = accu(posterior % theta);
			
			sigma2.col(g) += posterior;
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				r(pi[indx]).col(px[indx]) += posterior;
		}
	}

	for(int g=0; g<ng;g++)
		sumsig2[g] = accu(sigma2.col(g) % square(theta));
}

// [[Rcpp::export]]
long double loglikelihood_nrm(const arma::imat& a, const arma::mat& b, const arma::ivec& ncat, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup)
{
	const int nit = ncat.n_elem, nt = theta.n_elem, np = pcni.n_elem-1, ng = mu.n_elem;
	const int max_a = a.max();
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = nrm_trace(theta, a.col(i), b.col(i), ncat[i], exp_at);
	
	long double ll=0;

#pragma omp parallel
	{
		vec posterior(nt);
#pragma omp for reduction(+: ll)
		for(int p=0; p<np;p++)
		{
			posterior = posterior0.col(pgroup[p]);
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				posterior %= itrace(pi[indx]).col(px[indx]);

			ll += std::log(accu(posterior)); 		
		}
	}
	return ll;
}


// if fixed parameters, set ref_group to a negative nbr
// [[Rcpp::export]]
Rcpp::List estimate_nrm(arma::imat& a, const arma::mat& b_start, const arma::ivec& ncat,
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const arma::ivec& item_fixed,
						const int ref_group=0, const int max_iter = 200, const int pgw=80)
{
	const int nit = a.n_cols, nt = theta.n_elem, np = pni.n_elem, ng=gn.n_elem;	
	const int max_a = a.max();
	
	progress_est prog(max_iter, pgw);
	
	//lookup table
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	mat b = b_start;

	field<mat> r(nit);
	for(int i=0; i<nit; i++)
		r(i) = mat(nt,ncat[i]);
	
	vec thetabar(np,fill::zeros);
	
	vec sigma = sigma_start, mu=mu_start;
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const double tol = 1e-10;
	int iter = 0,min_error=0,stop=0;
	long double ll, old_ll=-std::numeric_limits<long double>::max();
	double maxdif_b;
	
	for(; iter<max_iter; iter++)
	{

		estep_nrm(a, b, exp_at, ncat, pni, pcni, pi, px, 
					theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);
		
		if(ll < old_ll)
		{
			stop += 2;
			break;
		}
		
		maxdif_b=0;
#pragma omp parallel for reduction(max: maxdif_b) reduction(+:min_error)
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
				maxdif_b = std::max(maxdif_b, std::abs(b.at(k,i) - pars[k-1]));
				b.at(k,i) = pars[k-1];
			}
		}
		if(min_error>0)
		{
			stop += 1;
			break;
		}
		for(int g=0;g<ng;g++)
		{			
			if(g==ref_group)
			{
				mu[g] = 0;
				sigma[g] = std::sqrt(sum_sigma2[g]/gn[g]);
			}
			else
			{
				mu[g] = sum_theta[g]/gn[g];		
				sigma[g] = std::sqrt(sum_sigma2[g]/gn[g] - mu[g] * mu[g]);
			}
		}
		
		//printf("iter: % 4i, logl: %.6f,  max b: %.8f\n", iter, (double)ll, maxdif_b);
		//fflush(stdout);

		prog.update(maxdif_b, iter);
		
		if(maxdif_b < .0001)
			break;
		
		old_ll = ll;
		
	}
	if(iter>=max_iter-1)
		stop += 4;
	
	prog.close();
	
	return Rcpp::List::create(Named("a")=a, Named("b")=b, Named("thetabar") = thetabar, Named("mu") = mu, Named("sd") = sigma, 
								Named("r")=r, Named("LL") = ll, Named("niter")=iter,
								Named("err")=stop, Named("maxdif_b")=maxdif_b); 
}

