#include <limits>
#include <RcppArmadillo.h>
#include "minimize.h"
#include "pl2_item.h"
#include "shared.h"
#include "posterior.h"

using namespace arma;
using Rcpp::Named;


// without the item prior part
// [[Rcpp::export]]
long double loglikelihood_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, 
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup)
{
	const int nit = ncat.n_elem, nt = theta.n_elem, ng=mu.n_elem;
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);

	return loglikelihood(itrace, pcni, pi, px, theta, mu, sigma, pgroup);
}


template <class T>
vec m_step_2pl(T &f, const double A, const vec& b, const int ncat, const double tol, int& min_error)
{
	vec pars = b.head(ncat);
	pars[0] = A;
	int itr=0,err=0;
	double ll_itm=0;
			
	nlm(pars, tol, itr, ll_itm, f, err);	
			
	if(f.A_prior!=1 && (std::abs(A) < .05 || max(abs(pars.tail(ncat-1))) > 50))
	{
		// 2pl can be poorly identified with local minima
		// on opposite sides of A=0, attempt to break out with a restart of nlm
		int err2=0;
		itr=0;
		double ll_itm2=0;
		vec pars2 = -b.head(ncat)/5;
		pars2[0] = -2*A;
		nlm(pars2, tol, itr, ll_itm2, f, err2);	
		if(err2==0 && ll_itm2<ll_itm)
		{
			min_error=err2;
			return pars2;
		}			
	}
	min_error=err;
	return pars;	
}


void identify_2pl(vec& mu, vec& sigma, const int ref_group, vec& A, mat& b, const int A_prior)
{
	const double mm = mu[ref_group], ss = sigma[ref_group];
	if(A_prior == 0)
	{
		mu = (mu - mm)/ss;
		sigma /= ss;
		A *= ss;
		b = (b - mm)/ss;
	} else
	{
		mu -= mm;
		b -= mm;
	}
	b.row(0).zeros();
}


// [[Rcpp::export]]
Rcpp::List estimate_pl2(arma::imat& a, const arma::vec& A_start, const arma::mat& b_start, const arma::ivec& ncat,
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						const arma::vec& theta_start, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const arma::ivec& item_fixed,
						const arma::ivec ip, const arma::ivec& inp, const arma::ivec& icnp,
						const int ref_group=0, const int A_prior=0, const double A_mu=0, const double A_sigma=0.5, 
						const int use_m2 = 150, const int max_iter=200, const int pgw=80, const int max_pre=10)
{
	const int nit = a.n_cols, nt = theta_start.n_elem, np = pni.n_elem, ng=gn.n_elem;
	
	const bool any_m2 = any(inp < use_m2);
	
	progress_est prog(max_iter, pgw);
	
	mat b = b_start;
	vec A = A_start, sigma = sigma_start, mu=mu_start, theta = theta_start;;

	field<mat> r(nit), itrace(nit);
	for(int i=0; i<nit; i++)
	{
		itrace(i) = pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);
		r(i) = mat(nt,ncat[i]);
	}
	vec thetabar(np,fill::zeros);
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const double tol = 1e-10;
	int iter = 0, min_error=0, stop=0;
	long double ll, old_ll=-std::numeric_limits<long double>::max(), prior_part=0;
	double maxdif_A=0, maxdif_b=0;
	
	vec h_ll(max_iter, fill::zeros);
	
	cube store_b(b.n_rows, nit, 2, fill::zeros);
	mat store_A(nit,2,fill::zeros), store_mu(ng,2,fill::zeros), store_sigma(ng,2,fill::zeros);
	int store_i=0;	
	
	vec pre_ll(max_pre,fill::zeros);
	mat mu_hist(ng,max_iter,fill::zeros), sd_hist(ng,max_iter,fill::zeros);
	

	for(int itr=0; itr<max_pre; itr++)
	{
		estep(itrace, pni, pcni, pi, px, theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);
		pre_ll[itr] = ll;	
			
		for(int g=0;g<ng;g++)
		{		
			mu[g] = sum_theta[g]/gn[g];		
			sigma[g] = std::sqrt(sum_sigma2[g]/gn[g] - mu[g] * mu[g]);
		}
			
		if(ref_group >= 0) identify_2pl(mu, sigma, ref_group, A, b, A_prior);
		
		if(ng>1 || ref_group < 0)
			scale_theta(mu, sigma, gn, theta_start, theta);
			
		for(int i=0; i<nit; i++)
			pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i], itrace(i));		

		if( ll < old_ll + 0.1 ) break;
		old_ll=ll;			
	}

	vec pre_mu=mu,pre_sigma=sigma;
	vec pre_A = A;
	mat pre_b = b;

	old_ll=-std::numeric_limits<long double>::max(); // reset as prior part was not taken into account in pre
	
	for(; iter<max_iter; iter++)
	{
		estep(itrace, pni, pcni, pi, px, 
					theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);
		
		h_ll[iter] = ll + prior_part;
		
		if(ll + prior_part < old_ll)
		{
			stop += 2;
			break;
		}
		old_ll = prior_part+ll;
		maxdif_A=0; maxdif_b=0; prior_part=0;
		store_A.col(store_i)=A;
		store_b.slice(store_i)=b;
		store_mu.col(store_i)=mu;
		store_sigma.col(store_i)=sigma;
		store_i=1-store_i;
		
		
#pragma omp parallel for reduction(max: maxdif_A, maxdif_b) reduction(+:min_error, prior_part)
		for(int i=0; i<nit; i++)
		{	
			if(item_fixed[i] == 1 || inp[i] < use_m2)
				continue;
			int err=0;	
			ll_pl2 f(a.colptr(i), theta.memptr(), r(i), A_prior, A_mu, A_sigma);

			vec pars = m_step_2pl(f, A[i], b.col(i), ncat[i], tol, err);
			
			min_error += err;
			if(err==0)
			{
				maxdif_A = std::max(maxdif_A, std::abs(A[i] - pars[0]));
				A[i] = pars[0];
				for(int k=1;k<ncat[i];k++)
				{
					maxdif_b = std::max(maxdif_b, std::abs(b.at(k,i) - pars[k]));
					b.at(k,i) = pars[k];
				}		
				prior_part -= f.prior_part_ll(A[i]);
			}
		}
		
		for(int g=0;g<ng;g++)
		{		
			mu[g] = sum_theta[g]/gn[g];		
			sigma[g] = std::sqrt(sum_sigma2[g]/gn[g] - mu[g] * mu[g]);
		}
		
		mu_hist.col(iter) = mu;
		sd_hist.col(iter) = sigma;
		
		if(ref_group >= 0) identify_2pl(mu, sigma, ref_group, A, b, A_prior);
		
		if(any_m2)
		{
			for(int i=0; i<nit; i++)
				pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i], itrace(i));

#pragma omp parallel for reduction(max: maxdif_A, maxdif_b) reduction(+:min_error, prior_part)
			for(int i=0; i<nit; i++)
			{	
				if(item_fixed[i] == 1 || inp[i] >= use_m2)
					continue;	
				
				ll_pl2_v2 f(itrace, theta, ip, pi, pcni, px, 
								pgroup, inp, icnp, mu, sigma, i, a.col(i),
								A_prior, A_mu, A_sigma);
				int err=0;
				vec pars = m_step_2pl(f, A[i], b.col(i), ncat[i], tol, err);
				
				min_error += err;
				if(err==0)
				{
					maxdif_A = std::max(maxdif_A, std::abs(A[i] - pars[0]));
					A[i] = pars[0];
					for(int k=1;k<ncat[i];k++)
					{
						maxdif_b = std::max(maxdif_b, std::abs(b.at(k,i) - pars[k]));
						b.at(k,i) = pars[k];
					}		
					prior_part -= f.prior_part_ll(A[i]);
				}
			}
		}
		
		if(ng>1 || ref_group<0)
			scale_theta(mu, sigma, gn, theta_start, theta);
		
		for(int i=0; i<nit; i++)
			pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i],itrace(i));
		
		if(min_error > 0)
		{
			stop += 1;
			break;
		}
		if(iter > 100 && min(abs(A)) < 0.005)
		{
			stop += 8;
			break;
		}	
	
		prog.update(std::max(maxdif_b, maxdif_A), iter);
				
		if(maxdif_b < .0001 && maxdif_A < .0001)
			break;
		
	}
	if(iter>=max_iter-1)
		stop += 4;
		
	prog.close();		
	
	return Rcpp::List::create(Named("A")=A, Named("b")=b, Named("thetabar") = thetabar, Named("mu") = mu, Named("sigma") = sigma, 
							  Named("niter")=iter, Named("theta")=theta, Named("prior_part") = prior_part, 
							  /*Named("r")=r,*/		
		Named("debug")=Rcpp::List::create( 	Named("error")=stop, Named("maxdif_A")=maxdif_A, Named("maxdif_b")=maxdif_b,
											Named("ll_history") = h_ll,	Named("store_A")=store_A, Named("store_b")=store_b,  
											Named("store_mu") = store_mu, Named("store_sigma") = store_sigma, Named("store_i") = 2-store_i,
											Named("pre_ll") = pre_ll,
											Named("pre_mu")=pre_mu,Named("pre_sigma")=pre_sigma,
											Named("pre_A")=pre_A,Named("pre_b")=pre_b,
											Named("mu_hist") = mu_hist, Named("sd_hist")=sd_hist)); 
}






