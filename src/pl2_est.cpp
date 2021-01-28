#include <RcppArmadillo.h>
#include "minimize.h"
#include "pl2_item.h"
#include "shared.h"

using namespace arma;
using Rcpp::Named;

void estep_pl2(field<mat>& itrace, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, field<mat>& r, vec& thetabar, vec& sumtheta, vec& sumsig2, const vec& mu, const vec& sigma, const ivec& pgroup, long double& ll)
{
	const int nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem, nit=r.n_elem;
		
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	for(int i=0; i<nit; i++)
		r(i).zeros();

	
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
double loglikelihood_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, 
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup)
{
	const int nit = ncat.n_elem, nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;	
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);

	long double ll=0;

#pragma omp parallel
	{
		vec posterior(nt);
#pragma omp for reduction(+: ll)
		for(int p=0; p<np;p++)
		{
			int g = pgroup[p];
			posterior = posterior0.col(g);
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				posterior %= itrace(pi[indx]).col(px[indx]);

			ll += std::log(accu(posterior)); 
		}
	}
	return (double)ll;
}


// stop estimation because of decreasing likelihood
// exceedingly rare now that we treat items with small nbr of obs differently
void est_stop(imat& a, const ivec& ncat, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
		 const ivec& pgroup, vec& theta, const arma::ivec& item_fixed, const vec& h_ll, const int iter, const int ref_group, 
		 const int A_prior, const double A_mu, const double A_sigma, const long double old_ll, const int store_i,
		 vec& A, mat& store_A, mat& b, cube& store_b, vec& mu, mat& store_mu, vec& sigma, mat& store_sigma,
		 long double ll)
{
	const int nit = a.n_cols;
	if(iter>3)
	{
		//quadratic approximation based on likelihood
		mat qx = {{1,0,0},{1,1,1},{1,2,4}};
		vec qy = {h_ll[iter-2], h_ll[iter-1], h_ll[iter]};
		vec qs = solve(qx, qy);
		double qpnt = -qs[1]/(2*qs[2]);
		//store_i, 1-store_i, current
		if(qpnt<=1)
		{
			A = store_A.col(1-store_i) * qpnt + (1-qpnt) * store_A.col(store_i);
			b = store_b.slice(1-store_i) * qpnt + (1-qpnt) * store_b.slice(store_i);
			mu = store_mu.col(1-store_i) * qpnt + (1-qpnt) * store_mu.col(store_i);
			sigma = store_sigma.col(1-store_i) * qpnt + (1-qpnt) * store_sigma.col(store_i);
		}
		else
		{
			qpnt-=1;
			A =  A * qpnt + (1-qpnt) * store_A.col(1-store_i);
			b = b * qpnt + (1-qpnt) * store_b.slice(1-store_i);		
			mu = mu * qpnt + (1-qpnt) * store_mu.col(1-store_i);
			sigma = sigma * qpnt + (1-qpnt) * store_sigma.col(1-store_i);	
		}
		if(ref_group>=0) //prevent rounding errors
		{
			mu[ref_group]=0;
			sigma[ref_group]=1;
		}
		//check if anything is won
		//compute likelihood
		double ll_new = loglikelihood_2pl(a, A, b, ncat, pni, pcni, pi, px, 
											theta, mu, sigma, pgroup);
		double prior_part=0;
				
		if(A_prior>0)
			for(int i=0;i<nit;i++)
				if(item_fixed[i] != 1)
				{
					ll_pl2_base f(A_prior, A_mu, A_sigma);
					prior_part -= f.prior_part_ll(A[i]);
				}
				
		if(ll_new + prior_part < old_ll)
		{
			//reject quadratic approx
			A=store_A.col(1-store_i);
			b=store_b.slice(1-store_i);
			mu=store_mu.col(1-store_i);
			sigma=store_sigma.col(1-store_i);
			// to do: get the old ll without prior part?
		} 
		else
		{
			ll=ll_new;
		}
	}
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


// [[Rcpp::export]]
Rcpp::List estimate_pl2(arma::imat& a, const arma::vec& A_start, const arma::mat& b_start, const arma::ivec& ncat,
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const arma::ivec& item_fixed,
						const arma::ivec ip, const arma::ivec& inp, const arma::ivec& icnp,
						const int ref_group=0, const int A_prior=0, const double A_mu=0, const double A_sigma=0.5, 
						const int use_m2 = 150, const int max_iter=200, const int pgw=80)
{
	const int nit = a.n_cols, nt = theta.n_elem, np = pni.n_elem, ng=gn.n_elem;
	
	progress_est prog(max_iter, pgw);
	
	mat b = b_start;
	vec A = A_start;

	field<mat> r(nit), itrace(nit);
	for(int i=0; i<nit; i++)
	{
		itrace(i) = pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);
		r(i) = mat(nt,ncat[i]);
	}
	vec thetabar(np,fill::zeros);
	
	vec sigma = sigma_start, mu=mu_start;
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const double tol = 1e-10;
	int iter = 0, min_error=0, stop=0;
	long double ll, old_ll=-std::numeric_limits<long double>::max(), prior_part=0;
	double maxdif_A, maxdif_b;
	vec h_ll(max_iter, fill::zeros);
	
	cube store_b(b.n_rows, nit, 2, fill::zeros);
	mat store_A(nit,2,fill::zeros), store_mu(ng,2,fill::zeros), store_sigma(ng,2,fill::zeros);
	int store_i=0;
	
	mat tmp1(2,max_iter,fill::zeros),tmp2(2,max_iter,fill::zeros);
	
	for(; iter<max_iter; iter++)
	{
		estep_pl2(itrace, pni, pcni, pi, px, 
					theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);
		
		h_ll[iter] = ll + prior_part;
		
		if(ll + prior_part < old_ll)
		{
			est_stop(a, ncat, pni, pcni, pi, px, pgroup, theta, item_fixed, h_ll, iter, ref_group, 
					 A_prior, A_mu, A_sigma, old_ll, store_i,
					A, store_A, b, store_b, mu, store_mu, sigma, store_sigma, ll);
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
			if(g==ref_group)
			{
				mu[g] = 0;
				sigma[g] = 1;
			}
			else
			{
				mu[g] = sum_theta[g]/gn[g];		
				sigma[g] = std::sqrt(sum_sigma2[g]/gn[g] - mu[g] * mu[g]);
			}
		}		
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
		for(int i=0; i<nit; i++)
			if(inp[i] < use_m2 && item_fixed[i] != 1)
				pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i],itrace(i));
		
		if(min_error>0)
		{
			stop += 1;
			break;
		}
		if(min(abs(A)) < 1e-4)
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
	
	return Rcpp::List::create(Named("A")=A, Named("b")=b, Named("thetabar") = thetabar, Named("mu") = mu, Named("sd") = sigma, 
									Named("r")=r, Named("LL") = (double)ll, Named("niter")=iter, Named("prior_part") = (double)prior_part,
									Named("err")=stop, Named("maxdif_A")=maxdif_A,Named("maxdif_b")=maxdif_b,
									Named("ll_history") = h_ll,
									Named("store_A")=store_A, Named("store_b")=store_b); 
}

