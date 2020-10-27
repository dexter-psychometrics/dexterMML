
#include <RcppArmadillo.h>
#include "minimize.h"
#include "d2pl_item.h"
#include "shared.h"

using namespace arma;
using Rcpp::Named;



// [[Rcpp::export]]
double test_ll_d2(arma::mat& r, arma::vec& theta, const arma::vec& par)
{
	ll_2pl_dich f(r.colptr(1), r.colptr(0), theta.memptr(), theta.n_elem);
	
	return f(par);
} 

// [[Rcpp::export]]
arma::vec test_gradient_d2(arma::mat& r, arma::vec& theta, const arma::vec& par)
{
	ll_2pl_dich f(r.colptr(1), r.colptr(0), theta.memptr(), theta.n_elem);
	vec g(par.n_elem);
	f.df(par,g);
	return g;
}

// [[Rcpp::export]]
arma::mat test_hess_d2(arma::mat& r, arma::vec& theta, const arma::vec& par)
{
	ll_2pl_dich f(r.colptr(1), r.colptr(0), theta.memptr(), theta.n_elem);
	mat h(par.n_elem,par.n_elem);
	f.hess(par,h);
	return h;
}



// this is a reduced version of the estep, can be use to compute the final ll
long double LL_2pl_dich(const vec& a, const vec& b, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, const vec& mu, const vec& sigma, const ivec& pgroup)

{
	const int nit = a.n_elem, nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;
	mat itrace(nt,nit);
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);
	
	for(int i=0; i<nit; i++)
		itrace.col(i) = 1/(1+exp(-a[i]*(theta-b[i])));
	
	long double ll=0;
	
#pragma omp parallel
	{
		vec posterior(nt);
# pragma omp for reduction(+:ll)
		for(int p=0; p<np;p++)
		{
			int g = pgroup[p];
			posterior = posterior0.col(g);
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
			{
				if(px[indx] == 1)
					posterior %= itrace.col(pi[indx]);
				else
					posterior %= 1-itrace.col(pi[indx]);
			}	
			ll += std::log(accu(posterior));
		}
	}
	return ll;
}


void estep_2pl_dich(const vec& a, const vec& b, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, mat& r0, mat& r1, vec& thetabar, vec& sumtheta, vec& sumsig2, const vec& mu, const vec& sigma, const ivec& pgroup, double& ll)
{
	const int nit = a.n_elem, nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;
	mat itrace(nt,nit);
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	
	r0.zeros();
	r1.zeros();
	sumtheta.zeros();
	
	for(int i=0; i<nit; i++)
		itrace.col(i) = 1/(1+exp(-a[i]*(theta-b[i])));
		
	mat sigma2(nt, ng, fill::zeros);
	
	ll=0;
	
#pragma omp parallel
	{
		vec posterior(nt);
# pragma omp for reduction(+:r0,r1,sigma2, sumtheta,ll)
		for(int p=0; p<np;p++)
		{
			int g = pgroup[p];
			posterior = posterior0.col(g);
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
			{
				if(px[indx] == 1)
					posterior %= itrace.col(pi[indx]);
				else
					posterior %= 1-itrace.col(pi[indx]);
			}	
			double sp = accu(posterior);
			// LL according to Bock/Aitkin 1981 eq (5) and (6), omitting constant C casue I don't know what C is
			ll += std::log(sp); 
			posterior = posterior / sp;
			sumtheta[g] += thetabar[p] = accu(posterior % theta);
			
			sigma2.col(g) += posterior;
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
			{
				if(px[indx] == 1)
					r1.col(pi[indx]) += posterior;
				else
					r0.col(pi[indx]) += posterior;
			}		
		}
	}

	for(int g=0; g<ng;g++)
		sumsig2[g] = accu(sigma2.col(g) % square(theta));
}



// [[Rcpp::export]]
Rcpp::List estimate_2pl_dich_multigroup(const arma::vec& a_start, const arma::vec& b_start, 
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const int ref_group=0, const int max_iter=200)
{
	const int nit = a_start.n_elem, nt = theta.n_elem, np = pni.n_elem, ng=gn.n_elem;;
	
	vec a(a_start.memptr(),nit), b(b_start.memptr(),nit);
	
	mat r0(nt,nit, fill::zeros), r1(nt,nit, fill::zeros);	
	
	vec thetabar(np,fill::zeros);
	
	vec sigma = sigma_start, mu=mu_start;
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const double tol = 1e-10;
	int iter = 0;
	double ll;
	
	for(; iter<max_iter; iter++)
	{
		estep_2pl_dich(a, b, pni, pcni, pi, px, 
						theta, r0, r1, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);
		
		
		double maxdif_a=0, maxdif_b=0;
#pragma omp parallel for reduction(max: maxdif_a, maxdif_b)
		for(int i=0; i<nit; i++)
		{				
			ll_2pl_dich f(r1.colptr(i), r0.colptr(i), theta.memptr(), nt);
			vec pars(2);
			pars[0] = a[i];
			pars[1] = b[i];
			int itr=0;
			double ll_itm=0;
			
			// minimize, still need to tweak lnsrch a bit of replace by better line search algo
			dfpmin(pars, tol, itr, ll_itm, f);
			
			maxdif_a = std::max(maxdif_a, std::abs(a[i]-pars[0]));
			maxdif_b = std::max(maxdif_b, std::abs(b[i]-pars[1]));
			
			a[i] = pars[0];
			b[i] = pars[1];
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
		
		//printf("\r% 3i", iter);
		printf("iter: % 4i, logl: %.6f, max a: %.8f, max b: %.8f\n", iter, ll, maxdif_a, maxdif_b);
		fflush(stdout);
		
		
		if(maxdif_a < .0001 && maxdif_b < .0001)
			break;
		
	}
	
	printf("\n");
	fflush(stdout);
	
	return Rcpp::List::create(Named("a")=a, Named("b")=b, Named("thetabar") = thetabar, Named("mu") = mu, Named("sd") = sigma, 
									Named("LL") = ll, Named("niter")=iter,Named("r0")=r0, Named("r1")=r1 ); 
}




