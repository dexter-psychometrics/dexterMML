
#include <RcppArmadillo.h>
#include "minimize.h"
#include "item_ll.h"

using namespace arma;
using Rcpp::Named;

// += for fields with all equal dimensons for use in omp reduction
void field_plus(field<mat>& a, const field<mat>& b) 
{
    const int n = b.n_slices;

	for(int i=0; i<n; i++)
		a(i) += b(i);
}

#pragma omp declare reduction( + : arma::field<mat> : field_plus(omp_out, omp_in)) \
initializer( omp_priv = omp_orig )

#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) \
initializer( omp_priv = omp_orig )

#pragma omp declare reduction( + : arma::vec : omp_out += omp_in ) \
initializer( omp_priv = omp_orig )





//a, b matrix, maxcat rows, nit columns


mat nrm_trace(const arma::vec& theta, int* ap, double* bp, const int ncat)
{
	const int nt = theta.n_elem;
	mat out(nt,ncat);		
	ivec a(ap, ncat, false, true);
	vec b(bp, ncat, false, true);
	
	for(int j=0; j< nt; j++)
	{
		double sm=0;
		for(int k=0;k<ncat;k++)
			sm += b[k]*std::exp(theta[j]*a[k]);
		for(int k=0;k<ncat;k++)
			out.at(j,k) = b[k]*std::exp(theta[j]*a[k])/sm;
	}
		
	return out;
}

// kan nog sneller gemaakt omdar er sufstats zijn voor deze 
// x moet gehercodeerd zijn naar categorie nummers, zonder gaten
void estep_nrm(imat& a, mat& b, const ivec& ncat, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, field<mat>& r, vec& thetabar, vec& sumtheta, vec& sumsig2, const vec& mu, const vec& sigma, const ivec& pgroup, double& ll)
{
	const int nit = ncat.n_elem, nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;
	
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
	{
		itrace(i) = nrm_trace(theta, a.colptr(i), b.colptr(i), ncat[i]);
		r(i).zeros();
	}
	
	mat sigma2(nt, ng, fill::zeros);
	sumtheta.zeros();
	
	ll=0;
	
/* pragma omp parallel */
	{
		vec posterior(nt);
/*  pragma omp for reduction(+:r,sigma2, sumtheta,ll) */
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


//a,b ncol items
// [[Rcpp::export]]
Rcpp::List estimate_nrm(arma::imat& a, const arma::mat& b_start, const arma::ivec& ncat,
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, const int ref_group=0)
{
	const int nit = a.n_cols, nt = theta.n_elem, np = pni.n_elem, ng=gn.n_elem;
	
	mat b = b_start;

	field<mat> r(nit);
	for(int i=0; i<nit; i++)
		r(i) = mat(nt,ncat[i]);
	
	vec thetabar(np,fill::zeros);
	
	vec sigma = sigma_start, mu=mu_start;
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const int max_iter = 100;
	const double tol = 1e-8;
	int iter = 0;
	double ll;
	
	for(; iter<max_iter; iter++)
	{

		estep_nrm(a, b, ncat, pni, pcni, pi, px, 
					theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);
		
		
		double maxdif_b=0;
/*  omp parallel for reduction(max: maxdif_b) */
		for(int i=0; i<nit; i++)
		{	
			ll_nrm f(a.col(i), theta, r(i));
			vec pars(b.colptr(i)+1,ncat[i]-1);
			pars=log(pars);
			int itr=0;
			double ll_itm=0;

			dfpmin(pars, tol, itr, ll_itm, f);

			pars = exp(pars);
			for(int k=1;k<ncat[i];k++)
			{
				maxdif_b = std::max(maxdif_b, std::abs(b.at(k,i) - pars[k-1]));
				b.at(k,i) = pars[k-1];
			}
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
		
		//printf("\r% 3i", iter);
		printf("iter: % 4i, logl: %.6f,  max b: %.8f\n", iter, ll, maxdif_b);
		fflush(stdout);
		
		
		if(maxdif_b < .0001)
			break;
		
	}
	
	printf("\n");
	fflush(stdout);
	
	return Rcpp::List::create(Named("a")=a, Named("b")=b, Named("thetabar") = thetabar, Named("mu") = mu, Named("sd") = sigma, 
									Named("LL") = ll, Named("niter")=iter,Named("r")=r); 
}

