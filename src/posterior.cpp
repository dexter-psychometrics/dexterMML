
#include <RcppArmadillo.h>
#include "shared.h"

using namespace arma;


void estep(field<mat>& itrace, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
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


// normalized posterior for persons
mat normalized_posterior(field<mat>& itrace, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
						const vec& theta, const vec& mu, const vec& sigma, const ivec& pgroup)
{
	const int nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;
		
	mat posterior0(nt,ng), posterior(nt,np);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);	

#pragma omp parallel for
	for(int p=0; p<np;p++)
	{
		posterior.col(p) = posterior0.col(pgroup[p]);
		
		for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
			posterior.col(p) %= itrace(pi[indx]).col(px[indx]);
		long double s=0;
		for(int t=0;t<nt;t++)
			s+= posterior.at(t,p);
		posterior.col(p) /= s;
		
	}
	return posterior;
}


// general likelhood
// without prior parts
long double loglikelihood(field<mat>& itrace, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup)
{
	const int nt = theta.n_elem, np = pcni.n_elem-1, ng = mu.n_elem;
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);
	
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

