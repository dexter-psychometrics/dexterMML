#include <RcppArmadillo.h>
#include <xoshiro.h>
#include <dqrng_distribution.h>
#include <omp.h>
#include "shared.h"
using namespace arma;
using Rcpp::Named;
#define SEED std::round(R::runif(0,1) * 2147483647)

/*
# Given samples of plausible values from an hierarchical model with >=2 groups
# this function samples means and variances of plausible values from their posterior
# It updates the prior used in sampling plausible values. 
#
# Based on code provided by Gelman, A., & Hill, J. (2006). Data analysis using regression and 
# multilevel/hierarchical models. Cambridge university press.
#
*/

void sample_musig(const arma::ivec& popn, std::vector<long double>& mu, 
						long double& sigma, long double& sigma_old, double& mu_all, double& sigma_all)
{
	const int J = popn.n_elem;
	const int n = accu(popn);
	
	if(J>1)
	{
		double smu=0;
		for(int j=0; j<J; j++)
		{
			double v_aj = 1/(popn[j]/SQR(sigma_old) + 1/SQR(sigma_all));
			double ahat_j = v_aj * ((popn[j]/SQR(sigma_old)) * mu[j] + (1/SQR(sigma_all))*mu_all);
			mu[j] = R::rnorm(ahat_j, std::sqrt(v_aj));	
			smu+=mu[j];
		}
		mu_all = R::rnorm(smu/n, sigma_all/std::sqrt((double)J));
		sigma = std::sqrt(n*SQR(sigma)/R::rchisq((double)n-1));
		sigma_all=0;
		for(int j=0; j<J; j++) sigma_all += SQR(mu[j]-mu_all);
		sigma_all = std::sqrt(sigma_all/R::rchisq((double)J -1));
	}
	else
	{
		sigma = std::sqrt(1/R::rgamma((n-1.0)/2, 1/(((n-1.0)/2)*(double)SQR(sigma))));
		mu[0] = R::rnorm(mu[0],sigma/std::sqrt((double)n));	
	}
	sigma_old = sigma;	
}

// this is an adaptation of the dexter algorithm
// works fine, however is unacceptably slow for adaptive or random tests
// [[Rcpp::export]]
arma::mat plausible_values_c(const arma::ivec& booklet_id, const arma::ivec& pop, const arma::ivec& pbn, const arma::ivec& pbcn, 
							const arma::ivec& pbnp, const arma::ivec& pbcnp,
							const arma::ivec& scoretab, 
							const arma::ivec& popn, 
							const arma::ivec& dsg_item_id, const arma::ivec& bnit, const arma::ivec& bcnit,
							const arma::vec& A, const arma::imat& a, const arma::mat& b, const arma::ivec& ncat,
							const int npv)
{
	const int nbp = pbn.n_elem, max_cat = a.n_rows, npop=popn.n_elem, nit=a.n_cols;
	const int np = accu(popn);
	const int max_thr = omp_get_max_threads();

	std::vector<long double> mu(npop, 0);
	long double sigma=4,sigma_old=4;
	double mu_all=0, sigma_all=1;
	
	mat out(np,npv);
	
	dqrng::xoshiro256plus rng(SEED); 	
	dqrng::uniform_distribution prl_runif(0, 1);
	
	mat aA(a.n_rows,nit,fill::zeros);
	for(int i=0; i<nit; i++)
		for(int k=1;k<ncat[i];k++)
			aA.at(k,i) = a.at(k,i)* A[i];
	
	const ivec cscoretab = cumsum(scoretab);
	
	const int n_prior_updates = 10;
	int pvcol = 0; 
	
	for(int iter=0; iter < npv + n_prior_updates; iter++)
	{
		if(iter>n_prior_updates) pvcol++;
#pragma omp parallel
		{
			dqrng::xoshiro256plus lrng(rng);      		
			lrng.long_jump(omp_get_thread_num() + 1);
			vec p(max_cat);
			p[0]=1;
			
#pragma omp for	
			for(int bpop=0; bpop<nbp; bpop++)
			{
				const int bk = booklet_id[bpop];
				const int popnr = pop[bpop];
				int np = pbnp[bpop], k;
				ivec stb = scoretab.subvec(pbcn[bpop], pbcn[bpop+1]-1); // check				

				dqrng::normal_distribution prl_rnorm((double)mu[popnr], (double)sigma);
				
				// this is fine for people organized in booklets, but especially for n=1 booklets we should use a version with min/max
				while(np>0)
				{
					//const double theta = R::rnorm((double)mu[popnr], (double)sigma);
					const double theta = prl_rnorm(lrng);
					int x=0;
					for(int ii=bcnit[bk]; ii< bcnit[bk+1]; ii++)
					{
						const int i = dsg_item_id[ii];
						for(k=1; k<ncat[i]; k++)
							p[k] = p[k-1] + std::exp(aA.at(k,i)*(theta-b.at(k,i)));

						//const double u=p[k-1]*R::runif(0,1);
						const double u=p[k-1]*prl_runif(lrng);
						k=0;
						while(u>p[k]) k++;		
				
						x += a.at(k,i);						
					}
					
					//check for interruptible omp with progress bar
					//https://cran.r-project.org/web/packages/RcppProgress/index.html
					if(stb[x] > 0)
					{
						out(cscoretab[pbcn[bpop] + x] - stb[x], pvcol) = theta;
						stb[x]--;
						np--;
					}					
				}				
			}
		}
		// compute mu, sigma
		std::fill(mu.begin(), mu.end(),0);
		for(int bpop=0; bpop<nbp; bpop++)
		{
			const int popnr = pop[bpop]; 
			for(int prs = pbcnp[bpop]; prs < pbcnp[bpop+1]; prs++)
				mu[popnr] += out.at(prs,pvcol);
		}

		for(int ps=0;ps<npop;ps++)
			mu[ps] /= popn[ps];
		sigma=0;
		for(int bpop=0; bpop<nbp; bpop++)
		{
			const int popnr = pop[bpop];
			for(int prs = pbcnp[bpop]; prs < pbcnp[bpop+1]; prs++)
				sigma += SQR(out.at(prs,pvcol)-mu[popnr]);
		}
		sigma = std::sqrt(sigma/np);

		
		sample_musig(popn, mu, sigma, sigma_old, mu_all, sigma_all);		
		rng.long_jump(max_thr+1);
	}
	return out;
}

// using SVE
// [[Rcpp::export]]
arma::mat plausible_values_c2(const arma::vec& A, const arma::imat& a, const arma::mat& b, const arma::ivec& ncat,
							const arma::ivec& pni, const arma::ivec& pcni, arma::ivec& pi, const arma::ivec& px, const arma::ivec& pop, const arma::ivec& popn, 
							const int npv, const arma::vec& starting_values,
							const int n_prior_updates=10, const int thin=10)
{
	const int np = pni.n_elem, nit = A.n_elem, npop=popn.n_elem, max_cat = a.n_rows;
	
	dqrng::xoshiro256plus rng(SEED); 	
	dqrng::uniform_distribution prl_runif(0, 1);
	dqrng::normal_distribution prl_rnorm(0, 1);
	
	const int max_thr = omp_get_max_threads();
	
	std::vector<long double> mu(npop, 0);
	long double sigma=4,sigma_old=4;
	double mu_all=0, sigma_all=1;
		
	mat out(np,npv);
		
	out.col(0) = starting_values;
	mat aA(a.n_rows,nit,fill::zeros);
	for(int i=0; i<nit; i++)
		for(int k=1;k<ncat[i];k++)
			aA.at(k,i) = a.at(k,i)* A[i];
	
	vec ws(np,fill::zeros);
	
	for(int p=0;p<np;p++)
		for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
			ws[p] += aA.at(px[indx],pi[indx]);
	
	
	int pvcol=0, prevcol;
	for(int iter = -n_prior_updates;;iter++)
	{	
		prevcol=pvcol;
		if(iter>=0 && iter % thin ==0)
			pvcol++;
		if(pvcol == npv)
			break;
			
		std::fill(mu.begin(), mu.end(),0);
		for(int p=0;p<np;p++)
			mu[pop[p]] += out.at(p,pvcol);
	
		for(int ps=0;ps<npop;ps++)
			mu[ps] /= popn[ps];
		sigma=0;
		for(int p=0;p<np;p++)
			sigma += SQR(out.at(p,pvcol) - mu[pop[p]]);
		sigma = std::sqrt(sigma/np);
		
#pragma omp parallel
		{
			// check if repeated declaration of this does not cause problems
			// maybe we need to jump/use rng after parallel loop?			
			dqrng::xoshiro256plus lrng(rng);      		
			lrng.long_jump(omp_get_thread_num() + 1);			
			vec P(max_cat);
			P[0]=1;
			int k;
#pragma omp for			
			for(int p=0;p<np;p++)
			{		
				//double theta = R::rnorm((double)mu[pop[p]], (double)sigma);
				double theta = prl_rnorm(lrng) * sigma + mu[pop[p]];
				double x=0;
				for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				{
					const int i = pi[indx];
					for(k=1; k<ncat[i]; k++)
						P[k] = P[k-1] + std::exp(aA.at(k,i)*(theta-b.at(k,i)));

					const double u=P[k-1]*prl_runif(lrng);
					k=0;
					while(u>P[k]) k++;		
					x += aA.at(k,i);	
				}
				double acc = std::exp((theta-out.at(p, prevcol))*(ws[p]-x));
				if(prl_runif(lrng) < acc)
					out.at(p,pvcol)=theta;
				else
					out.at(p,pvcol)=out.at(p,prevcol);					
			}
		}
		
		
		sample_musig(popn, mu, sigma, sigma_old, mu_all, sigma_all);
		rng.long_jump(max_thr+1);
		iter++;
	}
	return out;
}
