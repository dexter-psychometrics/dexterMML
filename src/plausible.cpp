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


// using SVE
// [[Rcpp::export]]
arma::mat plausible_values_c(const arma::vec& A, const arma::imat& a, const arma::mat& b, const arma::ivec& ncat,
							const arma::ivec& pni, const arma::ivec& pcni, arma::ivec& pi, const arma::ivec& px, const arma::ivec& pop, const arma::ivec& popn, 
							const int npv, const arma::vec& starting_values,
							const int n_prior_updates=10, const int thin=10, const int pgw=80)
{
	const int np = pni.n_elem, nit = A.n_elem, npop=popn.n_elem, max_cat = a.n_rows;
	progress prog(n_prior_updates+(npv-1)*thin+1, pgw,15);
	int tick=0;
	
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
		prog.update(++tick);
		iter++;
	}
	prog.close();
	return out;
}


// simulation
// [[Rcpp::export]]
arma::imat sim_2plc(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat,
					const arma::vec& theta)
{
	const int nit=a.n_cols, np=theta.n_elem,m=a.n_rows;
	mat aA(a.n_rows,nit,fill::zeros);
	for(int i=0; i<nit; i++)
		for(int k=1;k<ncat[i];k++)
			aA.at(k,i) = a.at(k,i)* A[i];
	imat out(np,nit);
	dqrng::xoshiro256plus rng(SEED); 	
	dqrng::uniform_distribution prl_runif(0, 1);
#pragma omp parallel
	{	
		vec P(m);	
		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(omp_get_thread_num() + 1);
		int k;
	#pragma omp for
		for(int i=0;i<nit;i++)
		{		
			P[0]=1;
			for(int p=0;p<np;p++)
			{
				for(k=1; k<ncat[i]; k++)
					P[k] = P[k-1] + std::exp(aA.at(k,i)*(theta[p]-b.at(k,i)));
				const double u=P[k-1]*prl_runif(lrng);
				k=0;
				while(u>P[k]) k++;	
				out.at(p,i)=k;
			}
		}
	}
	return out;
}


