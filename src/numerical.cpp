#include <stdexcept>
#include <RcppArmadillo.h>
#include "shared.h"
#include "pl2_est.h"
#include "nrm_est.h"

using namespace arma;


// Numerical hessians
// all functions assume: ref_group=0, nothing fixed, no priors
// just for testing


struct parvec_2pl
{
	arma::vec A;
	arma::mat b;
	arma::vec mu;
	arma::vec sigma;
	arma::ivec ncat,cncat;
	int nit, ng;
	
	parvec_2pl(const ivec& ncat_, const int ng_)
	{
		ncat = ncat_;
		nit = ncat.n_elem; 
		ng = ng_;
		cncat = arma::ivec(nit+1);
		cncat[0] = 0;
		std::partial_sum(ncat.begin(),ncat.end(),cncat.begin()+1);
		
		int max_cat = ncat.max();
		A = arma::vec(nit);
		b = arma::mat(max_cat,nit);
		mu = arma::vec(ng);
		sigma = arma::vec(ng);
	}
	void re_fill(const arma::vec& A_in, const arma::mat& b_in, const arma::vec& mu_in, const arma::vec& sigma_in)
	{
		A=A_in;
		b=b_in;
		mu=mu_in;
		sigma=sigma_in;		
	}
	double& operator[](const int q)
	{
		int p=q;
		if(p >= cncat[nit])
		{
			p += 2; //ref group
			p -= cncat[nit];
			if(p % 2 == 0)
			{
				return mu[p/2];
			} else
			{
				return sigma[(p-1)/2];
			}
		}
		else
		{
			for(int i=0; i<nit;i++)
			{
				if(cncat[i]==p)
				{
					return A[i];
				} 
				if(p>cncat[i] && p<cncat[i+1])
				{
					return b.at(p-cncat[i],i);
				}
			}
		}
		throw std::out_of_range ("parvec");
	}
};

struct parvec_nrm
{
	ivec ncat,cncat;
	arma::mat b;
	arma::vec mu;
	arma::vec sigma;
	int nit,ng;
	
	parvec_nrm(const ivec& ncat_, const int ng_)
	{
		ncat = ncat_ - 1;
		nit = ncat.n_elem; 
		ng = ng_;
		cncat = arma::ivec(nit+1);
		cncat[0] = 0;
		std::partial_sum(ncat.begin(),ncat.end(),cncat.begin()+1);
		
		b = arma::mat(1+ncat.max(),nit);
		mu = arma::vec(ng);
		sigma = arma::vec(ng);
	}
	void re_fill(const arma::mat& b_in, const arma::vec& mu_in, const arma::vec& sigma_in)
	{
		b=b_in;
		mu=mu_in;
		sigma=sigma_in;		
	}
	double& operator[](const int q)
	{
		int p=q;
		if(p >= cncat[nit])
		{
			p += 1; //ref group
			p -= cncat[nit];
			if(p % 2 == 0)
			{
				return mu[p/2];
			} else
			{
				return sigma[(p-1)/2];
			}
		}
		else
		{
			for(int i=0; i<nit;i++)
			{
				if(p >= cncat[i])
				{
					return b.at(p-cncat[i]+1,i);
				}
			}
		}
		throw std::out_of_range ("parvec");
	}
};




// [[Rcpp::export]]
arma::mat num_hessian_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, 
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup,
				const int A_prior=0, const double A_mu=0, const double A_sigma=0.5,
				const double ddelta = 1e-5)
{
	const long double delta = ddelta;
	
	const int ng = mu.n_elem, nit=a.n_cols;
	
	int npar = 2*ng-1 + accu(ncat)-nit;

	mat hess(npar,npar,fill::zeros);
	parvec_2pl p(ncat,ng);	
	
	long double s1,s2,s3,s4;
	long double fx = loglikelihood_2pl(a, A, b, ncat, pni, pcni, pi, px, theta, mu, sigma, pgroup);
	for(int i=0; i<npar; i++)
	{
		p.re_fill(A,b,mu,sigma);
        p[i] = p[i] + 2 * delta; 
		s1 = loglikelihood_2pl(a, p.A, p.b, ncat, pni, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
        p[i] = p[i] - 4 * delta; 
		s3 = loglikelihood_2pl(a, p.A, p.b, ncat, pni, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
        hess(i, i) = (s1 - 2*fx + s3) / (4 * SQR(delta));
		for(int j=i+1; j<npar; j++)
		{
			p.re_fill(A,b,mu,sigma);
			p[i] = p[i] + delta; 
			p[j] = p[j] + delta; 
			s1 = loglikelihood_2pl(a, p.A, p.b, ncat, pni, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
			p[j] = p[j] - 2*delta; 
			s2 = loglikelihood_2pl(a, p.A, p.b, ncat, pni, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
			p[i] = p[i] - 2*delta; 
			s4 = loglikelihood_2pl(a, p.A, p.b, ncat, pni, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
			p[j] = p[j] + 2*delta; 
			s3 = loglikelihood_2pl(a, p.A, p.b, ncat, pni, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
			hess(i,j)  = (s1 - s2 - s3 + s4) / (4 * SQR(delta));
			hess(j,i) = hess(i,j);
		}
	}
	return hess;
}



// [[Rcpp::export]]
arma::mat num_hessian_nrm(const arma::imat& a, const arma::mat& b, const arma::ivec& ncat, 
				const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup,
				const double ddelta = 1e-5)
{
	const long double delta = ddelta;
	
	const int ng = mu.n_elem, nit=a.n_cols;
	
	int npar = 2*ng-2;

	for(int i=0; i<nit; i++) 
		npar += ncat[i];	

	mat hess(npar,npar,fill::zeros);
	parvec_nrm p(ncat,ng);	
	
	long double s1,s2,s3,s4;
	long double fx = loglikelihood_nrm(a, b, ncat, pcni, pi, px, theta, mu, sigma, pgroup);
	
	for(int i=0; i<npar; i++)
	{
		p.re_fill(b,mu,sigma);
        p[i] = p[i] + 2 * delta; 
		s1 = loglikelihood_nrm(a, p.b, ncat, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
        p[i] = p[i] - 4 * delta; 
		s3 = loglikelihood_nrm(a, p.b, ncat, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
        hess(i, i) = (s1 - 2*fx + s3) / (4 * SQR(delta));
		for(int j=i+1; j<npar; j++)
		{
			p.re_fill(b,mu,sigma);
			p[i] = p[i] + delta; 
			p[j] = p[j] + delta; 
			s1 = loglikelihood_nrm(a, p.b, ncat, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
			p[j] = p[j] - 2*delta; 
			s2 = loglikelihood_nrm(a, p.b, ncat, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
			p[i] = p[i] - 2*delta; 
			s4 = loglikelihood_nrm(a, p.b, ncat, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
			p[j] = p[j] + 2*delta; 
			s3 = loglikelihood_nrm(a, p.b, ncat, pcni, pi, px, theta, p.mu, p.sigma, pgroup);
			hess(i,j)  = (s1 - s2 - s3 + s4) / (4 * SQR(delta));
			hess(j,i) = hess(i,j);
		}
	}
	return hess;
}
