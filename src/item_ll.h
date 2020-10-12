#ifndef DXM_ITEM_LIKELIHOOD_
#define DXM_ITEM_LIKELIHOOD_

#include <RcppArmadillo.h>
#include "shared.h"

// contains negative likelihood, gradient and hessian for single dichotomous 2pl item given r0 and r1
// for use in minimizing functions
struct ll_2pl_dich
{
	arma::vec r0,r1,theta;
	int n;

	ll_2pl_dich(double* r1p, double* r0p, double* thetap, const int ni)
	{
		n=ni;
		r1 = arma::vec(r1p,n,false,true);		
		r0 = arma::vec(r0p,n,false,true);
		theta = arma::vec(thetap,n,false,true);	
	}
	
	//returns minus log likelihood
	double operator()(const arma::vec& ab)
	{
		double ll=0;
		const double a = ab[0], b = ab[1];
		for(int i=0;i<n;i++)
		{
			double p = 1/(1+std::exp(-a*(theta[i]-b)));
			ll -= r1[i] * std::log(p) + r0[i] * std::log(1-p);
		}
		if(std::isinf(ll))
		{
			printf("infinity, a: %f, b: %f", a, b);
			fflush(stdout);
			Rcpp::stop("inf ll");
		}
		return ll;	
	}
	
	//returns gradient of minus ll
	void df(const arma::vec& ab, arma::vec& g)
	{	
		g.zeros();
		const double a = ab[0], b = ab[1];
		for(int i=0;i<n;i++)
		{
			double e = std::exp(a*(b-theta[i]));
			g[0] -= (b-theta[i]) * (r0[i] - r1[i]*e)/(e+1);
			g[1] -= a * (r0[i]-r1[i]*e)/(e+1);
		}
		if(!g.is_finite())
		{
			printf("infinity in gradient, a: %f, b: %f", a, b);
			fflush(stdout);
			Rcpp::stop("inf gradient");
		}
	}
	
	void hess(const arma::vec& ab, arma::mat& h)
	{
		h.zeros();
		const double a = ab[0], b = ab[1];
		for(int i=0;i<n;i++)
		{
			double e = std::exp(a*(b-theta[i])), t=theta[i];
			h.at(0,0) += (r0[i]*(e+1) - r0[i] - r1[i]*(e+1)*e - (2*r0[i]-r1[i]*e)*e) \
					* SQR(b-t)/SQR(e+1);
			
			h.at(0,1) += (a*r0[i]*(t-b) - a*(b-t)*(2*r0[i]-r1[i]*e)*e \
							+ r0[i]*(a*(b-t)+1)*(e+1) \
							- r1[i]*(a*(b-t) + 1)*(e+1)*e) \
							/ SQR(e+1);			
			
			e = std::exp(a*(b+t));
			
			h.at(1,1) -= SQR(a)*(r0[i]+r1[i])*e/(std::exp(2*a*b) + std::exp(2*a*t) + 2*e);
		}
		h.at(1,0) = h.at(0,1);
	}
};

// polytomous nrm (can also be used for dichotmous nrm or rasch)
// this uses Timo's b parametrisation
struct ll_nrm
{
	arma::mat r; //rows for categories, columns for theta
	arma::mat exp_theta_a;
	
	int nt;
	int ncat;
	
	// a should include 0 cat
	ll_nrm(const arma::ivec& a, const arma::vec& theta, double* rp) 
	{
		ncat = a.n_elem;
		nt = theta.n_elem;
		r = arma::mat(rp,ncat,nt);
		exp_theta_a = arma::mat(ncat,nt);

		for(int t=0;t<nt;t++)
			for(int k=0;k<ncat;k++)
				exp_theta_a.at(k,t) = std::exp(a[k] * theta[t]);		
	}
	
	double operator()(const arma::vec& b)
	{
		double ll=0;
		for(int t=0;t<nt;t++)
			ll -= accu(r.col(t) % arma::log(exp_theta_a.col(t) % b/accu(exp_theta_a.col(t) % b)));

		return ll;
	}

	void df(const arma::vec& b, arma::vec& g)
	{
		g.zeros();
		for(int t=0; t<nt; t++)
		{
			double s = accu(exp_theta_a.col(t) % b);
			g += r.col(t) % exp_theta_a.col(t)/s + 1/b;
		}
	}
	
	void hess(const arma::vec& b, arma::mat& h)
	{
		h.zeros();
		for(int t=0; t<nt; t++)
		{
			double s = SQR(accu(exp_theta_a.col(t) % b));
			for(int i=0; i<ncat;i++)
			{
				h.at(i,i) += SQR(exp_theta_a.at(i,t))/s;
				for(int j=i+1;j<ncat; j++)
					h.at(i,j) += exp_theta_a.at(i,t) * exp_theta_a.at(j,t)/s;
			}
		}
		for(int i=0;i<ncat;i++)
			for(int j=i+1; j<ncat; j++)
				h.at(j,i) = h.at(i,j);	
	}

};

#endif
