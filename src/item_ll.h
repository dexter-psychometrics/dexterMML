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
	arma::mat r;
	arma::mat eat;
	arma::vec p;
	arma::mat exp_at;
	arma::ivec a;
	
	int nt;
	int ncat;
	
	// a includes 0 cat
	ll_nrm(int* a_ptr, arma::mat& exp_at_, arma::mat& r_) 
	{
		ncat = r_.n_cols;
		nt = exp_at_.n_cols;
		p = arma::vec(ncat);
		a = arma::ivec(a_ptr, ncat,false,true);
		r = arma::mat(r_.memptr(),nt,ncat,false,true);
		exp_at = arma::mat(exp_at_.memptr(),exp_at_.n_rows, exp_at_.n_cols,false,true);
	}
	
	// b should not include 0 score
	double operator()(const arma::vec& beta)
	{
		double ll=0;
		const arma::vec b = arma::exp(beta);
		
		for(int t=0; t<nt; t++)
		{
			double s=1;
			p[0] = 1;
			for(int k=1; k<ncat; k++)
			{
				p[k] = b[k-1] * exp_at.at(a[k],t);
				s += p[k];
			}

			for(int k=0; k<ncat; k++)
				ll -= r.at(t,k) * std::log(p[k]/s);
		}		
		
		if(std::isinf(ll))
		{
			b.print("b:");
			fflush(stdout);
			Rcpp::stop("inf ll");
		}

		return ll;
	}

	void df(const arma::vec& beta, arma::vec& g)
	{
		
		g.zeros();
		const arma::vec b = arma::exp(beta);
		
		for(int t=0; t<nt; t++)
		{
			double s=1,ss = r.at(t,0);
			for(int k=1; k<ncat;k++)
			{
				s += b[k-1] * exp_at.at(a[k],t);
				ss += r.at(t,k);
			}
			for(int k=1; k<ncat;k++)
			{
				/* 
					df towards beta, denom is: 	s
					df towards b, denom is:		b[k-1]*s
				*/
				g[k-1] += (-r.at(t,k) * (s-exp_at.at(a[k],t)*b[k-1]) + b[k-1] * exp_at.at(a[k],t) * (ss-r.at(t,k)))/s;
			}
		}
		
		if(!g.is_finite())
		{
			b.print("b:");
			fflush(stdout);
			Rcpp::stop("inf gradient");
		}
	}

	void hess(const arma::vec& beta, arma::mat& h)
	{
		h.zeros();
		const arma::vec b = arma::exp(beta);
		
		
		for(int t=0; t<nt; t++)
		{
			// sums for theta[t]
			double s=1,ss = r.at(t,0);
			for(int k=1; k<ncat;k++)
			{
				s += b[k-1] * exp_at.at(a[k],t);
				ss += r.at(t,k);
			}
			for(int i=1;i<ncat;i++)
			{
				//diagonal
				double s_min = s-b[i-1]*exp_at.at(a[i],t);
				double num = r.at(t,i) * s_min * s;
				num += b[i-1] * r.at(t,i) * s_min * exp_at.at(a[i],t);
				num -= SQR(b[i-1]) * (ss - r.at(t,i)) * SQR(exp_at.at(a[i],t));
				
				h.at(i-1,i-1) += num/( SQR(b[i-1])*SQR(s));
				//upper tri
				for(int j=i+1;j<ncat;j++)
					h.at(i-1,j-1) -= exp_at.at(a[i])*exp_at.at(a[j])*ss/SQR(s);
			}
		}	
		for(int i=0;i<ncat-1;i++)
			for(int j=i+1;j<ncat-1;j++)
				h.at(j,i) = h.at(i,j);
		
	}

};

#endif
