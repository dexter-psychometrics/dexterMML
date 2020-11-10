#ifndef DXM_NRM_ITEM_
#define DXM_NRM_ITEM_

#include <RcppArmadillo.h>
#include "shared.h"

arma::mat nrm_trace(const arma::vec& theta, const arma::ivec& a, const arma::vec& b, const int ncat, const arma::mat& exp_at);

// polytomous nrm (can also be used for dichotmous nrm or rasch)

struct ll_nrm
{
	arma::mat r;
	arma::mat eat;
	arma::vec p, typical_size;
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
		typical_size = arma::vec(ncat-1, arma::fill::ones);
	}
	
	// b should not include 0 score
	// minus LL
	double operator()(const arma::vec& b)
	{
		double ll=0;
		const arma::vec eb = arma::exp(b);
		
		for(int t=0; t<nt; t++)
		{
			double s=1;
			p[0] = 1;
			for(int k=1; k<ncat; k++)
			{
				p[k] = eb[k-1] * exp_at.at(a[k],t);
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
	//gradient of minus LL
	void df(const arma::vec& b, arma::vec& g)
	{
		
		g.zeros();
		const arma::vec eb = arma::exp(b);
		
		for(int t=0; t<nt; t++)
		{
			double s=1,ss = r.at(t,0);
			for(int k=1; k<ncat;k++)
			{
				s += eb[k-1] * exp_at.at(a[k],t);
				ss += r.at(t,k);
			}
			for(int k=1; k<ncat;k++)
			{
				/* 
					df towards b, denom is: 	s
					df towards eb, denom is:	eb[k-1]*s
				*/
				g[k-1] += (-r.at(t,k) * (s-exp_at.at(a[k],t)*eb[k-1]) + eb[k-1] * exp_at.at(a[k],t) * (ss-r.at(t,k)))/s;
			}
		}
		
		if(!g.is_finite())
		{
			b.print("b:");
			fflush(stdout);
			Rcpp::stop("inf gradient");
		}
	}
	//negative=true -> hessian of negative ll
	void hess(const arma::vec& b, arma::mat& h, const bool negative=true)
	{
		h.zeros();
		const arma::vec eb = arma::exp(b);
		
		arma::vec d(ncat);
		
		for(int t=0; t<nt; t++)
		{
			// sums for theta[t]
			double s=1,ss = r.at(t,0);
			for(int k=1; k<ncat;k++)
			{
				s += eb[k-1] * exp_at.at(a[k],t);
				ss += r.at(t,k);		
			}
			for(int i=1;i<ncat;i++)
			{
				//diagonal
				double s_min = s-eb[i-1]*exp_at.at(a[i],t);				
				double num = ss*s_min*eb[i-1]*exp_at.at(a[i],t);
				h.at(i-1,i-1) -= num/SQR(s);

				//upper tri
				for(int j=i+1;j<ncat;j++)
					h.at(i-1,j-1) += (ss*eb[i-1]*eb[j-1]*exp_at.at(a[i],t)*exp_at.at(a[j],t))/SQR(s);
			}
		}
		if(!negative)
			h *= -1;
	
		for(int i=0;i<ncat-1;i++)
			for(int j=i+1;j<ncat-1;j++)
				h.at(j,i) = h.at(i,j);
		
	}

};

#endif


