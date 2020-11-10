#ifndef DXM_POLY2_ITEM_
#define DXM_POLY2_ITEM_

#include <RcppArmadillo.h>
#include "shared.h"

arma::mat poly2_trace(const arma::vec& theta, const arma::ivec& a, const double A,  const arma::vec& b, const int ncat);

// 2pl polytomous with fixed category weights
// only tested with dich yet
struct ll_poly2
{
	arma::mat r;
	arma::vec p, theta, epb, typical_size;
	arma::ivec a;	
	
	int nt;
	int ncat;
	
	// a includes 0 cat
	ll_poly2(int* a_ptr, double* theta_ptr, arma::mat& r_) 
	{
		ncat = r_.n_cols;
		nt = r_.n_rows;
		theta = arma::vec(theta_ptr,nt,false,true);
		p = arma::vec(ncat);
		epb = arma::vec(ncat);
		p[0]=1;
		epb[0]=1;
		a = arma::ivec(a_ptr, ncat,false,true);
		r = arma::mat(r_.memptr(),nt,ncat,false,true);
		typical_size = arma::vec(ncat, arma::fill::ones);
	}
	
	// par = A,b2,b3, etc.
	// minus LL
	double operator()(const arma::vec& par)
	{
		double ll=0;
		const double A = par[0];
		
		for(int t=0; t<nt; t++)
		{
			double s=1;
			p[0] = 1;
			for(int k=1; k<ncat; k++)
			{
				p[k] = std::exp(A*a[k]*(theta[t]-par[k])); 				
				s += p[k];
			}

			for(int k=0; k<ncat; k++)
				ll -= r.at(t,k) * std::log(p[k]/s);
		}

		return ll;
	}
	//gradient of minus LL
	void df(const arma::vec& par, arma::vec& g)
	{
		
		g.zeros();
		const double A = par[0];
		for(int t=0; t<nt; t++)
		{
			double s=1;
			for(int k=1;k<ncat;k++)
			{
				p[k] = std::exp(A*a[k]*(theta[t]-par[k]));
				s += p[k];
			}
			double s1=0;
			for(int k=1; k<ncat; k++)
				for(int j=0; j<ncat; j++) if(k!=j)
				{
					s1 += a[k]*par[k]*r.at(t,k)*p[j];
					s1 -= a[k]*par[k]*r.at(t,j)*p[k];
					s1 -= a[k]*theta[t]*r.at(t,k)*p[j];
					s1 += a[k]*theta[t]*r.at(t,j)*p[k];
				}
			g[0] += s1/s;
			
			for(int k=1; k<ncat; k++)
			{
				double s3=0;
				for(int j=0; j<ncat; j++) if(k!=j)
				{
					s3+=r.at(t,j)*p[k];
					s3-=r.at(t,k)*p[j];
				}
				g[k] -= A*a[k]*s3/s;
			}		
				
				
		
		}
	}
	//negative=true -> hessian of negative ll
	void hess(const arma::vec& par, arma::mat& h, const bool negative=true)
	{
		h.zeros();
		const double A = par[0];
		const double A2 = A*A;
		const arma::ivec a2 = arma::square(a);
		arma::vec b2 = arma::square(par);
		
		b2[0] = 0;
		for(int t=0; t<nt; t++)
		{
			//AA
			const double t2 = SQR(theta[t]);
			double s=1,s2=0,s3=0,s4=0,s5=0,s6=0,sr=r.at(t,0);
			for(int k=1;k<ncat;k++)
			{
				p[k] = std::exp(A*a[k]*(theta[t]-par[k]));
				s += p[k];
				s2 += t2 * a2[k] * p[k];
				s3 += a2[k] * b2[k] * p[k];
				s4 -= 2*theta[t]*a2[k]*par[k]*p[k];
				s5 += theta[t]*a[k]*p[k];
				s6 -= a[k]*par[k]*p[k];
				sr += r.at(t,k);
			}
			h.at(0,0) += sr * (s*s2+s*s3+s*s4-SQR(s5)-2*s5*s6-SQR(s6))/SQR(s);
			//Ab
			for(int k=1;k<ncat;k++)
			{
				double ss=0;
				for(int i=0; i<ncat; i++)
				{
					double ss1=0;
					for(int j=1; j<ncat; j++)
					{
						double ss2=A*theta[t]*p[k]*a[k] - A*par[j]*p[k]*a[k];
						ss1 += p[j]*a[j]*(ss2/s + kron(j,k)*(A*a[j]*par[j]-A*a[j]*theta[t]-1));
					}					
					ss += r.at(t,i) * (kron(i,k)*a[i] + ss1/s);
				}
				h.at(0,k) += ss;	
			}
			//bb
			for(int k=1;k<ncat;k++)
			{
				h.at(k,k) -= (SQR(p[k]*a[k])/s - SQR(a[k])*p[k])*A2*sr/s;
				for(int j=k+1;j<ncat;j++)
					h.at(k,j) -= A2*a[k]*a[j]*(p[k]*p[j]*sr)/SQR(s);
			}		
		}	
		if(!negative)
			h *= -1;
		
		for(int i=0; i<ncat; i++)
			for(int j=i+1; j<ncat; j++)
				h.at(j,i) = h.at(i,j);
	}
};




#endif


