#ifndef DXM_PL2_ITEM_
#define DXM_PL2_ITEM_

#include <RcppArmadillo.h>
#include "shared.h"

arma::mat pl2_trace(const arma::vec& theta, const arma::ivec& a, const double A,  const arma::vec& b, const int ncat);
void pl2_trace(const arma::vec& theta, const arma::ivec& a, const double A,  const arma::vec& b, const int ncat, arma::mat& out);

arma::cube pl2_trace_GH(const arma::mat& theta, const arma::ivec& a, const double A, const arma::vec& b, const int ncat);

void pl2_icc(const arma::vec& theta, const arma::ivec& a, const double A, const arma::vec& b, const int ncat, 
				arma::mat& itrace, double* nc_ptr, double* nca_ptr, double* ncab_ptr);

struct ll_pl2_base
{
	int A_prior;
	double amu,asig,asig2;
	ll_pl2_base(const int a_prior=0, const double a_mu=0, const double a_sigma=0.5) : A_prior(a_prior), amu(a_mu), asig(a_sigma)
	{
		asig2=SQR(a_sigma);
	}
	
	double prior_part_ll(const double A)
	{
		if(A_prior==1) //lognormal
		{
			if(A<=0) return std::numeric_limits<double>::infinity();
			return std::log(asig) - std::log(std::exp(-SQR(amu-std::log(A))/(2*asig2))/A) + std::log(2.0)/2 + std::log(arma::datum::pi)/2; 
		} 
		else if(A_prior==2) //normal
		{
			return -R::dnorm(A,amu,asig,true);
		}
		return 0;
	}
	
	double prior_part_df(const double A)
	{		
		if(A_prior==1) //lognormal
		{
			if(A<=0) return -std::numeric_limits<double>::infinity();
			return (asig2 + std::log(A) - amu)/(A*asig2);
		} 
		else if(A_prior==2) //normal
		{
			return (A-amu)/SQR(asig);
		}
		return 0;
	}
	
	double prior_part_hess(const double A)
	{
		if(A_prior==1) //lognormal
		{
			if(A<=0) return std::numeric_limits<double>::infinity();
			return (amu-asig2 - std::log(A) +1)/(SQR(A)*asig2);
		} 
		else if(A_prior==2) //normal
		{
			return 1/SQR(asig);
		}
		return 0;
	}
};

struct ll_pl2 : ll_pl2_base
{
	arma::mat r;
	arma::vec p, theta, epb, typical_size;
	arma::ivec a;	
		
	int nt, ncat;
	
	// a includes 0 cat
	ll_pl2(int* a_ptr, double* theta_ptr, arma::mat& r_, const int a_prior=0, const double a_mu=0, const double a_sigma=0.5) :
		ll_pl2_base(a_prior, a_mu, a_sigma)
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
				ll -= r.at(t,k) * (std::log(p[k])-std::log(s));
		}
		ll += prior_part_ll(A);
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
		g[0] += prior_part_df(A);
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
			//h.at(0,0) += sr * (s*s2+s*s3+s*s4-SQR(s5)-2*s5*s6-SQR(s6))/SQR(s);
			//re-arrange to prevent overflow
			h.at(0,0) += sr*((s2+s3+s4)/s-SQR(s5/s)-SQR(s6/s)-2*((s5/s)*(s6/s)));
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
				//h.at(k,k) -= (SQR(p[k]*a[k])/s - SQR(a[k])*p[k])*A2*sr/s;
				//re-arrange to prevent overflow
				h.at(k,k) -= SQR(a[k])*p[k]*(p[k]/s-1)*A2*sr/s;
				for(int j=k+1;j<ncat;j++)
					h.at(k,j) -= A2*a[k]*a[j]*(p[k]*p[j]*sr)/SQR(s);
					
			}		
		}	

		for(int i=0; i<ncat; i++)
			for(int j=i+1; j<ncat; j++)
				h.at(j,i) = h.at(i,j);
		
		h.at(0,0) += prior_part_hess(A);
		
		if(!negative)
			h *= -1;		
	}
};

// far more costly but hopefully more stable for small n or other problems
struct ll_pl2_v2 : ll_pl2_base
{
	arma::ivec x,a;
	arma::vec theta, typical_size;
	arma::mat posterior;
	int np,ncat,nt;	

	ll_pl2_v2(const arma::field<arma::mat>& itrace, const arma::vec& theta_, 
			const arma::ivec& ip, const arma::ivec& pi, const arma::ivec& pcni, const arma::ivec& px, 
			const arma::ivec& pgroup, const arma::ivec& inp, const arma::ivec& icnp,
			const arma::vec& mu, const arma::vec& sigma, const int item, 
			const arma::ivec& a_item, const int a_prior=0, const double a_mu=0, const double a_sigma=0.5) :
		ll_pl2_base(a_prior, a_mu, a_sigma)
	{
		const int ng = mu.n_elem;
		nt = theta_.n_elem;
		a = a_item;
		ncat = a.n_elem;
		np = inp[item];
		theta = theta_;
		posterior = arma::mat(nt,np);
		x = arma::ivec(np);
		arma::mat posterior0(nt,ng);
		typical_size = arma::vec(ncat, arma::fill::ones);
		
		for(int g=0; g<ng; g++)
			posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);
		
		for(int ii=icnp[item],pp=0; ii < icnp[item+1]; ii++,pp++)
		{
			const int p=ip[ii];
			posterior.col(pp) = posterior0.col(pgroup[p]);			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
			{
				if(pi[indx]!=item)
					posterior.col(pp) %= itrace(pi[indx]).col(px[indx]);
				else
					x[pp] = px[indx];
			}
		}		
		
	}
	//negative ll
	double operator()(const arma::vec& pars)
	{
		double ll=0;
		arma::vec b = pars;
		b[0]=0;
		double A=pars[0];
		
		arma::mat itrace = pl2_trace(theta, a, A, b, ncat);
		
		for(int p=0; p<np;p++)
			ll -= std::log(arma::accu(posterior.col(p) % itrace.col(x[p])));

		return ll + prior_part_ll(A);;
	}
	void df(const arma::vec& pars, arma::vec& g)
	{
		g.zeros();
		arma::vec b = pars;
		b[0]=0;
		double A=pars[0];
		
		arma::vec ptr(nt);
		
		arma::mat itr = pl2_trace(theta, a, A, b, ncat);
		arma::vec C = arma::sum(itr,1);
		arma::mat atb(nt,ncat,arma::fill::zeros);
		for(int k=1; k<ncat;k++)
			atb.col(k) = a[k] * (theta-b[k]); 
		
		for(int p=0; p<np;p++)
		{
			ptr = itr.col(x[p]) % posterior.col(p);
			double dnm = arma::accu(ptr);
			// of course this can be largely precomputed
			g[0] -=arma::accu((atb.col(x[p]) - arma::sum(atb % itr,1)) % ptr)/dnm;
			
			for(int k=1;k<ncat;k++)
			{
				g[k] -= arma::accu(A * a[k] * (itr.col(k) - kron(k,x[p])) % ptr)/dnm;
			}
		}
		g[0] += prior_part_df(A);
	}
	void hess(const arma::vec& pars, arma::mat& h, const bool negative=true)
	{
		h.zeros();
		arma::vec b = pars;
		b[0]=0;
		const double A=pars[0];
		
		arma::vec nconst(nt),nconst_a(nt),nconst_ab(nt);
		arma::mat itr(nt,ncat);
		arma::cube itr2(nt,ncat,ncat);
		for(int k=0; k<ncat; k++)
		{
			itr2.slice(k).col(k) = square(itr.col(k));
			pl2_icc(theta, a, A, b, ncat, itr, nconst.memptr(), nconst_a.memptr(), nconst_ab.memptr());
			for(int l=k+1; l<ncat; l++)
			{
				itr2.slice(k).col(l) = itr.col(k) % itr.col(l);
				itr2.slice(l).col(k) = itr2.slice(k).col(l);
			}
		}

		arma::vec ptr(nt);		
		arma::mat atb(nt,ncat,arma::fill::zeros);
		for(int k=1; k<ncat;k++)
			atb.col(k) = a[k] * (theta-b[k]); 
		
		arma::mat item_const1(nt, ncat);
		for(int x1=0; x1<ncat; x1++)
			item_const1.col(x1) = -(nconst_a % theta - nconst_ab)/nconst - a[x1] * ( b[x1] - theta);
		
		arma::mat AA(nt,ncat);
		arma::vec bb(ncat);
		arma::cube Ab(nt,ncat,ncat);
		for(int x1=0; x1<ncat; x1++)
		{
			AA.col(x1) = arma::square(item_const1.col(x1)) - arma::sum(arma::square(atb) % itr,1) + square(arma::sum(atb % itr,1));
			for(int k=1;k<ncat;k++)
				Ab.slice(k).col(x1) = a[k] * ( kron(k, x1)*(a[k]*(b[k] -theta)-1) \
											+ itr.col(k) % (a[x1]*(theta-b[x1]) +  a[k]*(theta-b[k]) + 1) \
											+ (2*itr.col(k)- kron(k, x1)) % (nconst_ab - theta%nconst_a)/nconst);
		}
		
		
		
		for(int p=0;p<np;p++)
		{
			ptr = itr.col(x[p]) % posterior.col(p);
			ptr = ptr/accu(ptr);
			const int x1 = x[p];
			
			h.at(0,0) -= accu( AA.col(x1) % ptr) - SQR(accu(item_const1.col(x1) % ptr));
				
			for(int k=1;k<ncat;k++)
				bb[k] = accu(-(kron(k, x1)*a[x1] - itr.col(k)*a[k]) % ptr);
		
			for(int k=1;k<ncat;k++)
			{
				for(int l=k; l<ncat; l++)
				{
					double bb2 = a[k]*a[l]*accu(ptr % (kron(l, x1)*kron(k, x1) - kron(l, x1)*itr.col(k) - kron(k, l)*itr.col(k) - kron(k, x1)*itr.col(l) + 2*itr2.slice(k).col(l) ));
																
					h.at(k,l) -= SQR(A) * (bb2 - bb[k] * bb[l]);
				}
				//ab
				h.at(0,k) -= A * (accu(Ab.slice(k).col(x1) % ptr) - bb[k] * accu(item_const1.col(x1) % ptr));

			}
		}	
			

		for(int i=0; i<ncat; i++)
			for(int j=i+1;j<ncat;j++)
				h.at(j,i)=h.at(i,j);
		if(!negative)
			h=-h;
	}
};


#endif


