#ifndef DXM_WEIRD_ITEMS_
#define DXM_WEIRD_ITEMS_

#include <RcppArmadillo.h>
#include "shared.h"


// 1plG, dichotomous, scored 0,1
// g + (1-g)/exp(beta-theta)

// define lg: g = exp(lg)/(1+exp(lg)

struct ll_1G
{
	arma::mat r;
	arma::vec theta;
	arma::vec et;
	arma::vec et2;
	arma::vec p, typical_size;
	
	int nt;
	//int ncat; dichotomous
	
	// a includes 0 cat
	ll_1G(arma::vec& theta_,arma::mat& r_) 
	{
		nt = theta_.n_elem;
		p = arma::vec(2);
		r = arma::mat(r_.memptr(),nt,2,false,true);
		theta=arma::vec(theta_.memptr(),nt,false,true);
		et=arma::exp(theta);
		et2 = arma::square(et);
		typical_size = arma::vec(2, arma::fill::ones);
	}
	
	// b should not include 0 score
	// minus LL
	double operator()(const arma::vec& par)
	{
		double ll=0;
		const double guess = std::exp(par[0])/(1+std::exp(par[0]));	

		for(int t=0; t<nt; t++)
		{
			const double p1 = guess + (1-guess)/(1+std::exp(par[1]-theta[t]));
			ll -= (r.at(t,0) * std::log(1-p1) + r.at(t,1) * std::log(p1)); 
		}		

		return ll;
	}
	//gradient of minus LL
	void df(const arma::vec& par, arma::vec& g)
	{
		g.zeros();

		const double guess = std::exp(par[0])/(1+std::exp(par[0])), eb = std::exp(par[1]);	
		const double dfguess = guess - SQR(guess),eb2=SQR(eb);
		const double gdenom1 = (SQR(guess)-guess)*eb; 
		for(int t=0; t<nt; t++)
		{
			g[0] -= (eb*(guess*(r.at(t,0)+ r.at(t,1))-r.at(t,1)) + r.at(t,0)*et[t])/(gdenom1 + (guess-1)*et[t]);
			const double ebt = eb/et[t];
			g[1] -= (guess*ebt*(r.at(t,0)+r.at(t,1))+r.at(t,0)-r.at(t,1)*ebt)/((ebt+1)*(guess*ebt+1));
		}
		g[0] *= dfguess;
	}

	void hess(const arma::vec& par, arma::mat& h)
	{
		h.zeros();

		const double guess = std::exp(par[0])/(1+std::exp(par[0])), eb = std::exp(par[1]);	
		const double dfguess = guess - SQR(guess),eb2=SQR(eb);
		const double elg = std::exp(par[0]);
		
		for(int t=0; t<nt; t++)
		{
			const double bt = eb/et[t];
			const double part1 = elg/(SQR(elg+1)*SQR(elg*bt+elg+1));
			h.at(0,0) += part1 * (r.at(t,0)*(SQR(elg*bt) + 2*SQR(elg)*bt + 2*elg*bt + SQR(elg) + 2*elg + 1) + r.at(t,1)*(SQR(elg*bt) + SQR(elg)*bt - bt));			
			h.at(0,1) -= r.at(t,1)*eb*et[t]/(SQR(guess*eb) +2*guess*eb*et[t] + et2[t]);
			const double part2 = bt/(SQR(bt+1)*SQR(guess*bt+1));
			h.at(1,1) += part2 * (r.at(t,0)*(SQR(guess*bt) + 2*guess*bt+1) + r.at(t,1)*(SQR(guess*bt)-guess*SQR(bt)-guess+1));
		}
		h.at(0,1) = h.at(0,1) * dfguess; 
		h.at(1,0) = h.at(0,1);
	}

};


struct ll_1AG
{
	arma::mat r;
	arma::vec theta;
	arma::vec typical_size;
	double a;
	
	int nt;
	//int ncat; dichotomous
	
	// a includes 0 cat
	ll_1AG(arma::vec& theta_,arma::mat& r_, const double a_fix) 
	{
		nt = theta_.n_elem;
		r = arma::mat(r_.memptr(),nt,2,false,true);
		theta=arma::vec(theta_.memptr(),nt,false,true);
		typical_size = arma::vec(2, arma::fill::ones);
		a=a_fix;
	}
	

	// b should not include 0 score
	// minus LL
	double operator()(const arma::vec& par)
	{
		// par = guess,alpha,beta
		double ll=0;
		const double g = par[0],b=par[1];		

		for(int t=0; t<nt; t++)
		{
			const double e = 1/(1+std::exp(b-theta[t])); 
			const double p1 = e + (1-e)/(1+std::exp(-g-a*theta[t]));
			ll -= r.at(t,0) * std::log(1-p1) + r.at(t,1) * std::log(p1); 
		}		

		return ll;
	}
	//alles hieronder overnieuw het was een -g
	//gradient of minus LL
	void df(const arma::vec& par, arma::vec& grad)
	{
		grad.zeros();
		const double g=par[0],b=par[1];
		const double eg = std::exp(g);
		
		for(int t=0; t<nt; t++)
		{
			const double part1 = std::exp(b+g-a*theta[t]-theta[t]);
			const double eatg = std::exp(g-a*theta[t]), ebt = std::exp(b-theta[t]);
			const double bp =  r.at(t,1) * part1/((eatg+1)*(ebt+eatg+1)) - r.at(t,0)/(eatg+1);
			
			grad[0] += bp;
			//grad[1] -= theta[t] * bp;
			grad[1] += r.at(t,1)*part1/((ebt+1)*(ebt+eatg+1)) - r.at(t,0)/(ebt+1);
		}
	}

	void hess(const arma::vec& par, arma::mat& h)
	{
		h.zeros();
		const double g=par[0],b=par[1];
		const double eb=std::exp(b);
		
		
		for(int t=0; t<nt; t++)
		{
			const double eatg = std::exp(g-a*theta[t]), ebt = std::exp(b-theta[t]);
			const double part1 = eatg*ebt*(SQR(eatg)-ebt-1)/(SQR(eatg+1)*SQR(ebt+eatg+1));

			h.at(0,0) -= r.at(t,1) * part1 - r.at(t,0)*eatg/SQR(eatg+1);
			//h.at(0,1) += r.at(t,1) * theta[t] * part1 - r.at(t,0)*theta[t]*eatg/SQR(eatg+1);
			h.at(0,1) += r.at(t,1) * ebt * eatg/SQR(ebt+eatg+1);
			
			//h.at(1,1) -= r.at(t,1) * SQR(theta[t])*part1 - r.at(t,0) * SQR(theta[t]) * eatg/SQR(eatg+1);
			//h.at(1,2) -= r.at(t,1) * theta[t] *eatg*ebt/SQR(eatg+ebt+1);

			h.at(1,1) += r.at(t,1) * eatg*ebt*(eatg-SQR(ebt)+1)/(SQR(ebt+1)*SQR(ebt+eatg+1)) + r.at(t,0)*ebt/SQR(ebt+1);

			
		}
		
		h.at(1,0) = h.at(0,1);
		//h.at(2,0) = h.at(0,2);
		//h.at(2,1) = h.at(1,2);
	}

};


struct ll_1AG_alpha
{
	arma::field<arma::mat> r;
	arma::vec theta;
	arma::vec typical_size;
	arma::vec g,b;
	
	int nt,nit;
	
	// a includes 0 cat
	ll_1AG_alpha(arma::vec& theta_,arma::vec& g_, arma::vec& b_, arma::field<arma::mat>& r_) 
	{
		nt = theta_.n_elem;
		nit = b_.n_elem;
		g = arma::vec(g_.memptr(),nit,false,true);
		b = arma::vec(b_.memptr(),nit,false,true);
		r = arma::field<arma::mat>(nit);
		for(int i=0;i<nit;i++)
			r(i) = arma::mat(r_(i).memptr(),nt,2);
		
		theta=arma::vec(theta_.memptr(),nt,false,true);
		typical_size = arma::vec(nit, arma::fill::ones);
	}
	

	// b should not include 0 score
	// minus LL
	double operator()(const arma::vec& a_)
	{
		// par = guess,alpha,beta
		double ll=0, a=a_[0];
		for(int i=0;i<nit;i++)
		{
			for(int t=0; t<nt; t++)
			{
				const double e = 1/(1+std::exp(b[i]-theta[t])); 
				const double p1 = e + (1-e)/(1+std::exp(g[i]-a*theta[t]));
				ll -= r(i).at(t,0) * std::log(1-p1) + r(i).at(t,1) * std::log(p1); 
			}		
		}
		return ll;
	}
	//gradient of minus LL
	void df(const arma::vec& a_, arma::vec& grad)
	{
		double gr=0, a=a_[0];

		for(int i=0;i<nit;i++)
			for(int t=0; t<nt; t++)
			{
				const double eatg = std::exp(g[i]-a*theta[t]), ebt = std::exp(b[i]-theta[t]);
				gr -= theta[t] *( r(i).at(t,1) * std::exp(b[i]+g[i]-a*theta[t]-theta[t])/((eatg+1)*(ebt+eatg+1)) - r(i).at(t,0)/(eatg+1));
			}
		grad[0]=gr;
	}

	void hess(const arma::vec& a_, arma::mat& h)
	{
		double a=a_[0],hs=0;
		
		
		for(int i=0;i<nit;i++) for(int t=0; t<nt; t++)
		{
			const double eatg = std::exp(g[i]-a*theta[t]), ebt = std::exp(b[i]-theta[t]);
			const double part1 = eatg*ebt*(SQR(eatg)-ebt-1)/(SQR(eatg+1)*SQR(ebt+eatg+1));

			hs -= r(i).at(t,1) * SQR(theta[t])*part1 - r(i).at(t,0) * SQR(theta[t]) * eatg/SQR(eatg+1);			
		}
		h.at(0,0) = hs;
	}

};



#endif


