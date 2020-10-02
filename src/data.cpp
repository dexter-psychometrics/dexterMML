#include <RcppArmadillo.h>


using namespace arma;
using Rcpp::Named;



// will assume no missing categories for now
// gives two data sets:
// person, item, score
// item, person, score
// [[Rcpp::export]]
Rcpp::List mat_pre(arma::imat& dat)
{
	const int nit = dat.n_cols, np = dat.n_rows;
	
	ivec isum(nit,fill::zeros), imax(nit,fill::zeros), inp(nit,fill::zeros);
	ivec psum(np,fill::zeros), pni(np,fill::zeros);

// margins
	for(int i=0; i<nit; i++)
	{
		const ivec rsp(dat.colptr(i), np, false, true);

		for(int p=0; p<np; p++)
		{
			if(rsp[p]>=0)
			{
				isum[i] += rsp[p];
				psum[p] += rsp[p];				
				inp[i]++;
				pni[p]++;
				imax[i] = std::max(imax[i],rsp[p]);				
			}		
		}	
	}

	// cumulative pointers	
	ivec icnp(nit+1), pcni(np+1);
	icnp[0] = 0;
	pcni[0] = 0;
	std::partial_sum(inp.begin(),inp.end(),icnp.begin()+1);
	std::partial_sum(pni.begin(),pni.end(),pcni.begin()+1);

	// response vectors organized two ways
	const int sz = icnp[nit];
	
	ivec ip(sz), ix(sz), pi(sz),px(sz);
	
	ivec pindx(pcni.memptr(), np, true, true);

	for(int i=0; i<nit; i++)
	{
		const ivec rsp(dat.colptr(i), np, false, true);
		int indx = icnp[i];

		for(int p=0; p<np; p++)
		{
			if(rsp[p]>=0)
			{
				ip[indx] = p;
				ix[indx++] = rsp[p];
				pi[pindx[p]] = i;
				px[pindx[p]++] = rsp[p]; 
			}
		}
	}
	return Rcpp::List::create(
		Named("pi") = pi, Named("px") = px, Named("ip") = ip, Named("ix") = ix,
		Named("inp") = inp, Named("icnp") = icnp, Named("pni") = pni, Named("pcni") = pcni,
		Named("imax") = imax, Named("isum") = isum, Named("psum") = psum);
}



inline double square(const double v){ return v*v;}


// Cohens prox algorithm for dichomous items
// [[Rcpp::export]]
Rcpp::List prox_dich(const arma::ivec& isum, const arma::ivec& psum,
					 const arma::ivec& inp, const arma::ivec& pni,
					 const arma::ivec& icnp, const arma::ivec& pcni,
					 const arma::ivec& ip, const arma::ivec& pi,
					 const int max_iter=20, const double min_change=0.01)
{
	const int np= psum.n_elem, nit = isum.n_elem;

	//logits, persons have extreme score protection, items don't
	vec plogit(np), ilogit(nit);
	
	for(int i=0; i<nit; i++)
		ilogit[i] = std::log( ((double)isum[i]) / (inp[i] - isum[i]));

#pragma omp parallel for	
	for(int p=0;p<np;p++)
	{
		if(psum[p]>0 && psum[p]<pni[p])
			plogit[p] = std::log(((double)(psum[p]))/(pni[p]-psum[p]));			
		else if(psum[p] == 0)
			plogit[p] = std::log(0.5/(pni[p]-0.5));
		else
			plogit[p] = std::log((pni[p]-0.5)/0.5);
	}
	
	vec beta(ilogit.memptr(),nit,true,true), theta(np, fill::zeros);
		
	for(int iter=0; iter<max_iter; iter++)
	{
		if(iter>0)
		{
			//items
#pragma omp parallel for
			for(int i=0; i<nit; i++) 
			{
				double m=0, v=0; 
				for(int p=icnp[i]; p<icnp[i+1]; p++)
					m += theta[ip[p]];
				
				m /= inp[i];
				
				for(int p=icnp[i]; p<icnp[i+1]; p++)
					v += square(theta[ip[p]] - m);
				
				v /= inp[i];
				beta[i] = m - std::sqrt(1 + v/2.9) * ilogit[i];
			}
			beta = beta - mean(beta);
		}	

		//persons
		double max_change = 0;
#pragma omp parallel for reduction(max: max_change)
		for(int p=0; p<np; p++)
		{
			double m=0, v=0;
			for(int i = pcni[p]; i < pcni[p+1]; i++)
				m += beta[pi[i]];
				
			m /= pni[p];
			
			for(int i = pcni[p]; i < pcni[p+1]; i++)
				v += square(beta[pi[i]] - m);
			
			v /= pni[p];
			double old = theta[p];
			theta[p] = m + std::sqrt(1 + v/2.9) * plogit[p];
			max_change = std::max(std::abs(theta[p]-old),max_change);
		}
		if(max_change < min_change)
			break;
	}
	return Rcpp::List::create(
		Named("theta") = theta, Named("beta") = beta);
}


// weights misschien nog doen
// dichotoom logistische regressie voor startwaardes 2PL
// [[Rcpp::export]]
Rcpp::List start_lr(const arma::vec& theta, 
					const arma::ivec& ip, const arma::ivec& ix,
					const arma::ivec& inp, const arma::ivec& icnp,
					const arma::vec& ibeta)
{
	
	const int nit = ibeta.n_elem;

	
	vec beta(ibeta.memptr(),nit);
	vec alpha(nit);
	
	const int max_iter = 20;
	const double convergence = 1e-6;

#pragma omp parallel for
	for(int i=0; i<nit; i++)
	{
		vec g(2);
		mat h(2,2);
		double a = 1;
		double b = -beta[i];
		double ll_old=0;
		for(int iter =0; iter<max_iter; iter++)
		{
			double ll =0;
			g.zeros();
			h.zeros();
			for(int j=icnp[i]; j<icnp[i+1]; j++)
			{
				int p = ip[j];
				double prob = 1/(1+std::exp(-a*theta[p]-b));
				if(ix[j]==1)
					ll += std::log(prob);
				else
					ll += std::log(1-prob);
				
				g[0] += (ix[j]-prob) * theta[p];
				g[1] += ix[j]-prob;
				
				h.at(0,0) += prob * (1-prob) * theta[p] * theta[p];
				h.at(1,0) += prob * (1-prob) * theta[p];
				h.at(1,1) += prob * (1-prob);
			}

			h.at(0,1) = h.at(1,0);
			
			vec delta = solve(h,g);
			
			a += delta.at(0,0);
			b += delta.at(1,0);
			
			if(iter>0 && std::abs(ll-ll_old)<convergence)
				break;
				
			ll_old = ll;
		}

		beta[i]=-b/a;
		alpha[i]=a;
	}
	return Rcpp::List::create(Named("alpha") = alpha, Named("beta") = beta);
}
//https://thelaziestprogrammer.com/sharrington/math-of-machine-learning/solving-logreg-newtons-method


