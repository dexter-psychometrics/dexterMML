#include <RcppArmadillo.h>
#include "shared.h"
using namespace arma;


double E_2plu(const double theta, const arma::vec& A, const arma::imat& a, const arma::mat& b, const arma::ivec& items, const arma::ivec& ncat)
{
	const int nit = items.n_elem;
	double ws=0;
	for(int ix=0; ix<nit; ix++)
	{
		int i = items[ix];
		double num=0,den=1;
		for(int k=1; k<ncat[i]; k++)
		{
			double p = std::exp(A[i]*a.at(k,i)*(theta-b.at(k,i)));
			num += a.at(k,i)*p;
			den += p;
		}
		ws += num/den;
	}
	return ws;	
}


double Ew_2plu(const double theta, const vec& A, const imat& a, const arma::mat& b, const arma::ivec& items, const arma::ivec& ncat)
{
	const int nit = items.n_elem;
	double E=0,I=0,J=0;
	for(int ix=0; ix<nit; ix++)
	{
		int i = items[ix];
		double SpA=1,SpA_a=0,SpA_a2=0,SpA_a3=0;
		for(int k=1; k<ncat[i]; k++)
		{
			double pA = std::exp(A[i]*a.at(k,i)*(theta-b.at(k,i)));
			SpA += pA;
			SpA_a += a.at(k,i)*pA;
			SpA_a2 += SQR(a.at(k,i))*pA;
			SpA_a3 += CUB(a.at(k,i))*pA;

		}
		E += SpA_a/SpA;
		I += (SpA*A[i]*SpA_a2 - SQR(SpA_a)*A[i])/SQR(SpA);
		J -= SQR(A[i])*(-SQR(SpA)*SpA_a3 + 3*SpA*SpA_a*SpA_a2 - 2*CUB(SpA_a))/CUB(SpA);
	} 
	return E-(J/(2*I));	
}


template<bool WTH>
double get_thetau(const double s, const vec& A, const imat& a, const mat& b, const ivec& items, const ivec& ncat, int &err)
{
	std::function<double(const double theta, const vec& A, const imat& a, const arma::mat& b, const arma::ivec& items, const arma::ivec& ncat)> func;
	
	if(WTH) func = Ew_2plu;
	else func = E_2plu;
	
	double xl = .5, rts = -.5;	
	double fl = func(xl, A,a, b, items, ncat),
		   f = func(rts, A,a, b, items, ncat);
	
	double dx;
	
	const int max_iter = 200;
	const double acc = 1e-8; // this might be a little too small but it seems to work
	
	
	for(int iter=0; iter<max_iter; iter++)
	{		
		if((xl > rts) != (fl > f)) //monotonicity
		{
			err=1;
			return rts;
		}
		dx = (xl-rts) * (f-s)/(f-fl);
		xl = rts;
		fl = f;
		rts += std::copysign(std::min(std::abs(dx),.99), dx); 
		f = func(rts, A, a, b, items, ncat);
			
		if(std::abs(dx) < acc)
			break;
	}
	return rts;
};


template< bool WTH>
vec templ_theta_2plu(const imat& a, const vec& A, const mat& b, const ivec& ncat,
			const ivec& pni, const ivec& pcni, ivec& pi, const ivec& px)
{
	const int np = pni.n_elem;
	vec theta(np);
	int errors=0;

#pragma omp parallel for reduction(+:errors)
	for(int p=0; p<np;p++)
	{
		const ivec items(pi.memptr()+pcni[p], pni[p], false, true);
		int ws=0;
		int err=0;
		bool ms=true;
		for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
		{
			ws += a.at(px[indx],pi[indx]);
			if(!WTH)
				ms = ms && (px[indx] == ncat[pi[indx]]-1);
		}
		if(!WTH && ws==0)
			theta[p] = -datum::inf;
		else if(!WTH && ms)
			theta[p] = datum::inf;
		else
			theta[p] = get_thetau<WTH>((double)ws, A, a, b, items, ncat,err);
		errors += err;
	}
	if(errors>0)
		Rcpp::stop("WLE estimates do not converge");
	return theta;
};	


// [[Rcpp::export]]
arma::vec theta_2plu(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat,
					const arma::ivec& pni, const arma::ivec& pcni, arma::ivec& pi, const arma::ivec& px,
					const bool WLE=false)
{
	if(WLE)
		return templ_theta_2plu<true>(a, A, b, ncat, pni, pcni, pi, px);
	else
		return templ_theta_2plu<false>(a, A, b, ncat, pni, pcni, pi, px);
}