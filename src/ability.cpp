#include <RcppArmadillo.h>
#include "shared.h"
using namespace arma;


//to do: check categories match

// [[Rcpp::export]]
double E_2pl(const double theta, const arma::mat& aA, const arma::mat& b, const arma::ivec& items, const arma::ivec& ncat)
{
	const int nit = items.n_elem;
	double ws=0;
	for(int ix=0; ix<nit; ix++)
	{
		int i = items[ix];
		double num=0,den=1;
		for(int k=1; k<ncat[i]; k++)
		{
			double p = std::exp(aA.at(k,i)*(theta-b.at(k,i)));
			num += aA.at(k,i)*p;
			den += p;
		}
		ws += num/den;
	}
	return ws;	
}


// [[Rcpp::export]]
double Ew_2pl(const double theta, const arma::mat& aA, const arma::mat& b, const arma::ivec& items, const arma::ivec& ncat)
{
	const int nit = items.n_elem;
	double E=0,I=0,J=0;
	for(int ix=0; ix<nit; ix++)
	{
		int i = items[ix];
		double S=1,Sa=0, Sa2=0, Sa3=0;
		for(int k=1; k<ncat[i]; k++)
		{
			double p = std::exp(aA.at(k,i)*(theta-b.at(k,i)));
			S += p;
			Sa += p*aA.at(k,i);
			Sa2 += p*SQR(aA.at(k,i));
			Sa3 += p*CUB(aA.at(k,i));
		}
		I += (S*Sa2-SQR(Sa))/SQR(S);
		J -= (-SQR(S)*Sa3+3*S*Sa*Sa2-2*CUB(Sa))/CUB(S);
		E += Sa/S;
	}
	return E-(J/(2*I));	
}


template<bool WTH>
double get_theta(const double s, const mat& aA, const mat& b, const ivec& items, const ivec& ncat, int &err)
{
	std::function<double(const double theta, const mat& aA, const mat& b, const ivec& items, const ivec& ncat)> func;
	
	if(WTH) func = Ew_2pl;
	else func = E_2pl;
	
	double xl = .5, rts = -.5;	
	double fl = func(xl, aA, b, items, ncat),
		   f = func(rts, aA, b, items, ncat);
	
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
		f = func(rts, aA, b, items, ncat);
			
		if(std::abs(dx) < acc)
			break;
	}
	return rts;
};

template< bool WTH>
vec templ_theta_2pl(const imat& a, const vec& A, const mat& b, const ivec& ncat,
			const ivec& pni, const ivec& pcni, ivec& pi, const ivec& px)
{
	const int np = pni.n_elem, nit = A.n_elem;
	vec theta(np);
	int errors=0;
	mat aA(a.n_rows,nit,fill::zeros);
	for(int i=0; i<nit; i++)
		for(int k=1;k<ncat[i];k++)
			aA.at(k,i) = a.at(k,i)* A[i];
#pragma omp parallel for reduction(+:errors)
	for(int p=0; p<np;p++)
	{
		const ivec items(pi.memptr()+pcni[p], pni[p], false, true);
		double ws=0;
		int err=0;
		bool ms=true;
		for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
		{
			ws += aA.at(px[indx],pi[indx]);
			if(!WTH)
				ms = ms && (px[indx] == ncat[pi[indx]]-1);
		}
		if(!WTH && ws==0)
			theta[p] = -datum::inf;
		else if(!WTH && ms)
			theta[p] = datum::inf;
		else
			theta[p] = get_theta<WTH>(ws, aA, b, items, ncat,err);
		errors += err;
	}
	if(errors>0)
		Rcpp::stop("WLE estimates do not converge");
	return theta;
};	


// [[Rcpp::export]]
arma::vec theta_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat,
					const arma::ivec& pni, const arma::ivec& pcni, arma::ivec& pi, const arma::ivec& px,
					const bool WLE=false)
{
	if(WLE)
		return templ_theta_2pl<true>(a, A, b, ncat, pni, pcni, pi, px);
	else
		return templ_theta_2pl<false>(a, A, b, ncat, pni, pcni, pi, px);
}
