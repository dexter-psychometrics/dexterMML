#include <RcppArmadillo.h>
#include "shared.h"
using namespace arma;


//to do: check categories match

template<bool WLE, bool USE_A>
double E_2pl(const double theta, const vec& A, const imat& a, const mat& b, const ivec& items, const ivec& ncat)
{
	const int nit = items.n_elem;
	if(WLE)
	{
		double E=0,I=0,J=0;
		if(USE_A)
		{
			for(int ix=0; ix<nit; ix++)
			{
				int i = items[ix];
				double S=1,Sa=0, Sa2=0, Sa3=0;
				for(int k=1; k<ncat[i]; k++)
				{
					double aA = A[i]*a.at(k,i);
					double p = std::exp(aA*(theta-b.at(k,i)));
					S += p;
					Sa += p*aA;
					Sa2 += p*SQR(aA);
					Sa3 += p*CUB(aA);
				}
				I += (S*Sa2-SQR(Sa))/SQR(S);
				J -= (-SQR(S)*Sa3+3*S*Sa*Sa2-2*CUB(Sa))/CUB(S);
				E += Sa/S;
			}		
		}
		else
		{
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
		}	
		return E-(J/(2*I));	
	}
	else
	{
		double ws=0;
		for(int ix=0; ix<nit; ix++)
		{
			int i = items[ix];
			double num=0,den=1;
			for(int k=1; k<ncat[i]; k++)
			{
				double p = std::exp(A[i]*a.at(k,i)*(theta-b.at(k,i)));
				if(USE_A) num += A[i]*a.at(k,i)*p;
				else num += a.at(k,i)*p;
				den += p;
			}
			ws += num/den;
		}
		return ws;		
	}
}

template<bool WLE, bool USE_A>
double get_theta(const double s, const vec& A, const imat& a, const mat& b, const ivec& items, const ivec& ncat, int &err)
{
	double xl = .5, rts = -.5;	
	double fl = E_2pl<WLE, USE_A>(xl, A, a, b, items, ncat),
		   f = E_2pl<WLE, USE_A>(rts, A, a, b, items, ncat);
	
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
		f = E_2pl<WLE, USE_A>(rts, A, a, b, items, ncat);
			
		if(std::abs(dx) < acc)
			break;
	}
	return rts;
};

template< bool WLE, bool USE_A>
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
			if(USE_A) ws += aA.at(px[indx],pi[indx]);
			else ws += a.at(px[indx],pi[indx]);
			
			if(!WLE) ms = ms && (px[indx] == ncat[pi[indx]]-1);
		}
		if(!WLE && ws==0)
			theta[p] = -datum::inf;
		else if(!WLE && ms)
			theta[p] = datum::inf;
		else
			theta[p] = get_theta<WLE, USE_A>(ws, A, a, b, items, ncat,err);
		errors += err;
	}
	if(errors>0)
		Rcpp::stop("WLE estimates do not converge");
	return theta;
};	


// [[Rcpp::export]]
arma::vec theta_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat,
					const arma::ivec& pni, const arma::ivec& pcni, arma::ivec& pi, const arma::ivec& px,
					const bool WLE=false, const bool USE_A=true)
{
	if(WLE && USE_A) return templ_theta_2pl<true, true>(a, A, b, ncat, pni, pcni, pi, px);
	if(!WLE && USE_A) return templ_theta_2pl<false, true>(a, A, b, ncat, pni, pcni, pi, px);
	if(WLE && !USE_A) return templ_theta_2pl<true, false>(a, A, b, ncat, pni, pcni, pi, px);
	//if(!WLE && !USE_A)
	return templ_theta_2pl<false, false>(a, A, b, ncat, pni, pcni, pi, px);
}
