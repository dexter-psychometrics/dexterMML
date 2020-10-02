
#include <RcppArmadillo.h>
#include "minimize.h"

#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) \
initializer( omp_priv = omp_orig )

#pragma omp declare reduction( + : arma::vec : omp_out += omp_in ) \
initializer( omp_priv = omp_orig )

using namespace arma;
using Rcpp::Named;



// no groups


vec gaussian_pts(const double mu, const double s, const vec& theta)
{
	vec out = exp(-0.5*square((theta - mu)/s));
	out = out / accu(out);
	return out;
}


void estep_2pl_dich(const vec& a, const vec& b, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, mat& r0, mat& r1, vec& thetabar, double& sumsig2, const double mu=0, const double sigma=1)
{
	const int nit = a.n_elem, nt = theta.n_elem, np = pni.n_elem;
	mat itrace(nt,nit);
	
	const vec posterior0 = gaussian_pts(mu,sigma,theta);
	
	r0.zeros();
	r1.zeros();
	
	for(int i=0; i<nit; i++)
		itrace.col(i) = 1/(1+exp(-a[i]*(theta-b[i])));
		
	vec sigma2(nt, fill::zeros);
	
#pragma omp parallel
	{
		vec posterior(nt);
# pragma omp for reduction(+:r0,r1,sigma2)
		for(int p=0; p<np;p++)
		{
			for(int i=0; i<nt;i++)
				posterior[i] = posterior0[i];
				
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
			{
				if(px[indx] == 1)
					posterior %= itrace.col(pi[indx]);
				else
					posterior %= 1-itrace.col(pi[indx]);
			}	
			
			posterior = posterior / accu(posterior);
			thetabar[p] = accu(posterior % theta);
			sigma2 += posterior; 
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
			{
				if(px[indx] == 1)
					r1.col(pi[indx]) += posterior;
				else
					r0.col(pi[indx]) += posterior;
			}		
		}
	}

	sumsig2 = accu(sigma2 % square(theta));
}


struct ll_2pl_dich
{
	vec r0,r1,theta;
	int n;

	ll_2pl_dich(double* r1p, double* r0p, double* thetap, const int ni)
	{
		n=ni;
		r1 = vec(r1p,n,false,true);		
		r0 = vec(r0p,n,false,true);
		theta = vec(thetap,n,false,true);	
	}
	
	//returns minus log likelihood
	double operator()(const vec& ab)
	{
		double ll=0;
		const double a = ab[0], b = ab[1];
		for(int i=0;i<n;i++)
		{
			double p = 1/(1+std::exp(-a*(theta[i]-b)));
			ll -= r1[i] * std::log(p) + r0[i] * std::log(1-p);
		}
		return ll;	
	}
	
	//returns gradient of minus ll
	void df(const vec& ab, vec& g)
	{	
		g.zeros();
		const double a = ab[0], b = ab[1];
		for(int i=0;i<n;i++)
		{
			double e = std::exp(a*(b-theta[i]));
			g[0] -= (b-theta[i]) * (r0[i] - r1[i]*e)/(e+1);
			g[1] -= a * (r0[i]-r1[i]*e)/(e+1);
		}
	}
};


// [[Rcpp::export]]
void test(arma::vec& r1, arma::vec& r0, arma::vec& theta, arma::vec& p)
{
	const int nt = theta.n_elem;
	vec g(2);
	
	ll_2pl_dich f(r1.memptr(), r0.memptr(), theta.memptr(), nt);
	
	printf("%f\n", f(p));
	fflush(stdout);
	f.df(p,g);
	g.print("g:");
	fflush(stdout);
	
}

// [[Rcpp::export]]
Rcpp::List estimate_2pl_dich(const arma::vec& a_start, const arma::vec& b_start, 
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const double mu=0, const double sigma=1)
{
	const int nit = a_start.n_elem, nt = theta.n_elem, np = pni.n_elem;
	
	vec a(a_start.memptr(),nit), b(b_start.memptr(),nit);
	
	mat r0(nt,nit, fill::zeros), r1(nt,nit, fill::zeros);
	
	vec pars(2);
	vec thetabar(np,fill::zeros);
	
	double sumsig2;
	
	const int max_iter = 1;
	const double tol = 1e-6;
	
	for(int iter=0; iter<max_iter; iter++)
	{
		estep_2pl_dich(a, b, pni, pcni, pi, px, 
						theta, r0, r1, thetabar, sumsig2, mu, sigma);

		for(int i=0; i<nit; i++)
		{		
			ll_2pl_dich f(r1.colptr(i), r0.colptr(i), theta.memptr(), nt);
			
			pars[0] = a[i];
			pars[1] = b[i];
			int itr=0;
			double ll=0;
			
			dfpmin(pars, tol, itr, ll, f);
			a[i] = pars[0];
			b[i] = pars[1];
		}
	}
	return Rcpp::List::create(Named("a")=a, Named("b")=b, Named("thetabar") = thetabar, Named("sumsig2") = sumsig2);

}

