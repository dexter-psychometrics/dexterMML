#include <RcppArmadillo.h>
#include "shared.h"
#include "weird_models.h"
#include "posterior.h"
#include "minimize.h"
#include <xoshiro.h>
#include <dqrng_distribution.h>
#include <omp.h>

using namespace arma;
using Rcpp::Named;

//in place
void g_trace(const vec& theta, const double beta, const double lguess, mat& out)
{
	const int nt = theta.n_elem;
	const double eg = std::exp(lguess);
	const double guess = eg/(1+eg);
	
	for(int t=0; t<nt; t++)
	{
		out.at(t,1) = guess + (1-guess)/(1+exp(beta-theta[t]));
	}
	out.col(0) = 1-out.col(1);
}

mat g_trace(const vec& theta, const double beta, const double lguess)
{
	mat out(theta.n_elem,2);	
	g_trace(theta, beta,lguess,out);
	return out;
}

// [[Rcpp::export]]
double test_ll_1plG(arma::vec theta, arma::mat& r, const arma::vec& par)
{

	ll_1G f(theta, r);
	
	return f(par);
} 

// [[Rcpp::export]]
arma::vec test_gradient_1plG(arma::vec theta, arma::mat& r, const arma::vec& par)
{

	ll_1G f(theta, r);
	vec g(par.n_elem);
	f.df(par,g);
	return g;
}

// [[Rcpp::export]]
arma::mat test_hess_1plG(arma::vec theta, arma::mat& r, const arma::vec& par)
{

	ll_1G f(theta, r);
	vec g(par.n_elem);
	f.df(par,g);
	mat h(par.n_elem,par.n_elem);
	f.hess(par,h);
	return h;
}



// [[Rcpp::export]]
Rcpp::List estimate_1plG(const arma::vec& b_start, 
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						const arma::vec& theta_start, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const int ref_group=0, const int max_iter = 200)
{
	const int nit = b_start.n_elem, nt = theta_start.n_elem, np = pni.n_elem, ng=gn.n_elem;	
	
	vec sigma = sigma_start, mu=mu_start, theta=theta_start;
	
	field<mat> itrace(nit);
	
	vec b = b_start;
	vec lg(nit);
	lg.fill(-2.2);
	vec old_b=b, old_lg = lg;

	
	field<mat> r(nit);
	for(int i=0; i<nit; i++)
		r(i) = mat(nt,2);
	
	vec thetabar(np,fill::zeros);
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const double tol = 1e-10;
	int iter = 0,min_error=0,stop=0;
	long double ll, old_ll=-std::numeric_limits<long double>::max();
	double maxdif_b=0,maxdif_g=0;
	
	bool adapt_theta = true;
	
	
	
	for(; iter<max_iter; iter++)
	{
		
		
		for(int i=0; i<nit; i++)
			itrace(i) = g_trace(theta, b[i], lg[i]);
		
		estep(itrace, pni, pcni, pi, px, theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);

		if(ll < old_ll)
		{
			if(adapt_theta) adapt_theta = false;
			else
			{
				stop += 2;
				break;
			}	
		}
		
		old_b=b;
		old_lg=lg;
#pragma omp parallel for reduction(+:min_error)
		for(int i=0; i<nit; i++)
		{	
			//if(item_fixed[i] == 1)
			//	continue;

			ll_1G f(theta, r(i));
			vec pars(2);
			pars[0] = lg[i];
			pars[1] = b[i];
			int itr=0,err=0;
			double ll_itm=0;
			nlm(pars, tol, itr, ll_itm, f, err);	
			min_error += err;
			b[i] = pars[1];
			lg[i] = pars[0];
		}
		if(min_error>0)
		{
			stop += 1;
			break;
		}
		for(int g=0;g<ng;g++)
		{			
			mu[g] = sum_theta[g]/gn[g];		
			sigma[g] = std::sqrt(sum_sigma2[g]/gn[g] - mu[g] * mu[g]);
		}
			
		if(ref_group >= 0) mu[ref_group] = 0;
		
		if(adapt_theta)	scale_theta(mu, sigma, gn, theta_start, theta);
	
		
		maxdif_b = abs(b-old_b).max();
		maxdif_g = abs(exp(lg)/(1+exp(lg)) - exp(old_lg)/(1+exp(old_lg))).max();
	
		if(maxdif_b < 0.0001 && maxdif_g < 0.0001)
			break;
		
		old_ll = ll;
		
	}
	if(iter>=max_iter-1)
		stop += 4;
	
	
	return Rcpp::List::create(Named("lg")=lg, Named("b")=b, Named("thetabar") = thetabar, Named("mu") = mu, Named("sigma") = sigma, 
								Named("niter")=iter, Named("err")=stop, Named("theta")=theta,
		Named("debug")=Rcpp::List::create( Named("error")=stop, Named("maxdif_g")=maxdif_b,Named("maxdif_b")=maxdif_b,Named("ll")=ll,Named("r")=r)); 
}


cube g_trace_GH(const mat& theta, const double b, const double lg)
{
	const int ng = theta.n_cols;
	cube out(theta.n_rows, 2, ng);
	
	for(int g=0;g<ng;g++)
		out.slice(g) = g_trace(theta.col(g), b, lg);

	return out;
}


// [[Rcpp::export]]
double loglikelihood_1plG_GH(const arma::vec& b, const arma::vec& lguess,
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& GH_theta, const arma::vec& w, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup)
{
	const int nit = b.n_elem, nt = GH_theta.n_elem, ng=mu.n_elem;
	ivec ncat(nit);
	ncat.fill(2);
	
	mat theta(nt,ng);
	for(int g=0;g<ng;g++)
		theta.col(g) = GH_theta * sigma[g] + mu[g];
	
	field<cube> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = g_trace_GH(theta, b[i], lguess[i]);

	return (double)(loglikelihood_GH(itrace, ncat, pni, pcni, pi, px, theta, w, mu, sigma, pgroup));
}

#define SEED std::round(R::runif(0,1) * 2147483647)

// simulation
// [[Rcpp::export]]
arma::imat sim_1plGc(const arma::vec& b, const arma::vec& guess,	const arma::vec& theta)
{
	const int nit=b.n_elem, np=theta.n_elem;
	
	imat out(np,nit);
	dqrng::xoshiro256plus rng(SEED); 	
	dqrng::uniform_distribution prl_runif(0, 1);
#pragma omp parallel
	{	
		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(omp_get_thread_num() + 1);
#pragma omp for
		for(int i=0;i<nit;i++)
		{		
			const double c1=guess[i], c2 = 1-guess[i];
			for(int p=0;p<np;p++)
			{
				out.at(p,i) = (int)(prl_runif(lrng) < c1+c2/(1+exp(b[i]-theta[p])));
			}
		}
	}
	return out;
}

// ****************************************** 1PL AG ******************************************

//in place
void ag_trace(const vec& theta, const double guess,const double alpha, const double beta,  mat& out)
{
	const int nt = theta.n_elem;

	for(int t=0; t<nt; t++)
	{
		const double e = 1/(1+std::exp(beta-theta[t])); 
		out.at(t,1) = e + (1-e)/(1+std::exp(-guess-alpha*theta[t]));
	}
	out.col(0) = 1-out.col(1);
}

mat ag_trace(const vec& theta, const double guess, const double alpha, const double beta)
{
	mat out(theta.n_elem,2);	
	ag_trace(theta, guess, alpha, beta ,out);
	return out;
}

// [[Rcpp::export]]
double test_ll_AG(arma::vec theta, arma::mat& r, const arma::vec& par,const double a)
{

	ll_1AG f(theta, r,a);
	
	return f(par);
} 

// [[Rcpp::export]]
arma::vec test_gradient_AG(arma::vec theta, arma::mat& r, const arma::vec& par,const double a)
{

	ll_1AG f(theta, r,a);
	vec g(par.n_elem);
	f.df(par,g);
	return g;
}

// [[Rcpp::export]]
arma::mat test_hess_AG(arma::vec theta, arma::mat& r, const arma::vec& par,const double a)
{

	ll_1AG f(theta, r,a);
	vec g(par.n_elem);
	f.df(par,g);
	mat h(par.n_elem,par.n_elem);
	f.hess(par,h);
	return h;
}

// [[Rcpp::export]]
double test_ll_AG_alpha(arma::vec& a, arma::vec theta, arma::vec& g, arma::vec& b, arma::field<arma::mat>& r)
{

	ll_1AG_alpha f(theta, g,b,r);
	
	return f(a);
} 

// [[Rcpp::export]]
double test_gradient_AG_alpha(arma::vec& a, arma::vec theta, arma::vec& g, arma::vec& b, arma::field<arma::mat>& r)
{

	ll_1AG_alpha f(theta, g,b,r);
	vec df(1);
	f.df(a,df);
	return df[0];
}

// [[Rcpp::export]]
double test_hess_AG_alpha(arma::vec& a, arma::vec theta, arma::vec& g, arma::vec& b, arma::field<arma::mat>& r)
{

	ll_1AG_alpha f(theta, g,b,r);
	mat h(1,1);
	f.hess(a,h);
	return h[0];
}


// [[Rcpp::export]]
Rcpp::List estimate_1plAG(const arma::vec& b_start, 
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						const arma::vec& theta_start, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const int ref_group=0, const int max_iter = 200)
{
	const int nit = b_start.n_elem, nt = theta_start.n_elem, np = pni.n_elem, ng=gn.n_elem;	
	
	vec sigma = sigma_start, mu=mu_start, theta=theta_start;
	
	field<mat> itrace(nit);
	double a= 0.1;
	double old_a=a;
	vec b = b_start;
	vec g(nit);
	g.fill(2);

	vec old_b=b, old_g = g;

	
	field<mat> r(nit);
	for(int i=0; i<nit; i++)
		r(i) = mat(nt,2);
	
	vec thetabar(np,fill::zeros);
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const double tol = 1e-10;
	int iter = 0,min_error=0,stop=0;
	long double ll, old_ll=-std::numeric_limits<long double>::max();
	double maxdif_b=0,maxdif_g=0,maxdif_a=0;
	
	bool adapt_theta = true;
	
	
	
	for(; iter<max_iter; iter++)
	{
		
		
		for(int i=0; i<nit; i++)
			itrace(i) = ag_trace(theta, g[i],a,b[i]);
		
		estep(itrace, pni, pcni, pi, px, theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);

		if(ll < old_ll)
		{
			if(adapt_theta) adapt_theta = false;
			else
			{
				stop += 2;
				break;
			}	
		}
		
		old_b=b;
		old_g=g;
		old_a=a;
//pragma omp parallel for reduction(+:min_error)
		for(int i=0; i<nit; i++)
		{	
			//if(item_fixed[i] == 1)
			//	continue;
			ll_1AG f(theta, r(i),a);
			vec pars(2);
			pars[0] = g[i];
			pars[1] = b[i];
			int itr=0,err=0;
			double ll_itm=0;
			nlm(pars, tol, itr, ll_itm, f, err);	
			min_error += err;
			b[i] = pars[1];
			g[i] = pars[0];
		}
		double lla=0;
		int ittr=0;
		
		vec a_par(1);
		a_par[0]=a;
		ll_1AG_alpha f_alpha(theta,g,b,r);
		D1min(a_par, tol, ittr, lla, f_alpha, min_error);
		a=a_par[0];
		if(min_error>0)
		{
			stop += 1;
			break;
		}
		
		
		for(int g=0;g<ng;g++)
		{			
			mu[g] = sum_theta[g]/gn[g];		
			sigma[g] = std::sqrt(sum_sigma2[g]/gn[g] - mu[g] * mu[g]);
		}
			
		if(ref_group >= 0) mu[ref_group] = 0;
		
		if(adapt_theta)	scale_theta(mu, sigma, gn, theta_start, theta);
	
		
		maxdif_g = abs(g-old_g).max();
		maxdif_a = std::abs(old_a-a);
		maxdif_b = abs(b-old_b).max();
		
	
		if(maxdif_b < 0.0001 && maxdif_g < 0.0001 && maxdif_a < 0.0001 )
			break;
		
		old_ll = ll;
	}
	if(iter>=max_iter-1)
		stop += 4;
	
	
	return Rcpp::List::create(Named("g")=g, Named("b")=b,Named("alpha")=a, Named("thetabar") = thetabar, Named("mu") = mu, Named("sigma") = sigma, 
								Named("niter")=iter, Named("err")=stop, Named("theta")=theta,
		Named("debug")=Rcpp::List::create( Named("error")=stop, Named("maxdif_g")=maxdif_b,Named("maxdif_b")=maxdif_b,Named("ll")=ll,Named("r")=r)); 
}


// simulation
// [[Rcpp::export]]
arma::imat sim_1plAGc(const arma::vec& guess, const arma::vec& a, const arma::vec& b, const arma::vec& theta)
{
	const int nit=b.n_elem, np=theta.n_elem;
	
	imat out(np,nit);
	dqrng::xoshiro256plus rng(SEED); 	
	dqrng::uniform_distribution prl_runif(0, 1);
#pragma omp parallel
	{	
		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(omp_get_thread_num() + 1);
#pragma omp for
		for(int i=0;i<nit;i++)
		{		
			for(int p=0;p<np;p++)
			{
				const double e = 1/(1+std::exp(b[i]-theta[p])) ;
				out.at(p,i) = (int)(prl_runif(lrng) < e + (1-e)/(1+std::exp(-guess[i]-a[i]*theta[p])));
			}
		}
	}
	return out;
}


cube ag_trace_GH(const mat& theta, const double g, const double a, const double b)
{
	const int ng = theta.n_cols;
	cube out(theta.n_rows, 2, ng);
	
	for(int g=0;g<ng;g++)
		out.slice(g) = ag_trace(theta.col(g), g,a,b);

	return out;
}


// [[Rcpp::export]]
double loglikelihood_1plAG_GH(const arma::vec& g, const double a, const arma::vec& b,
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& GH_theta, const arma::vec& w, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup)
{
	const int nit = b.n_elem, nt = GH_theta.n_elem, ng=mu.n_elem;
	ivec ncat(nit);
	ncat.fill(2);
	
	mat theta(nt,ng);
	for(int g=0;g<ng;g++)
		theta.col(g) = GH_theta * sigma[g] + mu[g];
	
	field<cube> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = ag_trace_GH(theta, g[i], a, b[i]);

	return (double)(loglikelihood_GH(itrace, ncat, pni, pcni, pi, px, theta, w, mu, sigma, pgroup));
}

