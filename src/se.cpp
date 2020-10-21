
#include <RcppArmadillo.h>
#include "est_2PL_dich.h"
#include "minimize.h"
#include "item_ll.h"
#include "shared.h"


using namespace arma;
using Rcpp::Named;


// helper function to efficiently compute empirical hessian, currently only for 2pl dichotomous
mat group_LL_2pl_dich(const vec& a, const vec& b, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, const mat& mu, const mat& sigma, const ivec& pgroup, const int ref_group)
{
	const int nptb = mu.n_rows, ng=mu.n_cols;
	mat ll(nptb, ng, fill::zeros);
	const int nit = a.n_elem, nt = theta.n_elem, np = pni.n_elem;
	mat itrace(nt,nit);
	
	cube posterior0(nt,nptb,ng);
	for(int g=0; g<ng; g++) if(g!=ref_group)
		for(int p=0; p<nptb; p++)
			posterior0.slice(g).col(p) = gaussian_pts(mu.at(p,g),sigma.at(p,g),theta);
	
	for(int i=0; i<nit; i++)
		itrace.col(i) = 1/(1+exp(-a[i]*(theta-b[i])));
	
	
#pragma omp parallel
	{
		vec posterior(nt);
# pragma omp for reduction(+:ll)
		for(int p=0; p<np;p++) 
		{
			int g = pgroup[p];
			if(g!=ref_group)
			{
				posterior.ones();

				for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				{
					if(px[indx] == 1)
						posterior %= itrace.col(pi[indx]);
					else
						posterior %= 1-itrace.col(pi[indx]);
				}	
				for(int i=0;i<nptb; i++)
					ll.at(i,g) += std::log(accu(posterior % posterior0.slice(g).col(i)));
			}
		}
	}
	return ll;
}



// r0 and r1 are from the last iteration in estimation routine
mat J_2pl_dich(const arma::vec& a_fixed, const arma::vec& b_fixed, 
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				arma::vec& theta, const arma::vec& mu_fixed, const arma::vec& sigma_fixed, const arma::ivec& gn, const arma::ivec& pgroup,
				const int ref_group=0)
{	
	const int nit = a_fixed.n_elem, nt = theta.n_elem, np = pni.n_elem, ng=gn.n_elem;
	
	mat a(nit,2),b(nit,2);
	mat mu(ng,2), sigma(ng,2);
	
	mat r0(nt,nit, fill::zeros), r1(nt,nit, fill::zeros);	
	
	vec thetabar(np,fill::zeros);	
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const double tol = 1e-8;
	double ll;
	
	const double delta = 1e-05;
	vec signed_delta(2);
	signed_delta[0] = -delta;
	signed_delta[1] = delta;

	
	const int npar = 2 * (nit+ng-1);
	mat jacob(npar, npar);


	for(int j=0; j<npar; j++)
	{
		a.col(0) = a_fixed; a.col(1) = a_fixed; 
		b.col(0) = b_fixed; b.col(1) = b_fixed;
		mu.col(0) = mu_fixed;mu.col(1) = mu_fixed;
		sigma.col(0) = sigma_fixed;	sigma.col(1) = sigma_fixed;		
		
		for(int d=0; d<=1; d++)
		{
			if(j<2*nit && j % 2 == 0)
				a.at(j/2,d) += signed_delta[d]; 
			else if(j < 2*nit)
				b.at(j/2, d) += signed_delta[d]; 
			else 
			{
				int g = j/2;
				if(g==ref_group)
					continue;
				if(j % 2 == 0)
					mu.at(g,d) += signed_delta[d]; 
				else
					sigma.at(g,d) += signed_delta[d]; 			
			}
			
			estep_2pl_dich(a.col(d), b.col(d), pni, pcni, pi, px, 
							theta, r0, r1, thetabar, sum_theta, sum_sigma2, mu.col(d), sigma.col(d), pgroup, ll);

#pragma omp parallel for
			for(int i=0; i<nit; i++)
			{				
				ll_2pl_dich f(r1.colptr(i), r0.colptr(i), theta.memptr(), nt);
				vec pars = {a.at(i,d), b.at(i,d)};
				int itr=0;
				double ll_itm=0;
				
				dfpmin(pars, tol, itr, ll_itm, f);
				
				a.at(i,d) = pars[0];
				b.at(i,d) = pars[1];
			}

			for(int g=0;g<ng;g++) 
			{					
				mu.at(g,d) = sum_theta[g]/gn[g];	
				sigma.at(g,d) = std::sqrt(sum_sigma2[g]/gn[g] - mu.at(g,d) * mu.at(g,d)); 
			}
		}
		for(int i=0; i<nit; i++)
		{
			jacob.at(i*2,j) = (a.at(i,1) - a.at(i,0))/(2*delta);
			jacob.at(i*2+1,j) = (b.at(i,1) - b.at(i,0))/(2*delta);
		}
		for(int g=0; g<ng; g++) if(g!= ref_group)
		{			
			int p = 2 * nit;
			if(g>ref_group) p -= 2;
			jacob.at(p+g*2,j) = (mu.at(g,1) - mu.at(g,0))/(2*delta);
			jacob.at(p+g*2+1,j) = (sigma.at(g,1) - sigma.at(g,0))/(2*delta);
		}		
	}
	
	return jacob;
}

// [[Rcpp::export]]
Rcpp::List Oakes_2pl_dich(const arma::vec& a, const arma::vec& b, arma::mat& r0, arma::mat& r1, 
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& gn, const arma::ivec& pgroup,
				const int ref_group=0)
{
	const int ng = mu.n_elem, nit=a.n_elem, nt=theta.n_elem;
	const int npar = (ng+nit-1)*2;
	// observed hessian
	mat obs(npar, npar, fill::zeros), h(2,2);

	// fill items with analytical hessian based on r0,r1 from last em iteration, ask Timo if that is correct
	for(int i=0; i<nit; i++)
	{				
		ll_2pl_dich f(r1.colptr(i), r0.colptr(i), theta.memptr(), nt);
		vec pars(2);
		pars[0] = a[i];
		pars[1] = b[i];
		f.hess(pars,h);
		obs.at(i*2,i*2) = h.at(0,0);
		obs.at(i*2+1,i*2+1) = h.at(1,1);
		obs.at(i*2,i*2+1) = h.at(1,0);
		obs.at(i*2+1,i*2) = h.at(1,0);
	}
	
	// fill groups with empirical hessian, central difference method
	if(ng>1)
	{
		const double d=1e-05;
		const vec delta_mu = {.0, 2*d, -2*d, d, d, -d, -d, .0, .0};
		const vec delta_sigma = {.0, .0, .0, d, -d, -d, d, 2*d, -2*d};
		
		mat mu_ptb(9,ng), sigma_ptb(9,ng);	
		
		for(int g=0; g<ng;g++) 
		{
			mu_ptb.col(g) = mu[g] + delta_mu;
			sigma_ptb.col(g) = sigma[g] + delta_sigma;	
		}
		
		mat ll = group_LL_2pl_dich(a, b, pni, pcni, pi, px, theta, mu_ptb, sigma_ptb, pgroup, ref_group);
		
		int pos = 2*nit;	
		for(int g=0; g<ng;g++) if(g!= ref_group)
		{			
			obs.at(pos,pos) = (ll.at(1,g) - 2* ll.at(0,g) + ll.at(2,g)) / (4*SQR(d));
			obs.at(pos,pos+1) = (ll.at(3,g) - ll.at(4,g) - ll.at(6,g) + ll.at(5,g)) / (8 * SQR(d)); // 2*4=8
			obs.at(pos+1,pos+1) = (ll.at(7,g) - 2*ll.at(0,g) + ll.at(8,g)) / (4 * SQR(d));		
			obs.at(pos+1,pos) = obs.at(pos,pos+1);
			
			pos += 2;
		}
	}
	mat J = J_2pl_dich(a, b, pni, pcni, pi, px, theta, mu, sigma, gn, pgroup, ref_group);
	// J for latent variables diagonal should be 0
	int pos = 2*nit;	
	for(int g=0; g<ng;g++) if(g!= ref_group)
	{
		J.submat(pos,pos,pos+1,pos+1).zeros();
		pos++;
	}
	
	mat hess = obs+(J+J.t())/2;	
	
	return Rcpp::List::create(Named("observed")=obs, Named("J")=J, Named("H")=hess);
}



/* example from chalmers
central_difference2 <- function(par, f, delta, ...){
        np <- length(par)
        hess <- matrix(0, np, np)
        fx <- f(par, ...)
        for(i in seq_len(np)){
            for(j in i:np){
                if(i == j){
                    p <- par
                    p[i] <- p[i] + 2 * delta; s1 <- f(p, ...)
                    p[i] <- p[i] - 4 * delta; s3 <- f(p, ...)
                    hess[i, i] <- (s1 - 2*fx + s3) / (4 * delta^2)
                } else {
                    p <- par
                    p[i] <- p[i] + delta; p[j] <- p[j] + delta; s1 <- f(p, ...)
                    p[j] <- p[j] - 2*delta; s2 <- f(p, ...)
                    p[i] <- p[i] - 2*delta; s4 <- f(p, ...)
                    p[j] <- p[j] + 2*delta; s3 <- f(p, ...)
                    hess[i,j] <- hess[j,i] <- (s1 - s2 - s3 + s4) / (4 * delta^2)
                }
            }
        }
        (hess + t(hess))/2
    }
*/

/*
s1 mu: 2, sigma: 0
s3 mu: -2, sigma: 0
compute i: 1, j: 1

s1 mu: 1, sigma: 1
s2 mu: 1, sigma: -1
s4 mu: -1, sigma: -1
s3 mu: -1, sigma: 1
compute i: 1, j: 2

s1 mu: 0, sigma: 2
s3 mu: 0, sigma: -2
compute i: 2, j: 2
*/

