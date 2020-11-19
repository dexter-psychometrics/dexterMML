
#include <RcppArmadillo.h>
#include "minimize.h"
#include "nrm_item.h"
#include "shared.h"

using namespace arma;
using Rcpp::Named;


// loglikelihood for groups based on matrix of perturbed mu and sigma values
// in order to compute empirical gradients or hessian
mat ll_group_nrm(const imat& a, const mat& b, const mat& exp_at, const ivec& ncat, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, const mat& mu, const mat& sigma, const ivec& pgroup)
{
	const int nit = ncat.n_elem, nt = theta.n_elem, np = pni.n_elem, nptb = mu.n_rows, ng=mu.n_cols;

	mat ll(nptb, ng, fill::zeros);
	
	cube posterior0(nt,nptb,ng);
	for(int g=0; g<ng; g++) 
		for(int p=0; p<nptb; p++)
			posterior0.slice(g).col(p) = gaussian_pts(mu.at(p,g),sigma.at(p,g),theta);
	
	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = nrm_trace(theta, a.col(i), b.col(i), ncat[i], exp_at);

	
#pragma omp parallel
	{
		vec posterior(nt);
#pragma omp for reduction(+:ll)
		for(int p=0; p<np;p++)
		{
			int g = pgroup[p];
			posterior.ones();
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				posterior %= itrace(pi[indx]).col(px[indx]);

			for(int i=0;i<nptb; i++)
				ll.at(i,g) += std::log(accu(posterior % posterior0.slice(g).col(i)));
			
		}
	}

	return ll;
}


// estep where not to be update items and groups are skipped according to the design matrix
// this is both faster and prevents rounding errors in the jacobian for 'almost 0' chanhes

void se_estep_nrm(imat& a, mat& b, const mat& exp_at, const ivec& ncat, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, field<mat>& r, vec& sumtheta, vec& sumsig2, const vec& mu, const vec& sigma, const ivec& pgroup,
				const ivec& update_i, const ivec& update_g)
{
	const int nit = ncat.n_elem, nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;
	
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
	{
		itrace(i) = nrm_trace(theta, a.col(i), b.col(i), ncat[i], exp_at);
		r(i).zeros();
	}
	
	mat sigma2(nt, ng, fill::zeros);
	sumtheta.zeros();
	
	
#pragma omp parallel
	{
		vec posterior(nt);
#pragma omp for reduction(+:r,sigma2, sumtheta)
		for(int p=0; p<np;p++)
		{
			int g = pgroup[p];
			posterior = posterior0.col(g);
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				posterior %= itrace(pi[indx]).col(px[indx]);

			posterior = posterior / accu(posterior);			

			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				if(update_i[pi[indx]])
					r(pi[indx]).col(px[indx]) += posterior;

			if(update_g[g])
			{
				sumtheta[g] += accu(posterior % theta);
				sigma2.col(g) += posterior;
			}			
		}
	}

	for(int g=0; g<ng;g++) 
		if(update_g[g])
			sumsig2[g] = accu(sigma2.col(g) % square(theta));
}



mat J_nrm(arma::imat& a, const arma::mat& b_fixed, mat& exp_at, const arma::ivec& ncat,
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const arma::vec& mu_fixed, const arma::vec& sigma_fixed, const arma::ivec& gn, const arma::ivec& pgroup, 
						const arma::imat dsg_ii,const arma::imat dsg_gi, const int ref_group=0)
{
	const int nit = a.n_cols, nt = theta.n_elem, ng=gn.n_elem;
	
	imat dsg_ig = dsg_gi.t();
	
	cube b(b_fixed.n_rows, b_fixed.n_cols, 2);
	
	
	mat mu(ng,2), sigma(ng,2);
	vec sum_theta(ng), sum_sigma2(ng);
	
	field<mat> r(nit);
	for(int i=0; i<nit; i++)
		r(i) = mat(nt,ncat[i]);
	

	const double tol = 1e-10;
	
	const double delta = 1e-05;
	vec signed_delta(2);
	signed_delta[0] = -delta;
	signed_delta[1] = delta;
	
	const ivec gdummy(ng,fill::zeros);
	
	
	const int npar = accu(ncat) - nit + 2*ng-1;
	mat jacob(npar, npar, fill::zeros);
	
	int p=0;
	for(int i=0; i<nit; i++)
	{
		for(int k=1; k<ncat[i]; k++)
		{
			for(int d=0; d<=1; d++)
			{
				b.slice(d) = b_fixed;
				mu.col(d) = mu_fixed;
				sigma.col(d) = sigma_fixed;
				
				b.at(k,i,d) += signed_delta[d];
				
				se_estep_nrm(a, b.slice(d), exp_at, ncat, pni, pcni, pi, px, 
					theta, r, sum_theta, sum_sigma2, mu.col(d), sigma.col(d), pgroup, dsg_ii.col(i), dsg_gi.col(i));

#pragma omp parallel for
				for(int j=0; j<nit; j++) if(dsg_ii.at(j,i) == 1)
				{				
					ll_nrm f(a.colptr(j), exp_at, r(j));
					vec pars(b.slice(d).colptr(j)+1,ncat[j]-1);

					int itr=0,err;
					double ll_itm=0;

					if(ncat[j] == 2)
						D1min(pars, tol, itr, ll_itm, f, err); // 1 dimensional minimization
					else
						nlm(pars, tol, itr, ll_itm, f, err);	

					for(int kj=1;kj<ncat[j];kj++)
						b.at(kj,j,d) = pars[kj-1];
				}
				for(int g=0;g<ng;g++) if(dsg_gi.at(g,i) == 1)
				{
					if(g!= ref_group)
						mu.at(g,d) = sum_theta[g]/gn[g];
					sigma.at(g,d) = std::sqrt(sum_sigma2[g]/gn[g] - mu.at(g,d) * mu.at(g,d)); 
				}
			}
			int q=0;
			for(int j=0; j<nit; j++) 
			{
				if(dsg_ii.at(j,i) == 1)
				{
					for(int kj=1; kj<ncat[j]; kj++)
					{
						jacob.at(q+kj-1,p) = (b.at(kj,j,1) - b.at(kj,j,0))/(2*delta); 
					}
				}
				q += ncat[j]-1;
			}
			for(int g=0; g<ng; g++) 
			{	
				if(dsg_gi.at(g,i) == 1)
				{
					if(g!=ref_group)
						jacob.at(q++,p) = (mu.at(g,1) - mu.at(g,0))/(2*delta);
					jacob.at(q++,p) = (sigma.at(g,1) - sigma.at(g,0))/(2*delta);
				}
				else 
				{
					if(g!=ref_group) q++;
					q++;
				}
			}
			p++;
		}
	}
	// groups, almost a verbatim copy of the above, have to find a neater way maybe
	for(int g=0; g<ng; g++)
	{
		for(int gp = 0; gp<=1; gp++)
		{
			if(gp==0 && g==ref_group)
				continue;
			for(int d=0; d<=1; d++)
			{
				mu.col(d) = mu_fixed;
				sigma.col(d) = sigma_fixed;
				b.slice(d) = b_fixed;
				
				if(gp==0)
					mu.at(g,d) += signed_delta[d];
				else
					sigma.at(g,d) += signed_delta[d];
				
				se_estep_nrm(a, b.slice(d), exp_at, ncat, pni, pcni, pi, px, 
					theta, r, sum_theta, sum_sigma2, mu.col(d), sigma.col(d), pgroup, dsg_ig.col(g), gdummy);

#pragma omp parallel for
				for(int j=0; j<nit; j++) if(dsg_gi.at(g,j) == 1)
				{				
					ll_nrm f(a.colptr(j), exp_at, r(j));
					vec pars(b.slice(d).colptr(j)+1,ncat[j]-1);
					int itr=0,err;
					double ll_itm=0;

					if(ncat[j] == 2)
						D1min(pars, tol, itr, ll_itm, f, err); // 1 dimensional minimization
					else
						nlm(pars, tol, itr, ll_itm, f, err);	

					for(int kj=1;kj<ncat[j];kj++)
						b.at(kj,j,d) = pars[kj-1];
				}
			}
			int q=0;
			for(int j=0; j<nit; j++) 
			{
				if(dsg_ig.at(j,g) == 1)
				{
					for(int kj=1; kj<ncat[j]; kj++)
					{
						jacob.at(q+kj-1,p) = (b.at(kj,j,1) - b.at(kj,j,0))/(2*delta); 
					}
				}
				q += ncat[j]-1;
			}
			p++;
		}	
	}
	
	return jacob;
}		



// [[Rcpp::export]]
Rcpp::List Oakes_nrm(arma::imat& a, const arma::mat& b, const arma::ivec& ncat, arma::field<arma::mat>& r, 
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& gn, const arma::ivec& pgroup,
				const arma::imat& dsg_ii, const arma::imat& dsg_gi,	const int ref_group=0)
{
	const int ng = mu.n_elem, nit=a.n_cols, nt=theta.n_elem;
	const int npar = accu(ncat) - nit + 2*ng-1;
	
	const int max_a = a.max();
	
	//lookup table
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	
	//Jacobian
	mat J = J_nrm(a, b, exp_at, ncat, pni, pcni, pi, px, theta, mu, sigma, gn, pgroup, dsg_ii, dsg_gi, ref_group);
	
	// observed hessian
	const int max_cat = ncat.max();
	mat obs(npar, npar, fill::zeros), h(max_cat-1,max_cat-1);
	
	// fill items with analytical hessian based on r from last em iteration, ask Timo if that is correct
	int p=0;

	for(int i=0; i<nit; i++)
	{				
		ll_nrm f(a.colptr(i), exp_at, r(i));
		vec pars(b.colptr(i)+1,ncat[i]-1);
		
		f.hess(pars,h,false);
		
		for(int j=0; j<ncat[i]-1; j++)
			for(int k=0; k<ncat[i]-1; k++)
				obs.at(p+j,p+k) = h.at(j,k);
		p += ncat[i]-1;
	}

	// empirical observed hessian for groups
	const double d=1e-05;
	const vec delta_mu = {.0, 2*d, -2*d, d, d, -d, -d, .0, .0};
	const vec delta_sigma = {.0, .0, .0, d, -d, -d, d, 2*d, -2*d};
		
	mat mu_ptb(9,ng), sigma_ptb(9,ng);	
		
	for(int g=0; g<ng;g++) 
	{
		mu_ptb.col(g) = mu[g] + delta_mu;
		sigma_ptb.col(g) = sigma[g] + delta_sigma;	
	}

	// dit gaat in 1 keer voor alle groepen
	mat ll = ll_group_nrm(a, b, exp_at, ncat, pni, pcni, pi, px, theta, mu_ptb, sigma_ptb, pgroup);

	for(int g=0; g<ng;g++)
	{
		if(g==ref_group)
		{
			obs.at(p,p) = (ll.at(7,g) - 2*ll.at(0,g) + ll.at(8,g)) / (4 * SQR(d));	
			J.at(p,p)=0;
			p++;
		}
		else
		{
			obs.at(p,p) = (ll.at(1,g) - 2* ll.at(0,g) + ll.at(2,g)) / (4*SQR(d));
			obs.at(p,p+1) = (ll.at(3,g) - ll.at(4,g) - ll.at(6,g) + ll.at(5,g)) / (8 * SQR(d)); // 2*4=8
			obs.at(p+1,p+1) = (ll.at(7,g) - 2*ll.at(0,g) + ll.at(8,g)) / (4 * SQR(d));		
			obs.at(p+1,p) = obs.at(p,p+1);		
			J.submat(p,p,p+1,p+1).zeros();
			p+=2;
		}	
	}


	mat hess = obs+(J+J.t())/2;	

	return Rcpp::List::create(Named("observed")=obs, Named("J")=J, Named("H")=hess);
}


