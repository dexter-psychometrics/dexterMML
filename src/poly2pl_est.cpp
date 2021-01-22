#include <RcppArmadillo.h>
#include "minimize.h"
#include "poly2pl_item.h"
#include "shared.h"

using namespace arma;
using Rcpp::Named;

void estep_poly2(const imat& a, const vec& A, const mat& b, const ivec& ncat, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
				const vec& theta, field<mat>& r, vec& thetabar, vec& sumtheta, vec& sumsig2, const vec& mu, const vec& sigma, const ivec& pgroup, long double& ll)
{
	const int nit = ncat.n_elem, nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;
	
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
	{
		itrace(i) = poly2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);
		r(i).zeros();
	}
	
	mat sigma2(nt, ng, fill::zeros);
	sumtheta.zeros();
	
	ll=0;

#pragma omp parallel
	{
		vec posterior(nt);
#pragma omp for reduction(+: r, sigma2, sumtheta, ll)
		for(int p=0; p<np;p++)
		{
			int g = pgroup[p];
			posterior = posterior0.col(g);
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				posterior %= itrace(pi[indx]).col(px[indx]);

			double sp = accu(posterior);
			// LL according to Bock/Aitkin 1981 eq (5) and (6)
			ll += std::log(sp); 
			posterior = posterior / sp;
			sumtheta[g] += thetabar[p] = accu(posterior % theta);
			
			sigma2.col(g) += posterior;
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				r(pi[indx]).col(px[indx]) += posterior;

		}
	}

	for(int g=0; g<ng;g++)
		sumsig2[g] = accu(sigma2.col(g) % square(theta));
}

// [[Rcpp::export]]
double loglikelihood_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, 
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup)
{
	const int nit = ncat.n_elem, nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;	
	
	mat posterior0(nt,ng);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);

	field<mat> itrace(nit);
	
	for(int i=0; i<nit; i++)
		itrace(i) = poly2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);

	long double ll=0;

#pragma omp parallel
	{
		vec posterior(nt);
#pragma omp for reduction(+: ll)
		for(int p=0; p<np;p++)
		{
			int g = pgroup[p];
			posterior = posterior0.col(g);
			
			for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
				posterior %= itrace(pi[indx]).col(px[indx]);

			ll += std::log(accu(posterior)); 
		}
	}
	return (double)ll;
}




//more direct objective function


struct ll2_poly
{
	ivec x,ai;
	mat posterior;
	int np;
	vec theta;
	vec bi;	

	ll2_poly(const imat& a, const vec& A, const mat& b, const vec& theta_, const ivec& ncat, 
			const ivec& ip, const ivec& pi, const ivec& pcni, const ivec& px, 
			const ivec& pgroup, const ivec& inp, const ivec& icnp,
			const vec& mu, const vec& sigma, const int item)
	{
		const int nit = b.n_cols, nt=theta_.n_elem, ng = mu.n_elem;;
		np = inp[item];
		field<mat> itrace(nit);
		theta = theta_;
		for(int i=0; i<nit; i++)
		{
			itrace(i) = poly2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);
		}
		posterior = mat(nt,np);
		x = ivec(np);
		mat posterior0(nt,ng);
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
		ai = a.col(item).head(ncat[item]);
		bi = b.col(item);
	}
	double operator()(const arma::vec& pars)
	{
		double ll=0;
		vec bpars = pars;
		bpars[0]=0;
		double A=pars[0];
		
		mat itrace = poly2_trace(theta, ai, A, bpars, ai.n_elem);
		for(int p=0; p<np;p++)
			ll -= std::log(accu(posterior.col(p) % itrace.col(x[p])));
		return ll;
	}
	// can be used to estimate b's the traditional way
	mat r(const arma::vec& pars)
	{
		vec bpars = pars;
		bpars[0]=0;
		double A=pars[0];
		vec pst(theta.n_elem);
		mat itrace = poly2_trace(theta, ai, A, bpars, ai.n_elem);
		mat ipst(theta.n_elem,pars.n_elem,fill::zeros);
		for(int p=0; p<np;p++)
		{
			pst = posterior.col(p) % itrace.col(x[p]);
			ipst.col(x[p]) += pst/accu(pst);
		}	
		return ipst;
	}
};

struct ll2_astep : ll2_poly
{
	using ll2_poly::ll2_poly;
	double operator()(const double A)
	{
		vec par = bi;
		par[0]=A;
		return ll2_poly::operator()(par);
	}
};

//a_prior vergeten we nog even
vec ab_step(imat& a, const vec& A, const mat& b, vec& theta, const ivec& ncat, const ivec& ip, const ivec& pi, const ivec& pcni, const ivec& px, 
			const ivec& pgroup,
			const ivec& inp, const ivec& icnp, const vec& mu, const vec& sigma, const int item)
{
	
	const double tol = 1e-10;
	ll2_astep f(a, A, b, theta, ncat, ip, pi, pcni, px, pgroup, inp, icnp, mu, sigma,item);
	Brent brent;
	
	double alpha = A[item];
	
	brent.bracket(alpha-0.5, alpha+0.5, f);
	
	brent.minimize(f);

	alpha = brent.xmin;
	
	vec pars = b.col(item).head(ncat[item]);
	pars[0] = alpha;
	
	mat r = f.r(pars);
	ll_poly2_b fb(alpha, a.colptr(item), theta.memptr(), r);
	
	vec bpars = pars.tail(pars.n_elem-1);
	int itr=0,err=0;
	double ll_itm=0;			
			
	if(ncat[item]==2)
	{
		D1min(bpars, tol, itr, ll_itm, fb, err);
	}
	else
	{
		nlm(bpars, tol, itr, ll_itm, fb, err);
	}
	pars.tail(pars.n_elem-1)=bpars;
	return pars;
}





// stop estimation because of decreasing likelihood
// quite rare
void est_stop(imat& a, const ivec& ncat, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
		 const ivec& pgroup, vec& theta, const arma::ivec& item_fixed, const vec& h_ll, const int iter, const int ref_group, 
		 const int A_prior, const double A_mu, const double A_sigma, const long double old_ll, const int store_i,
		 /* in & out */
		 vec& A, mat& store_A, mat& b, cube& store_b, vec& mu, mat& store_mu, vec& sigma, mat& store_sigma,
		 field<mat>& r, vec& thetabar, vec& sum_theta, vec& sum_sigma2, long double ll)
{
	const int nit = a.n_cols;
	if(iter>3)
	{
		//quadratic approximation based on likelihood
		mat qx = {{1,0,0},{1,1,1},{1,2,4}};
		vec qy = {h_ll[iter-2], h_ll[iter-1], h_ll[iter]};
		vec qs = solve(qx, qy);
		double qpnt = -qs[1]/(2*qs[2]);
		//store_i, 1-store_i, current
		if(qpnt<=1)
		{
			A = store_A.col(1-store_i) * qpnt + (1-qpnt) * store_A.col(store_i);
			b = store_b.slice(1-store_i) * qpnt + (1-qpnt) * store_b.slice(store_i);
			mu = store_mu.col(1-store_i) * qpnt + (1-qpnt) * store_mu.col(store_i);
			sigma = store_sigma.col(1-store_i) * qpnt + (1-qpnt) * store_sigma.col(store_i);
		}
		else
		{
			qpnt-=1;
			A =  A * qpnt + (1-qpnt) * store_A.col(1-store_i);
			b = b * qpnt + (1-qpnt) * store_b.slice(1-store_i);		
			mu = mu * qpnt + (1-qpnt) * store_mu.col(1-store_i);
			sigma = sigma * qpnt + (1-qpnt) * store_sigma.col(1-store_i);	
		}
		if(ref_group>=0) //prevent rounding errors
		{
			mu[ref_group]=0;
			sigma[ref_group]=1;
		}
		//check if anything is won
		//compute likelihood
		long double ll_new;
		estep_poly2(a, A, b, ncat, pni, pcni, pi, px, 
					theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll_new);
		double prior_part=0;
				
		if(A_prior>0)
			for(int i=0;i<nit;i++)
				if(item_fixed[i] != 1)
				{
					ll_poly2 f(a.colptr(i), theta.memptr(), r(i), A_prior, A_mu, A_sigma);
					prior_part -= f.prior_part_ll(A[i]);
				}
				
		if(ll_new + prior_part < old_ll)
		{
			//reject quadratic approx
			A=store_A.col(1-store_i);
			b=store_b.slice(1-store_i);
			mu=store_mu.col(1-store_i);
			sigma=store_sigma.col(1-store_i);
			// to do: get the old ll without prior part?
		} 
		else
		{
			ll=ll_new;
		}
	}
}




// [[Rcpp::export]]
Rcpp::List estimate_poly2(arma::imat& a, const arma::vec& A_start, const arma::mat& b_start, const arma::ivec& ncat,
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const arma::ivec& item_fixed,
						const arma::ivec ip, const arma::ivec& inp, const arma::ivec& icnp,
						const int ref_group=0, const int A_prior=0, const double A_mu=0, const double A_sigma=0.5, 
						const int max_iter=200, const int pgw=80)
{
	const int nit = a.n_cols, nt = theta.n_elem, np = pni.n_elem, ng=gn.n_elem;
	
	progress_est prog(max_iter, pgw);
	
	mat b = b_start;
	vec A = A_start;

	field<mat> r(nit);
	for(int i=0; i<nit; i++)
		r(i) = mat(nt,ncat[i]);
	
	vec thetabar(np,fill::zeros);
	
	vec sigma = sigma_start, mu=mu_start;
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	const double tol = 1e-10;
	int iter = 0, min_error=0, stop=0;
	long double ll, old_ll=-std::numeric_limits<long double>::max(), prior_part=0;
	double maxdif_A, maxdif_b;
	vec h_ll(max_iter, fill::zeros);
	
	cube store_b(b.n_rows, nit, 2, fill::zeros);
	mat store_A(nit,2,fill::zeros), store_mu(ng,2,fill::zeros), store_sigma(ng,2,fill::zeros);
	int store_i=0;
	
	for(; iter<max_iter; iter++)
	{
		estep_poly2(a, A, b, ncat, pni, pcni, pi, px, 
					theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);
		
		h_ll[iter] = ll + prior_part;
		
		if(ll + prior_part < old_ll)
		{
			est_stop(a, ncat, pni, pcni, pi, px, pgroup, theta, item_fixed, h_ll, iter, ref_group, 
					 A_prior, A_mu, A_sigma, old_ll, store_i,
					/* in & out */
					A, store_A, b, store_b, mu, store_mu, sigma, store_sigma,
					r,  thetabar,  sum_theta,  sum_sigma2, ll);
			stop += 2;
			break;
		}
		old_ll = prior_part+ll;
		maxdif_A=0; maxdif_b=0; prior_part=0;
		store_A.col(store_i)=A;
		store_b.slice(store_i)=b;
		store_mu.col(store_i)=mu;
		store_sigma.col(store_i)=sigma;
		store_i=1-store_i;
		
#pragma omp parallel for reduction(max: maxdif_A, maxdif_b) reduction(+:min_error, prior_part)
		for(int i=0; i<nit; i++)
		{	
			if(item_fixed[i] == 1)
				continue;
			if(inp[i] < 150 && A[i]>0)
			{
				vec pars = ab_step(a, store_A.col(1-store_i), store_b.slice(1-store_i), theta, ncat, ip, pi, pcni, px, pgroup, inp,icnp, mu, sigma, i);
				maxdif_A = std::max(maxdif_A, std::abs(A[i] - pars[0]));
				A[i] = pars[0];
				for(int k=1;k<ncat[i];k++)
				{
					maxdif_b = std::max(maxdif_b, std::abs(b.at(k,i) - pars[k]));
					b.at(k,i) = pars[k];
				}	
				continue;
			} 
			ll_poly2 f(a.colptr(i), theta.memptr(), r(i), A_prior, A_mu, A_sigma);
			vec pars = b.col(i).head(ncat[i]);
			pars[0] = A[i];
			int itr=0,err=0;
			double ll_itm=0;
			
			nlm(pars, tol, itr, ll_itm, f, err);	
			if(A_prior!=1 && (std::abs(A[i]) < .05 || max(abs(pars)) > 50))
			{
				// 2pl can be poorly identified with local minima
				// on opposite sides of A=0, attempt to break out with a restart of nlm
				int err2=0;
				itr=0;
				double ll_itm2=0;
				vec pars2 = -b.col(i).head(ncat[i])/5;
				pars2[0] = -2*A[i];
				nlm(pars2, tol, itr, ll_itm2, f, err2);	
				if(err2==0 && ll_itm2<ll_itm)
				{
					pars=pars2;
					err=err2;
				}			
			}

			min_error += err;
			if(err==0)
			{
				maxdif_A = std::max(maxdif_A, std::abs(A[i] - pars[0]));
				A[i] = pars[0];
				for(int k=1;k<ncat[i];k++)
				{
					maxdif_b = std::max(maxdif_b, std::abs(b.at(k,i) - pars[k]));
					b.at(k,i) = pars[k];
				}		
				prior_part -= f.prior_part_ll(A[i]);
			}
		}
					
		if(min_error>0)
		{
			stop += 1;
			break;
		}
		if(min(abs(A)) < 1e-4)
		{
			stop += 8;
			break;
		}

		for(int g=0;g<ng;g++)
		{		
			if(g==ref_group)
			{
				mu[g] = 0;
				sigma[g] = 1;
			}
			else
			{
				mu[g] = sum_theta[g]/gn[g];		
				sigma[g] = std::sqrt(sum_sigma2[g]/gn[g] - mu[g] * mu[g]);
			}
		}
			
		prog.update(std::max(maxdif_b, maxdif_A), iter);
				
		if(maxdif_b < .0001 && maxdif_A < .0001)
			break;

	}
	if(iter>=max_iter-1)
		stop += 4;
	
	prog.close();
	
	return Rcpp::List::create(Named("A")=A, Named("b")=b, Named("thetabar") = thetabar, Named("mu") = mu, Named("sd") = sigma, 
									Named("r")=r, Named("LL") = (double)ll, Named("niter")=iter, Named("prior_part") = (double)prior_part,
									Named("err")=stop, Named("maxdif_A")=maxdif_A,Named("maxdif_b")=maxdif_b,
									Named("ll_history") = h_ll); 
}

