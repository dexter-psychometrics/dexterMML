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
Rcpp::List estimate_poly2(arma::imat& a, const arma::vec& A_start, const arma::mat& b_start, const arma::ivec& ncat,
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const arma::ivec& item_fixed,
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
	int iter = 0, min_error=0,stop=0;
	long double ll, old_ll=-std::numeric_limits<long double>::max(), prior_part=0;
	double maxdif_A, maxdif_b;
	
	for(; iter<max_iter; iter++)
	{
		estep_poly2(a, A, b, ncat, pni, pcni, pi, px, 
					theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);
		
		if(ll + prior_part < old_ll)
		{
			stop += 2;
			break;
		}
		old_ll = prior_part+ll;
		maxdif_A=0; maxdif_b=0; prior_part=0;
		
#pragma omp parallel for reduction(max: maxdif_A, maxdif_b) reduction(+:min_error, prior_part)
		for(int i=0; i<nit; i++)
		{	
			if(item_fixed[i] == 1)
				continue;
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
		if(min(abs(A)) < 1e-6)
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
									Named("err")=stop, Named("maxdif_A")=maxdif_A,Named("maxdif_b")=maxdif_b); 
}




// [[Rcpp::export]]
Rcpp::List estimate_poly2_test(arma::imat& a, const arma::vec& A_start, const arma::mat& b_start, const arma::ivec& ncat,
						const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, 
						const arma::ivec& item_fixed,
						const int ref_group=0, const int A_prior=0, const double A_mu=0, const double A_sigma=0.5, 
						const int max_iter=200, const int pgw=80)
{
	const int nit = a.n_cols, nt = theta.n_elem, np = pni.n_elem, ng=gn.n_elem;
	const double stop_crit = .0001;
	
	progress_est prog(max_iter, pgw);
	
	mat b = b_start;
	vec A = A_start;
	mat hb(b.n_cols,max_iter), hA(A.n_elem,max_iter);

	field<mat> r(nit);
	for(int i=0; i<nit; i++)
		r(i) = mat(nt,ncat[i]);
	
	vec thetabar(np,fill::zeros);
	
	vec sigma = sigma_start, mu=mu_start;
	
	vec sum_theta(ng), sum_sigma2(ng);
	
	vec h_ll(max_iter,fill::zeros);
	
	cube store_b(b.n_rows, nit, 2, fill::zeros);
	mat store_A(nit,2,fill::zeros);
	int store_i=0;
	
	const double tol = 1e-10;
	int iter = 0, min_error=0,stop=0;
	long double ll, old_ll=-std::numeric_limits<long double>::max();
	double maxdif_A, maxdif_b;
	bool backtrack=false;

	for(; iter<max_iter; iter++)
	{
		hA.col(iter) = A;
		hb.col(iter) = b.row(1).t();
		
		estep_poly2(a, A, b, ncat, pni, pcni, pi, px, 
					theta, r, thetabar, sum_theta, sum_sigma2, mu, sigma, pgroup, ll);
		
		h_ll[iter] = (double)ll;
		
		if(backtrack)
		{
			if(ll < old_ll)
			{
				ll = old_ll;
				A=store_A.col(1-store_i);
				b=store_b.slice(1-store_i);
			} else
			{
				maxdif_A = max(abs(A-store_A.col(1-store_i)));
				maxdif_b = (abs(b-store_b.slice(1-store_i))).max();
			}
			stop+=2;
			break;
		}		
		if(ll < old_ll)
		{
			backtrack=true;
			mat qx = {{1,0,0},{1,1,1},{1,2,4}};
			vec qy = {h_ll[iter-2], h_ll[iter-1], h_ll[iter]};
			vec qs = solve(qx, qy);
			double qpnt = -qs[1]/(2*qs[2]);
			//store_i, 1-store_i, current
			if(qpnt<=1)
			{
				A = store_A.col(1-store_i) * qpnt + (1-qpnt) * store_A.col(store_i);
				b = store_b.slice(1-store_i) * qpnt + (1-qpnt) * store_b.slice(store_i);
			}
			else
			{
				qpnt-=1;
				A =  A * qpnt + (1-qpnt) * store_A.col(1-store_i);
				b = b * qpnt + (1-qpnt) * store_b.slice(1-store_i);			
			}
			continue;
		}
		
		maxdif_A=0; maxdif_b=0;
		store_A.col(store_i)=A;
		store_b.slice(store_i)=b;
		store_i=1-store_i;
		
#pragma omp parallel for reduction(max: maxdif_A, maxdif_b) reduction(+:min_error)
		for(int i=0; i<nit; i++)
		{	
			if(item_fixed[i] == 1)
				continue;
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
				vec pars2 = -b.col(i).head(ncat[i])/5; //div by 5 is very arbitrary
				pars2[0] = -2*A[i]; // times 2 as well, can also try just flip?
				nlm(pars2, tol, itr, ll_itm2, f, err2);	
				if(err2==0 && ll_itm2<ll_itm)
				{
					pars=pars2;
					err=err2;
				}			
			}
			if(std::abs(pars[0])<1e-6)
				err=10;

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
			}
		}
	
		if(min_error>0)
		{
			stop += 1;
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
				
		if(maxdif_b < stop_crit && maxdif_A < stop_crit)
			break;
		
		old_ll = ll;		
	}
	if(iter>=max_iter-1)
		stop += 4;
	
	prog.close();
	
	return Rcpp::List::create(Named("A")=A, Named("b")=b, Named("thetabar") = thetabar, Named("mu") = mu, Named("sd") = sigma, 
									Named("r")=r, Named("LL") = (double)ll, Named("niter")=iter,
									Named("err")=stop, Named("maxdif_A")=maxdif_A,Named("maxdif_b")=maxdif_b,
									Named("h_ll")=h_ll,
									Named("ha")=hA,Named("hb") = hb); 
}


