#include <RcppArmadillo.h>
#include "minimize.h"
#include "shared.h"
#include "pl2_item.h"

using namespace arma;
using Rcpp::Named;



mat full_posterior_2pl(field<mat>& itrace, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
						const vec& theta, const vec& mu, const vec& sigma, const ivec& pgroup)
{
	const int nt = theta.n_elem, np = pni.n_elem, ng = mu.n_elem;
		
	mat posterior0(nt,ng), posterior(nt,np);
	for(int g=0; g<ng; g++)
		posterior0.col(g) = gaussian_pts(mu[g],sigma[g],theta);	

#pragma omp parallel for
	for(int p=0; p<np;p++)
	{
		posterior.col(p) = posterior0.col(pgroup[p]);
		
		for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
			posterior.col(p) %= itrace(pi[indx]).col(px[indx]);
	}
	return posterior;
}


// move to data.cpp
// ? remove dependence on dat with sorted mergejoin, but need ix for that (ix is just a flat dat without NA's)
void persons_ii(const int item1, const int item2, const imat& dat,
				const ivec& inp, const ivec& icnp, const ivec& ip,
				ivec& persons, int& np)
{
	np=0;
	for(int pp=icnp[item1]; pp<icnp[item1+1]; pp++)
	{
		const int p = ip[pp];
		if(dat.at(p,item2) >=0)
		{
			persons[np++] = p;
		}
	}		
}




void pl2_icc(const vec& theta, const ivec& a, const double A, const vec& b, const int ncat, 
				mat& itrace, double* nc_ptr, double* nca_ptr, double* ncab_ptr)
{
	const int nt = theta.n_elem;
	vec p(ncat);
	p[0] = 1;
	
	vec norm_const(nc_ptr,ncat,false,true), norm_const_a(nca_ptr,ncat,false,true), norm_const_ab(ncab_ptr,ncat,false,true);
	
	for(int t=0; t<nt; t++)
	{
		double s=1,sa=0,sab=0;		
		for(int k=1; k<ncat; k++)
		{
			p[k] = std::exp(A*a[k]*(theta[t]-b[k])); 
			s += p[k];
			sa += p[k]*a[k];
			sab += p[k]*a[k]*b[k];
		}

		for(int k=0; k<ncat; k++)
			itrace.at(t,k) = p[k]/s;			
		norm_const[t] = s;
		norm_const_a[t] = sa;
		norm_const_ab[t] = sab;
	}
}


// [[Rcpp::export]]
arma::mat full_hessian_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta, const arma::ivec& item_fixed,
						const arma::imat& dat, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, const arma::ivec& pgroup, const arma::ivec& gn,
						const arma::ivec& ip, const arma::ivec& inp, const arma::ivec& icnp,
						const arma::vec& mu, const arma::vec& sigma, const int ref_group,
						const arma::imat dsg_ii)
{
	const int ng = mu.n_elem, nit=a.n_cols, max_cat=ncat.max(), nt=theta.n_elem, np=pni.n_elem;
	
	int npar = 2*ng;
	if(ref_group >= 0) npar -= 2;
	for(int i=0; i<nit; i++) 
		if(item_fixed[i]==0)
			npar += ncat[i];	
			
	const vec theta2 = square(theta), sigma2 = square(sigma);
	
	field<mat> itrace(nit), itrace2(nit);
	mat nconst(nt,nit), nconst_a(nt,nit), nconst_ab(nt,nit);
	vec P2(nt);
	for(int i=0; i<nit; i++)
	{
		// pe compute necessary traces
		itrace(i) = mat(nt,ncat[i]);		
		pl2_icc(theta, a.col(i), A[i], b.col(i), ncat[i], itrace(i), nconst.colptr(i), nconst_a.colptr(i), nconst_ab.colptr(i));
		itrace2(i) = square(itrace(i));		
	}
		
	mat hess(npar, npar, fill::zeros);
	mat posterior = full_posterior_2pl(itrace, pni, pcni, pi, px, theta, mu, sigma, pgroup);
	vec sum_posterior = trans(sum(posterior,0));
	/*
	for(int p=0; p<np;p++)
	{
		posterior.col(p) = posterior.col(p)/sum_posterior[p];
		
	}
	sum_posterior.ones();
	*/
	int pr=0;

	// declare parallel
	ivec persons_ij(np);
	for(int p=0; p<np;p++)	// complete design for testing
		persons_ij[p]=p;
	
	std::vector<long double> AA_part(3), Ab(max_cat), bA(max_cat); 
	mat Ab_part(3,max_cat), bA_part(3,max_cat), bb(max_cat,max_cat);
	cube bb_part(max_cat,max_cat,3);
	
	

	int np_ij=np;// complete design for testing


	for(int i=0; i<nit; i++) if(item_fixed[i]==0)
	{				
		// block diagonal
		long double AA = 0;
		arma::mat atb(nt,ncat[i],arma::fill::zeros);
		for(int k=1; k<ncat[i];k++)
			atb.col(k) = a.at(k,i) * (theta-b.at(k,i)); 
		arma::vec D = -arma::sum(atb % itrace(i),1);
		arma::vec E = arma::sum(arma::square(atb) % itrace(i),1);
		
		for(int ii=icnp[i]; ii < icnp[i+1]; ii++)
		{
			const int p=ip[ii];
			
			const int x=dat.at(p,i);
			int x1=x;
			double dnm = sum_posterior[p];
			double part1 = arma::accu((arma::square(atb.col(x)) + 2*D%atb.col(x) + 2*arma::square(D) - E) % posterior.col(p));
			double part2 = SQR(arma::accu((atb.col(x) + D) % posterior.col(p)))/dnm;
	
			AA -= (part1-part2)/dnm;	
			
			for(int k=1;k<ncat[i];k++)
			{
				for(int l=k; l<ncat[i]; l++)
				{
					double bb0 = accu(-(kron(l, x1)*a.at(x1,i) - a.at(l,i)*itrace(i).col(l)) % posterior.col(p));
					double bb1 = accu(-(kron(k, x1)*a.at(x1,i) - itrace(i).col(k)*a.at(k,i)) % posterior.col(p));
					double bb2 = accu(posterior.col(p) % (kron(l, x1)*kron(k, x1)*SQR(a.at(x1,i)) - kron(l, x1)*a.at(k,i)*a.at(x1,i)*itrace(i).col(k) \
															- kron(k, l)*SQR(a.at(k,i))*itrace(i).col(k) \
															- kron(k, x1)*a.at(x1,i)*a.at(l,i)*itrace(i).col(l) + 2*a.at(k,i)*a.at(l,i)*itrace(i).col(l)%itrace(i).col(k)));
					
					hess.at(pr+k,pr+l) += SQR(A[i]) * (bb2/sum_posterior[p] - (bb0 * bb1)/SQR(sum_posterior[p]));
				}

				//ab
				
				// this one seems correct but is still the weakest in precision compared with num hess
				double Ab0 = accu((-kron(k, x1)*a.at(x1,i) + itrace(i).col(k)*a.at(k,i)) % posterior.col(p));

				double Ab1 = accu( ((theta - b.at(x1,i))*a.at(x1,i) + ( nconst_ab.col(i) - nconst_a.col(i)%theta)/nconst.col(i)) % posterior.col(p));

				//double Ab2 = accu(((kron(k, x1)*SQR(a.at(x1,i))*b.at(x1,i) - kron(k, x1)*SQR(a.at(x1,i))*theta - a.at(x1,i)*a.at(k,i)*b.at(x1,i)*itrace(i).col(k) + a.at(x1,i)*a.at(k,i)*itrace(i).col(k)%theta - SQR(a.at(k,i))*b.at(k,i)*itrace(i).col(k) + SQR(a.at(k,i))*itrace(i).col(k)%theta)%nconst.col(i) - kron(k, x1)*nconst_ab.col(i)*a.at(x1,i) + kron(k, x1)*a.at(x1,i)*theta%nconst_a.col(i) + 2*nconst_ab.col(i)*a.at(k,i)%itrace(i).col(k) - 2*a.at(k,i)*itrace(i).col(k)%theta%nconst_a.col(i))% posterior.col(p)/nconst.col(i));
				double Ab2 = accu(((kron(k, x1)*SQR(a.at(x1,i))*b.at(x1,i) - kron(k, x1)*SQR(a.at(x1,i))*theta - a.at(x1,i)*a.at(k,i)*b.at(x1,i)*itrace(i).col(k) + a.at(x1,i)*a.at(k,i)*itrace(i).col(k)%theta - SQR(a.at(k,i))*b.at(k,i)*itrace(i).col(k) + SQR(a.at(k,i))*itrace(i).col(k)%theta) \
									 + ( -kron(k, x1)*nconst_ab.col(i)*a.at(x1,i) + kron(k, x1)*a.at(x1,i)*theta%nconst_a.col(i) + 2*nconst_ab.col(i)*a.at(k,i)%itrace(i).col(k) - 2*a.at(k,i)*itrace(i).col(k)%theta%nconst_a.col(i))/nconst.col(i) \
									) % posterior.col(p));
				
				
				hess.at(pr,pr+k) -= A[i] * (Ab0*Ab1/SQR(sum_posterior[p]) -  (Ab2+Ab0)/sum_posterior[p]); // unsure about Ab0 ins he sum
	

			}
		}		
		hess.at(pr,pr) = -AA;
		
		

		// -------------------------------  off diagonal ----------------------------- //
		int qr = pr + ncat[i];
		for(int j=i+1; j<nit; j++) if(item_fixed[i]==0)
		{
			if(dsg_ii.at(j,i) == 1)
			{
				//persons_ii(i,j, dat, inp, icnp, ip, persons_ij, np_ij);

				AA = 0;
				bb.zeros();
				std::fill(Ab.begin(), Ab.end(), .0L);
				std::fill(bA.begin(), bA.end(), .0L);
				for(int pp=0;pp<np_ij; pp++)
				{
					const int p=persons_ij[pp];
					const int x1=dat.at(p,i), x2=dat.at(p,j);
					
					std::fill(AA_part.begin(), AA_part.end(), .0L);
					Ab_part.zeros();
					bA_part.zeros();
					bb_part.zeros();
					for(int t=0; t<nt; t++)
					{
						AA_part[0] += -(nconst_a.at(t,i)*theta[t] + b.at(x1,i)*nconst.at(t,i)*a.at(x1,i) - theta[t]*nconst.at(t,i)*a.at(x1,i) - nconst_ab.at(t,i))*posterior.at(t,p)/nconst.at(t,i);
						AA_part[1] += (theta[t]*nconst.at(t,i)*a.at(x1,i) - theta[t]*nconst_a.at(t,i) - nconst.at(t,i)*a.at(x1,i)*b.at(x1,i) + nconst_ab.at(t,i))*(theta[t]*a.at(x2,j)*nconst.at(t,j) \
									- theta[t]*nconst_a.at(t,j) - b.at(x2,j)*a.at(x2,j)*nconst.at(t,j) + nconst_ab.at(t,j))*posterior.at(t,p)/(nconst.at(t,i)*nconst.at(t,j));
						AA_part[2] += (a.at(x2,j)*theta[t]*nconst.at(t,j) - a.at(x2,j)*nconst.at(t,j)*b.at(x2,j) - theta[t]*nconst_a.at(t,j) + nconst_ab.at(t,j))*posterior.at(t,p)/nconst.at(t,j);
						
						for(int l=1;l<ncat[j];l++)
						{
							
							Ab_part.at(0,l) += -(kron(l, x2)*posterior.at(t,p)*a.at(x2,j) - posterior.at(t,p)*itrace(j).at(t,l)*a.at(l,j))*(a.at(x1,i)*theta[t]*nconst.at(t,i) - a.at(x1,i)*nconst.at(t,i)*b.at(x1,i) \
												- theta[t]*nconst_a.at(t,i)	+ nconst_ab.at(t,i))/nconst.at(t,i);
							Ab_part.at(1,l) += -(b.at(x1,i)*a.at(x1,i)*nconst.at(t,i) - a.at(x1,i)*theta[t]*nconst.at(t,i) + theta[t]*nconst_a.at(t,i) - nconst_ab.at(t,i))*posterior.at(t,p)/nconst.at(t,i);
							Ab_part.at(2,l) += -kron(l, x2)*a.at(x2,j)*posterior.at(t,p) + posterior.at(t,p)*itrace(j).at(t,l)*a.at(l,j);
							
						}
						for(int k=1;k<ncat[j];k++)
						{
							bA_part.at(0,k) += -(kron(k,x1)*posterior.at(t,p)*a.at(x1,i) - posterior.at(t,p)*itrace(i).at(t,k)*a.at(k,i))*(a.at(x2,j)*theta[t]*nconst.at(t,j) - a.at(x2,j)*nconst.at(t,j)*b.at(x2,j) - theta[t]*nconst_a.at(t,j)	+ nconst_ab.at(t,j))/nconst.at(t,j);
							bA_part.at(1,k) += -(b.at(x2,j)*a.at(x2,j)*nconst.at(t,j) - a.at(x2,j)*theta[t]*nconst.at(t,j) + theta[t]*nconst_a.at(t,j) - nconst_ab.at(t,j))*posterior.at(t,p)/nconst.at(t,j);
							bA_part.at(2,k) += -kron(k,x1)*a.at(x1,i)*posterior.at(t,p) + posterior.at(t,p)*itrace(i).at(t,k)*a.at(k,i);						
						}
						for(int k=1; k<ncat[i]; k++)
						{
							for(int l=1;l<ncat[j]; l++)
							{
								bb_part.at(k,l,0) += -(kron(l, x2)*a.at(x2,j) - itrace(j).at(t,l)*a.at(l,j))*posterior.at(t,p);
								bb_part.at(k,l,2) += -(kron(k, x1)*a.at(x1,i) - a.at(k,i)*itrace(i).at(t,k))*posterior.at(t,p);
								bb_part.at(k,l,1) += kron(l, x2)*kron(k, x1)*a.at(x1,i)*posterior.at(t,p)*a.at(x2,j) - kron(l, x2)*posterior.at(t,p)*itrace(i).at(t,k)*a.at(x2,j)*a.at(k,i) \
													- kron(k, x1)*posterior.at(t,p)*a.at(x1,i)*itrace(j).at(t,l)*a.at(l,j) + posterior.at(t,p)*itrace(i).at(t,k)*itrace(j).at(t,l)*a.at(l,j)*a.at(k,i);
								
							}
						}
						
					}
					
					// AA ~ r .99 rechte lijn
					AA += AA_part[1]/sum_posterior[p] - (AA_part[2]*AA_part[0])/SQR(sum_posterior[p]); 
					//Ab ~ r. 1
					for(int l=1;l<ncat[j];l++) 
						Ab[l] += Ab_part.at(0,l)/sum_posterior[p] - (Ab_part.at(1,l) * Ab_part.at(2,l))/SQR(sum_posterior[p]); 
					for(int k=1; k<ncat[i]; k++)
						bA[k] += bA_part.at(0,k)/sum_posterior[p] - (bA_part.at(1,k) * bA_part.at(2,k))/SQR(sum_posterior[p]);

					//bb ~r .99, vermeignvulding met ca 1.1 off
					for(int k=1; k<ncat[i]; k++)
					{
						for(int l=1;l<ncat[j]; l++)
						{
							bb.at(k,l) += bb_part.at(k,l,1)/sum_posterior[p] - (bb_part.at(k,l,2) * bb_part.at(k,l,0))/SQR(sum_posterior[p]);
						}
					}

				}
				//fill
				hess.at(pr,qr) = AA;
				
				for(int k=1; k<ncat[i]; k++)
					hess.at(pr+k,qr)	= bA[k] * A[i]; 		
				for(int l=1;l<ncat[j];l++) 
					hess.at(pr,qr+l) = Ab[l] * A[j]; 			
				
				for(int k=1;k<ncat[i];k++)
					for(int l=1; l<ncat[j]; l++)
						hess.at(pr+k,qr+l) = A[i]*A[j]*bb.at(k,l);				
			}	
			qr += ncat[j];
		}
		pr += ncat[i];
	}
	// pop
	if(ng>1 || ref_group != 0)
	{
		//blockdiagonal
		// running sums
		vec dmu(ng,fill::zeros), dsig(ng,fill::zeros), dmusig(ng,fill::zeros);
		
		// diagonal parts independent of posterior
		mat th_mu2(nt,ng), sigc0(nt,ng), sigc1(nt,ng), muc0(nt,ng), muc1(nt,ng), msc0(nt,ng), msc1(nt,ng) ;
		
		for(int g=0; g<ng; g++) if(g != ref_group)
		{
			
			double dnm = 2*sigma2[g];
			double m0 = accu(exp(-square(theta-mu[g])/dnm));
			double m1 = accu((mu-theta) % exp(-square(theta-mu[g])/dnm))/m0;			
			double m2 = accu((square((mu[g] - theta)/sigma) - 1.0) % exp(-square(theta-mu[g])/dnm))/m0;
			
			th_mu2.col(g) = square(theta-mu[g]);
			
			muc0.col(g) = -m2 - 1 + th_mu2.col(g)/sigma2[g] - 2*m1*(mu[g] - theta)/sigma2[g] + 2*SQR(m1)/sigma2[g];
			muc1.col(g) = theta-mu[g] + m1; 
			
			double c0 = accu(th_mu2.col(g) % exp(-th_mu2.col(g)/(2*sigma2[g]))) / (m0 * sigma[g]);//missch foutje
			double c1 = accu((-3 + th_mu2.col(g)/sigma2[g]) * th_mu2.col(g) * exp(-th_mu2.col(g)/(2*sigma2[g]))) / m0;
			
			sigc0.col(g) = -3*th_mu2.col(g) - c1 + square(th_mu2.col(g)/sigma[g]) - 2*th_mu2.col(g) * c0 + 2*SQR(c0);
			sigc1.col(g) = th_mu2.col(g) - c0;
			
			double ms0 = accu((2 - th_mu2.col(g)/sigma2[g]) % (mu - theta) % exp(-th_mu2.col(g)/(2*sigma2[g])))/m0;
			double ms1 = accu((mu - theta) % exp(-th_mu2.col(g)/(2*sigma2[g])))/m0;
			
			msc0.col(g) = theta + ms1-mu[g];
			msc1.col(g) = (2*mu[g] - 2*theta - ms0 - pow((mu[g] - theta),3)/sigma2[g] + th_mu2.col(g)*ms1/sigma2[g] + (mu[g] - theta)*c0/sigma[g] - 2*ms1*c0/sigma[g]);

		}
		
		for(int p=0; p<np; p++)
		{
			const int g = pgroup[p];
			if(g!=ref_group)
			{
				dmu[g] +=  accu(muc0.col(g) % posterior.col(p)) / sum_posterior[p];				
				dmu[g] -= SQR(accu(muc1.col(g) * posterior.col(p))/(sigma[g]*sum_posterior[p]));
				
				dsig[g] += accu(sigc0.col(g) % posterior.col(p))/sum_posterior[p];
				dsig[g] -= SQR(accu(sigc1.col(g)  % posterior.col(p))/sum_posterior[p]);
				
				dmusig[g] += accu(msc1.col(g) % posterior.col(p))/sum_posterior[p];
				dmusig[g] -= (accu(sigc1.col(g) % posterior.col(p)) * accu(msc0.col(g) % posterior.col(p)))/SQR(sigma[g]*sum_posterior[p]);	
			}
		}
		//fill
		for(int g=0; g<ng; g++) if(g != ref_group)
		{
			hess.at(pr,pr) = dmu[g]/sigma2[g];
			hess.at(pr+1,pr+1) = dsig[g]/SQR(sigma2[g]);
			hess.at(pr,pr+1) = dmusig[g]/CUB(sigma[g]);
			pr += 2;
		}
		
		// off diagonal
		// pop*pop is missing data since sadly no-one can be in two populations simultaneously
		/*
		for(int g=0; g<ng; g++) if(g != ref_group)
		{
			int qr=0;
			for(int i=0; i<nit; i++) if(item_fixed[i]==0)
			{
				if(dsg_gi.at(g,i) == 1)
				{
					
				
				}
				qr += ncat[i];
			}
		}
		*/
		
	}

	return hess;
}

// [[Rcpp::export]]
arma::vec gradient_2pl(arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat,
						const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup, 
						const arma::ivec ip, const arma::ivec& inp, const arma::ivec& icnp)
{
	const int max_cat = ncat.max();
	const int nit =ncat.n_elem;
	vec grd(accu(ncat),fill::zeros);
	
	field<mat> itrace(nit);
	for(int i=0; i<nit; i++)
	{
		itrace(i) = pl2_trace(theta, a.col(i), A[i], b.col(i), ncat[i]);
	}
	
	vec g(max_cat);
	int pr=0;
	for(int i=0; i<nit; i++)
	{	
		ll_pl2_v2 f(itrace, theta, ip, pi, pcni, px, 
							pgroup, inp, icnp, mu, sigma, i, a.col(i));
		
		vec pars = b.col(i).head(ncat[i]);
		pars[0] = A[i];
		f.df(pars,g);
		for(int p=0;p<ncat[i];p++)
			grd[pr++] = g[p];
	}
	return grd;
}

