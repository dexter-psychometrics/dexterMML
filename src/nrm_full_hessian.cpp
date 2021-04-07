#include <RcppArmadillo.h>
#include "data.h"
#include "shared.h"
#include "nrm_item.h"

using namespace arma;


//actually, this is identical to 2pl, prbl goes for estep function as well if supply trace, might simplify some code
mat full_posterior_nrm(field<mat>& itrace, const ivec& pni, const ivec& pcni, const ivec& pi, const ivec& px, 
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
		long double s=0;
		for(int t=0;t<nt;t++)
			s+= posterior.at(t,p);
		posterior.col(p) /= s;
		
	}
	return posterior;
}


// [[Rcpp::export]]
arma::mat full_hessian_nrm(const arma::imat& a, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta, const arma::ivec& item_fixed,
						const arma::ivec& ix, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, const arma::ivec& pgroup, const arma::ivec& gn,
						const arma::ivec& ip, const arma::ivec& inp, const arma::ivec& icnp,
						const arma::vec& mu, const arma::vec& sigma, const int ref_group,
						const arma::imat dsg_ii, const arma::imat& dsg_gi, 
						const int prog_width=80)
{
	const int ng = mu.n_elem, nit=a.n_cols, max_cat=ncat.max(), nt=theta.n_elem, np=pni.n_elem;
	
	progress_prl progr((nit+ng-1)*(nit+ng)/2, prog_width);
	
	ivec ncat1 = ncat-1;
	int npar = 2*ng;
	if(ref_group >= 0) npar--;
	for(int i=0; i<nit; i++) 
		if(item_fixed[i]==0)
			npar += ncat1[i];	
	
	ivec cncat1(nit +1);
	cncat1[0] = 0;
	std::partial_sum(ncat1.begin(),ncat1.end(),cncat1.begin()+1);
	
	ivec g_indx(ng);	
	g_indx[0]= cncat1[nit];
	for(int g=1; g<ng; g++)
	{		
		g_indx[g] = g_indx[g-1]+1;
		if(g-1 != ref_group)
			g_indx[g]++;
	}

	const vec sigma2 = square(sigma);
	
	const int max_a = a.max();
	mat exp_at(max_a+1, nt, fill::ones);
	for(int t=0; t< nt; t++)
		for(int k=1; k<=max_a;k++)
			exp_at.at(k,t) = std::exp(k*theta[t]);
	
	field<mat> itrace(nit);
	field<cube> itrace2(nit);

	// pe compute necessary traces
	for(int i=0; i<nit; i++)
	{		
		itrace(i) = nrm_trace(theta, a.col(i), b.col(i), ncat[i], exp_at);	
		
		itrace2(i) = cube(nt,ncat[i],ncat[i]);
		for(int k=0; k<ncat[i]; k++)
		{
			itrace2(i).slice(k).col(k) = square(itrace(i).col(k));
			for(int l=k+1; l<ncat[i]; l++)
			{
				itrace2(i).slice(k).col(l) = itrace(i).col(k) % itrace(i).col(l);
				itrace2(i).slice(l).col(k) = itrace2(i).slice(k).col(l);
			}
		}
	}

	
	mat hess(npar, npar, fill::zeros);
	mat posterior = full_posterior_nrm(itrace, pni, pcni, pi, px, theta, mu, sigma, pgroup);

#pragma omp parallel
	{
		const bool _is_main_thread = omp_get_thread_num() == 0;
		int np_ij; 
		ivec persons_ij(np), x1_i(np), x2_j(np);

		mat bb(max_cat,max_cat);		
		vec d_bb(max_cat);

#pragma omp for	
		for(int i=0; i<nit; i++) if(item_fixed[i]==0)
		{				
			if(progr.interrupted) continue;

			// block diagonal
			const int pr = cncat1[i];
			
			for(int ii=icnp[i]; ii < icnp[i+1]; ii++)
			{
				const int p = ip[ii];			
				const int x1 = ix[ii]; 
				
				for(int k=1;k<ncat[i];k++)
					d_bb[k] = accu((kron(k, x1)- itrace(i).col(k)) % posterior.col(p));
		
				for(int k=1;k<ncat[i];k++)
				{
					for(int l=k; l<ncat[i]; l++)
					{
						double bb2 = accu(posterior.col(p) % (kron(l, x1)*kron(k, x1) - kron(l, x1)*itrace(i).col(k) - kron(k, l)*itrace(i).col(k) - kron(k, x1)*itrace(i).col(l) + 2*itrace2(i).slice(k).col(l) ));
																
						hess.at(pr+k-1,pr+l-1) += bb2 - d_bb[k] * d_bb[l];
					}
				}
			}		

			// -------------------------------  off diagonal ----------------------------- //

			for(int j=i+1; j<nit; j++) if(item_fixed[i]==0)
			{
				const int qr=cncat1[j];
				if(dsg_ii.at(j,i) == 1)
				{
					persons_ii(i,j, ix, inp, icnp, ip, persons_ij, x1_i, x2_j, np_ij);
					
					bb.zeros();

					for(int pp=0;pp<np_ij; pp++)
					{
						const int p=persons_ij[pp];
						const int x1 = x1_i[pp], x2 = x2_j[pp]; 
				
						vec sumj(ncat[j]);
						
						for(int l=1;l<ncat[j];l++) 
						{
							sumj[l] = accu((itrace(j).col(l) - kron(l,x2)) % posterior.col(p));
						}					
						
						for(int k=1; k<ncat[i]; k++)
						{
							double sumi = accu((itrace(i).col(k) - kron(k,x1)) % posterior.col(p));
							
							for(int l=1;l<ncat[j]; l++)
							{
								bb.at(k,l) += (kron(l, x2)*kron(k, x1) + accu(posterior.col(p) % (itrace(i).col(k) % itrace(j).col(l)- kron(l, x2)*itrace(i).col(k) - kron(k, x1)*itrace(j).col(l) )));
								bb.at(k,l) -= sumj[l] * sumi;
							}						
						}

					}
					//fill
					for(int k=1;k<ncat[i];k++)
						for(int l=1; l<ncat[j]; l++)
							hess.at(pr+k-1,qr+l-1) = bb.at(k,l);				
				}	
			}
			progr.update(nit-i, _is_main_thread);			
		}
	}
	if(progr.interrupted) Rcpp::stop("user interruption");
	
	// pop
	{
		//blockdiagonal
		// running sums
		vec dmu(ng,fill::zeros), dsig(ng,fill::zeros), dmusig(ng,fill::zeros);
		
		// diagonal parts independent of posterior
		mat th_mu2(nt,ng), sigc0(nt,ng), sigc1(nt,ng), muc0(nt,ng), muc1(nt,ng), msc0(nt,ng), msc1(nt,ng) ;
		
		for(int g=0; g<ng; g++)
		{
			
			double dnm = 2*sigma2[g];
			double m0 = accu(exp(-square(theta-mu[g])/dnm));
			double m1 = accu((mu[g]-theta) % exp(-square(theta-mu[g])/dnm))/m0;			
			double m2 = accu((square((mu[g] - theta)/sigma[g]) - 1.0) % exp(-square(theta-mu[g])/dnm))/m0;
			
			th_mu2.col(g) = square(theta-mu[g]);
			
			muc0.col(g) = -m2 - 1 + th_mu2.col(g)/sigma2[g] - 2*m1*(mu[g] - theta)/sigma2[g] + 2*SQR(m1)/sigma2[g];
			muc1.col(g) = theta-mu[g] + m1; 
			
			double c0 = accu(th_mu2.col(g) % exp(-th_mu2.col(g)/(2*sigma2[g]))) / (m0*sigma2[g]); // herschrijven naar ss0
			
			double ss0 = accu(th_mu2.col(g) % exp(-th_mu2.col(g)/(2*sigma2[g])));
			double ss1 = accu((-3 + th_mu2.col(g)/sigma2[g]) % th_mu2.col(g) % exp(-th_mu2.col(g)/(2*sigma2[g])));

			// klopt
			sigc0.col(g) = -3*th_mu2.col(g) - ss1/m0 + square(th_mu2.col(g))/sigma2[g] - 2*th_mu2.col(g)*ss0/(m0*sigma2[g]) + 2*SQR(ss0)/(sigma2[g]*SQR(m0));
			sigc1.col(g) = (th_mu2.col(g) - ss0/m0);  // wordt ook gebruikt in msc

			double ms0 = accu((2 - th_mu2.col(g)/sigma2[g]) % (mu[g] - theta) % exp(-th_mu2.col(g)/(2*sigma2[g])))/m0;
			double ms1 = accu((mu[g] - theta) % exp(-th_mu2.col(g)/(2*sigma2[g])))/m0;
			
			msc0.col(g) = theta + ms1-mu[g];
			msc1.col(g) = (2*mu[g] - 2*theta - ms0 - pow((mu[g] - theta),3)/sigma2[g] + th_mu2.col(g)*ms1/sigma2[g] + (mu[g] - theta)*c0/sigma[g] - 2*ms1*c0/sigma[g]);

		}
#pragma omp parallel for reduction(+:dmu,dsig,dmusig)
		for(int p=0; p<np; p++)
		{
			const int g = pgroup[p];
			
			dsig[g] += accu(sigc0.col(g) % posterior.col(p)) \
							- SQR(accu(sigc1.col(g)  % posterior.col(p))/sigma[g]);
			if(g!=ref_group)
			{
				dmu[g] +=  accu(muc0.col(g) % posterior.col(p));			
				dmu[g] -= SQR(accu(muc1.col(g) % posterior.col(p))/sigma[g]);				
				dmusig[g] += accu(msc1.col(g) % posterior.col(p));
				dmusig[g] -= (accu(sigc1.col(g) % posterior.col(p)) * accu(msc0.col(g) % posterior.col(p)))/sigma2[g];	
			}
		}
		
		//fill
		for(int g=0; g<ng; g++) 
		{
			const int pr = g_indx[g];
			if(g==ref_group)
			{
				hess.at(pr,pr) = dsig[g]/SQR(sigma2[g]);
			} else
			{
				hess.at(pr,pr) = dmu[g]/sigma2[g];
				hess.at(pr+1,pr+1) = dsig[g]/SQR(sigma2[g]);
				hess.at(pr,pr+1) = dmusig[g]/CUB(sigma[g]);
			}			
		}
		progr.update(ng,true);
		if(progr.interrupted) Rcpp::stop("user interruption");
		// off diagonal
#pragma omp parallel
		{
			const bool _is_main_thread = omp_get_thread_num() == 0;
			cube bmu1(nt,max_cat,max_cat);
			field<cube> bmu0(ng), bsig0(ng);
			for(int g=0; g<ng; g++)
			{
				bmu0(g) = cube(nt,max_cat,max_cat);
				bsig0(g) = cube(nt,max_cat,max_cat);
			}
			mat dbmu(ng,max_cat);
			mat dbsig(ng,max_cat);
			mat mu2(nt,ng),sig2(nt,ng);
#pragma omp for
			for(int i=0; i<nit; i++) if(item_fixed[i]==0)
			{
				if(progr.interrupted) continue;				
				dbmu.zeros();
				dbsig.zeros();
				
				for(int g=0; g<ng; g++) if(dsg_gi.at(g,i) == 1)
				{
					double dnm = 2*sigma2[g];
					double m1 = accu((mu[g]-theta) % exp(-square(theta-mu[g])/dnm));	
					mu2.col(g) = -mu[g] + theta + m1; //not dependent on item
					sig2.col(g) = square(mu[g]-theta) - m1;

					for(int x1=0; x1<ncat[i]; x1++)
					{
						for(int k=1; k<ncat[i]; k++)
						{
							bmu0(g).slice(x1).col(k) = kron(k,x1) * (mu[g]-theta) \
														- (mu[g]-theta) %  itrace(i).col(k)\
														-  kron(k,x1) * m1 \
														+ m1 * itrace(i).col(k);
						
							bsig0(g).slice(x1).col(k) = -(kron(k,x1)  * square(mu[g]-theta) \
														- square(mu[g]-theta) % itrace(i).col(k) \
														- kron(k,x1) * m1 \
														+ m1 * itrace(i).col(k));
						}
					}
				}
				for(int x1=0; x1<ncat[i]; x1++)
				{
					for(int k=1; k<ncat[i]; k++)
					{
						bmu1.slice(x1).col(k) = -kron(k,x1)  +  itrace(i).col(k);
					}				
				}				
				
				for(int ii=icnp[i]; ii < icnp[i+1]; ii++) 
				{
					const int p=ip[ii];	
					const int g = pgroup[p];
					const int x1 = ix[ii];	
					for(int k=1; k<ncat[i]; k++)
					{
						dbmu.at(g,k) +=  accu(bmu0(g).slice(x1).col(k) % posterior.col(p)) - (accu(bmu1.slice(x1).col(k) % posterior.col(p)) * accu(mu2.col(g) % posterior.col(p)));
						dbsig.at(g,k) +=  accu(bsig0(g).slice(x1).col(k) % posterior.col(p)) - (accu(bmu1.slice(x1).col(k) % posterior.col(p)) * accu(sig2.col(g) % posterior.col(p)));
					}

				}
				for(int g=0; g<ng; g++) if(dsg_gi.at(g,i) == 1)
				{
					if(g==ref_group)
					{
						for(int k=1; k<ncat[i]; k++)
							hess.at(cncat1[i]+k-1, g_indx[g]) = -dbsig.at(g,k)/CUB(sigma[g]);
					}
					else
					{
						for(int k=1; k<ncat[i]; k++)
						{
							hess.at(cncat1[i]+k-1, g_indx[g]) = -dbmu.at(g,k)/sigma2[g];
							hess.at(cncat1[i]+k-1, g_indx[g]+1) = -dbsig.at(g,k)/CUB(sigma[g]);
						}
					}
				}
				progr.update(ng-1, _is_main_thread);
			}
		}
	}
	if(progr.interrupted) Rcpp::stop("user interruption");
	//lower tri
	for(int i=0;i<npar;i++)
		for(int j=i+1; j<npar;j++)
			hess.at(j,i) = hess.at(i,j);
	
	progr.close();
	return hess;
}


