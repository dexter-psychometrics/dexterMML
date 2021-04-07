#include <RcppArmadillo.h>
#include "data.h"
#include "shared.h"
#include "pl2_item.h"
#include "posterior.h"

using namespace arma;


// to do: see where long doubles are needed


// as yet without the item prior part
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





// [[Rcpp::export]]
arma::mat full_hessian_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta, const arma::ivec& item_fixed,
						const arma::ivec& ix, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, const arma::ivec& pgroup, const arma::ivec& gn,
						const arma::ivec& ip, const arma::ivec& inp, const arma::ivec& icnp,
						const arma::vec& mu, const arma::vec& sigma, const int ref_group,
						const arma::imat dsg_ii, const arma::imat& dsg_gi, 
						const int A_prior=0, const double A_mu=0, const double A_sigma=0.5, const int prog_width=80)
{
	const int ng = mu.n_elem, nit=a.n_cols, max_cat=ncat.max(), nt=theta.n_elem, np=pni.n_elem;
	
	progress_prl progr((nit+ng-1)*(nit+ng)/2, prog_width);
	
	// pre-specify where to save for parallel process
	ivec i_indx(nit+1);
	i_indx[0]=0;
	for(int i=0; i<nit; i++)
	{
		i_indx[i+1] = i_indx[i];
		if(item_fixed[i]==0)
			i_indx[i+1] += ncat[i];
	}
	ivec g_indx(ng+1);	
	g_indx[0]= i_indx[nit];
	for(int g=0; g<ng; g++)
	{		
		g_indx[g+1] = g_indx[g];
		if(g != ref_group)
			g_indx[g+1] +=2;
	}
	const int npar = g_indx[ng];

	const vec sigma2 = square(sigma);
	
	field<mat> itrace(nit);
	field<cube> itrace2(nit);
	mat nconst(nt,nit), nconst_a(nt,nit), nconst_ab(nt,nit);
	
	// pe compute necessary traces
	for(int i=0; i<nit; i++)
	{		
		itrace(i) = mat(nt,ncat[i]);		
		pl2_icc(theta, a.col(i), A[i], b.col(i), ncat[i], itrace(i), nconst.colptr(i), nconst_a.colptr(i), nconst_ab.colptr(i));
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
	mat posterior = normalized_posterior(itrace, pni, pcni, pi, px, theta, mu, sigma, pgroup);

	
	//precompute often recurring constant part
	cube item_const1(nt, max_cat,nit);
	for(int i=0; i<nit; i++) if(item_fixed[i]==0)
	{
		for(int x1=0; x1<ncat[i]; x1++)
		{
			item_const1.slice(i).col(x1) = -(nconst_a.col(i) % theta - nconst_ab.col(i))/nconst.col(i) - a.at(x1,i) * ( b.at(x1,i) - theta);
		}									
	}
	
		
#pragma omp parallel
	{
		const bool _is_main_thread = omp_get_thread_num() == 0;
		int np_ij; 
		ivec persons_ij(np), x1_i(np), x2_j(np);

		std::vector<long double> Ab(max_cat), bA(max_cat); 
		mat bb(max_cat,max_cat);
		
		vec d_bb(max_cat);
		cube d_Ab(nt,max_cat,max_cat);
		mat d_AA(nt,max_cat);
		cube od_AA1(nt,max_cat,max_cat);
#pragma omp for	
		for(int i=0; i<nit; i++) if(item_fixed[i]==0)
		{				
			if(progr.interrupted) continue;
			
			// block diagonal
			const int pr = i_indx[i];
			long double AA = 0;
			arma::mat atb(nt,ncat[i],arma::fill::zeros);
			for(int k=1; k<ncat[i];k++)
				atb.col(k) = a.at(k,i) * (theta-b.at(k,i)); 
			
			for(int x1=0; x1<ncat[i]; x1++)
			{
				d_AA.col(x1) = square(item_const1.slice(i).col(x1)) - arma::sum(arma::square(atb) % itrace(i),1) + square(arma::sum(atb % itrace(i),1));
				for(int k=1;k<ncat[i];k++)
					d_Ab.slice(k).col(x1) = a.at(k,i) * ( kron(k, x1)*(a.at(k,i)*(b.at(k,i) -theta)-1) \
												+ itrace(i).col(k) % (a.at(x1,i)*(theta-b.at(x1,i)) +  a.at(k,i)*(theta-b.at(k,i)) + 1) \
												+ (2*itrace(i).col(k)- kron(k, x1)) % (nconst_ab.col(i) - theta%nconst_a.col(i))/nconst.col(i));
			}
			
			
			for(int ii=icnp[i]; ii < icnp[i+1]; ii++)
			{
				const int p=ip[ii];			
				const int x1=ix[ii]; 
				
				AA += accu( d_AA.col(x1) % posterior.col(p)) - SQR(accu(item_const1.slice(i).col(x1) % posterior.col(p)));
				
				for(int k=1;k<ncat[i];k++)
					d_bb[k] = accu(-(kron(k, x1)*a.at(x1,i) - itrace(i).col(k)*a.at(k,i)) % posterior.col(p));
		
				for(int k=1;k<ncat[i];k++)
				{
					for(int l=k; l<ncat[i]; l++)
					{
						double bb2 = a.at(k,i)*a.at(l,i)*accu(posterior.col(p) % (kron(l, x1)*kron(k, x1) - kron(l, x1)*itrace(i).col(k) - kron(k, l)*itrace(i).col(k) - kron(k, x1)*itrace(i).col(l) + 2*itrace2(i).slice(k).col(l) ));
																
						hess.at(pr+k,pr+l) += SQR(A[i]) * (bb2 - d_bb[k] * d_bb[l]);
					}
					//ab
					hess.at(pr,pr+k) += A[i] * (accu(d_Ab.slice(k).col(x1) % posterior.col(p)) - d_bb[k] * accu(item_const1.slice(i).col(x1) % posterior.col(p)));

				}
			}		
			hess.at(pr,pr) = AA;
			
			// -------------------------------  off diagonal ----------------------------- //

			for(int j=i+1; j<nit; j++) if(item_fixed[i]==0)
			{
				const int qr=i_indx[j];
				if(dsg_ii.at(j,i) == 1)
				{
					persons_ii(i,j, ix, inp, icnp, ip, persons_ij, x1_i, x2_j, np_ij);
					
					AA = 0;
					bb.zeros();
					std::fill(Ab.begin(), Ab.end(), .0L);
					std::fill(bA.begin(), bA.end(), .0L);
					
					for(int x1=0; x1<ncat[i]; x1++)
						for(int x2=0; x2<ncat[j]; x2++)
							od_AA1.slice(x1).col(x2) = item_const1.slice(i).col(x1) % item_const1.slice(j).col(x2);

					for(int pp=0;pp<np_ij; pp++)
					{
						const int p=persons_ij[pp];
						const int x1 = x1_i[pp], x2 = x2_j[pp]; 
				
						// AA ~ r 1
						// a number of sums can be saved as doubles
						double sum1 = accu(item_const1.slice(i).col(x1) % posterior.col(p));
						double sum2 = accu(item_const1.slice(j).col(x2) % posterior.col(p));
						vec sum4(ncat[j]);
						
						AA += accu(od_AA1.slice(x1).col(x2) % posterior.col(p)) - sum1*sum2; 
						
						for(int l=1;l<ncat[j];l++) 
						{
							sum4[l] = accu(a.at(l,j)*(itrace(j).col(l) - kron(l,x2)) % posterior.col(p));
							Ab[l] += accu(a.at(l,j)*(itrace(j).col(l) - kron(l,x2)) % item_const1.slice(i).col(x1) % posterior.col(p)) \
										- sum1 * sum4[l];							
						}					
						
						for(int k=1; k<ncat[i]; k++)
						{
							double sum3 = accu(a.at(k,i)*(itrace(i).col(k) - kron(k,x1)) % posterior.col(p));
							bA[k] += accu(a.at(k,i)*(itrace(i).col(k) - kron(k,x1)) % item_const1.slice(j).col(x2) % posterior.col(p)) - sum2 * sum3;
							
							for(int l=1;l<ncat[j]; l++)
							{
								bb.at(k,l) += a.at(k,i)*a.at(l,j) * (kron(l, x2)*kron(k, x1) + accu(posterior.col(p) % (itrace(i).col(k) % itrace(j).col(l)- kron(l, x2)*itrace(i).col(k) - kron(k, x1)*itrace(j).col(l) )));
								bb.at(k,l) -= sum4[l] * sum3;
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
			}
			progr.update(nit-i, _is_main_thread);			
		}
	}
	if(progr.interrupted) Rcpp::stop("user interruption");
	
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
			if(g!=ref_group)
			{
				dmu[g] +=  accu(muc0.col(g) % posterior.col(p));			
				dmu[g] -= SQR(accu(muc1.col(g) % posterior.col(p))/sigma[g]);
				
				dsig[g] += accu(sigc0.col(g) % posterior.col(p)) \
							- SQR(accu(sigc1.col(g)  % posterior.col(p))/sigma[g]);
				
				dmusig[g] += accu(msc1.col(g) % posterior.col(p));
				dmusig[g] -= (accu(sigc1.col(g) % posterior.col(p)) * accu(msc0.col(g) % posterior.col(p)))/sigma2[g];	
			}
		}
		
		//fill
		for(int g=0; g<ng; g++) if(g != ref_group)
		{
			int pr = g_indx[g];
			hess.at(pr,pr) = dmu[g]/sigma2[g];
			hess.at(pr+1,pr+1) = dsig[g]/SQR(sigma2[g]);
			hess.at(pr,pr+1) = dmusig[g]/CUB(sigma[g]);
		}
		progr.update(ng,true);
		if(progr.interrupted) Rcpp::stop("user interruption");
		// off diagonal
#pragma omp parallel
		{
			const bool _is_main_thread = omp_get_thread_num() == 0;
			field<cube> bmu0(ng), bsig0(ng);
			for(int g=0; g<ng; g++)
			{
				bmu0(g) = cube(nt,max_cat,max_cat);
				bsig0(g) = cube(nt,max_cat,max_cat);
			}
			cube amu0(nt,max_cat,ng), asig0(nt,max_cat,ng), bmu1(nt,max_cat,max_cat);
			mat amu1(nt,max_cat);
			mat amu2(nt,ng),asig2(nt,ng);
			vec damu(ng),dasig(ng);
			mat dbmu(ng,max_cat);
			mat dbsig(ng,max_cat);
#pragma omp for
			for(int i=0; i<nit; i++) if(item_fixed[i]==0)
			{
				if(progr.interrupted) continue;				
				damu.zeros();
				dasig.zeros();
				dbmu.zeros();
				dbsig.zeros();
				
				for(int g=0; g<ng; g++) if(g != ref_group && dsg_gi.at(g,i) == 1)
				{
					double dnm = 2*sigma2[g];
					double m1 = accu((mu[g]-theta) % exp(-square(theta-mu[g])/dnm));	
					amu2.col(g) = -mu[g] + theta + m1; //not dependent on item
					asig2.col(g) = square(mu[g]-theta) - m1;
					for(int x1=0; x1<ncat[i]; x1++)
					{
						amu0.slice(g).col(x1) =  a.at(x1,i) * (mu[g]-theta) % (b.at(x1,i)-theta) \
														- (mu[g]-theta) % (nconst_ab.col(i) - theta % nconst_a.col(i))/nconst.col(i) \
														- a.at(x1,i) * (b.at(x1,i)-theta) * m1 \
														+ m1 * (nconst_ab.col(i) - theta % nconst_a.col(i))/nconst.col(i);
												
						asig0.slice(g).col(x1) = -( a.at(x1,i) * square(mu[g]-theta) % (b.at(x1,i)-theta) \
														- square(mu[g]-theta) % (nconst_ab.col(i) - theta % nconst_a.col(i))/nconst.col(i) \
														- a.at(x1,i) * (b.at(x1,i)-theta) * m1 \
														+ m1 * (nconst_ab.col(i) - theta % nconst_a.col(i))/nconst.col(i));
						
						for(int k=1; k<ncat[i]; k++)
						{
							bmu0(g).slice(x1).col(k) = kron(k,x1) * a.at(x1,i) * (mu[g]-theta) \
														- (mu[g]-theta) % (a.at(k,i) * itrace(i).col(k))\
														- a.at(x1,i) * kron(k,x1) * m1 \
														+ m1 * a.at(k,i) * itrace(i).col(k);
						
							bsig0(g).slice(x1).col(k) = -(kron(k,x1) * a.at(x1,i) * square(mu[g]-theta) \
														- square(mu[g]-theta) % (a.at(k,i) * itrace(i).col(k)) \
														- a.at(x1,i) * kron(k,x1) * m1 \
														+ m1 * a.at(k,i) * itrace(i).col(k));
						}

					}
				}
				for(int x1=0; x1<ncat[i]; x1++)
				{
					amu1.col(x1) = a.at(x1,i) * (theta - b.at(x1,i)) + (nconst_ab.col(i) - theta % nconst_a.col(i))/nconst.col(i);
					for(int k=1; k<ncat[i]; k++)
					{
						bmu1.slice(x1).col(k) = -kron(k,x1) * a.at(x1,i) + a.at(k,i) * itrace(i).col(k);
					}				
				}
				
				
				for(int ii=icnp[i]; ii < icnp[i+1]; ii++) 
				{
					const int p=ip[ii];	
					const int g = pgroup[p];
					if(g!=ref_group)
					{
						const int x1 = ix[ii];	
						damu[g] +=  accu(amu0.slice(g).col(x1) % posterior.col(p)) - (accu(amu1.col(x1) % posterior.col(p)) * accu(amu2.col(g) % posterior.col(p)));
						dasig[g] += accu(asig0.slice(g).col(x1) % posterior.col(p)) - (accu(amu1.col(x1) % posterior.col(p)) * accu(asig2.col(g) % posterior.col(p)));
						for(int k=1; k<ncat[i]; k++)
						{
							dbmu.at(g,k) +=  accu(bmu0(g).slice(x1).col(k) % posterior.col(p)) - (accu(bmu1.slice(x1).col(k) % posterior.col(p)) * accu(amu2.col(g) % posterior.col(p)));
							dbsig.at(g,k) +=  accu(bsig0(g).slice(x1).col(k) % posterior.col(p)) - (accu(bmu1.slice(x1).col(k) % posterior.col(p)) * accu(asig2.col(g) % posterior.col(p)));
						}
					}
				}
				for(int g=0; g<ng; g++) if(g != ref_group && dsg_gi.at(g,i) == 1)
				{
					hess.at(i_indx[i],g_indx[g]) = damu[g]/sigma2[g];
					hess.at(i_indx[i],g_indx[g]+1) = dasig[g]/CUB(sigma[g]);
					for(int k=1; k<ncat[i]; k++)
					{
						hess.at(i_indx[i]+k, g_indx[g]) = A[i] * dbmu.at(g,k)/sigma2[g];
						hess.at(i_indx[i]+k, g_indx[g]+1) = A[i] * dbsig.at(g,k)/CUB(sigma[g]);
					}
				}
				progr.update(ng-1, _is_main_thread);
			}
		}
	}
	if(progr.interrupted) Rcpp::stop("user interruption");
	//prior part
	if(A_prior==1)
	{
		const double asig2 = SQR(A_sigma);
		for(int i=0; i<nit; i++) if(item_fixed[i]==0)
			hess.at(i_indx[i],i_indx[i]) -= (A_mu-asig2 - std::log(A[i]) + 1)/(SQR(A[i])*asig2); 
	}
	else if(A_prior==2)
	{
		const double a_part = 1/SQR(A_sigma);
		for(int i=0; i<nit; i++) if(item_fixed[i]==0)
			hess.at(i_indx[i],i_indx[i]) -= a_part; 
	}
	//lower tri
	for(int i=0;i<npar;i++)
		for(int j=i+1; j<npar;j++)
			hess.at(j,i) = hess.at(i,j);
	
	progr.close();
	return hess;
}


