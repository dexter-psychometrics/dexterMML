#include <stack>
#include <unordered_map>
#include <RcppArmadillo.h>

#include "shared.h"

using namespace arma;
using Rcpp::Named;


int count_not_NA(const imat& dat)
{
	const int sz = dat.n_cols*dat.n_rows;
	const int r4 = sz % 4;
	int nn=0;
	for(int i=0;i<r4;i++)
		nn += dat[i] >= 0;
#pragma omp parallel for reduction(+: nn)
	for(int i=r4;i<sz;i+=4)
	{
		nn += dat[i] >= 0;
		nn += dat[i+1] >= 0;
		nn += dat[i+2] >= 0;
		nn += dat[i+3] >= 0;	
	}
	return nn;
}


// [[Rcpp::export]]
Rcpp::List mat_pre(const arma::imat& dat, const int max_score, const arma::ivec& pgroup, const int ng)
{
	const int nit = dat.n_cols, np = dat.n_rows;
	
	icube icatg(max_score+1,nit, ng, fill::zeros);
	ivec inp(nit,fill::zeros);
	ivec psum(np,fill::zeros), pni(np,fill::zeros);
	
	const int sz = count_not_NA(dat);
	ivec pi(sz),px(sz);

	int pp=0;
	for(int p=0; p<np; p++)
	{
		for(int i=0; i<nit; i++)
		{
			if(dat.at(p,i) >=0)
			{
				int rsp = dat.at(p,i);
				psum[p] += rsp;				
				inp[i]++;
				pni[p]++;
				icatg.at(rsp,i,pgroup[p])++;	
				pi[pp] = i;
				px[pp++] = rsp; 
			}
		}
	}
	ivec ip(sz), ix(sz);
	int ii=0;
	for(int i=0;i<nit;i++)
		for(int p=0;p<np;p++)
			if(dat.at(p,i) >=0)
			{
				ix[ii] = dat.at(p,i);
				ip[ii++] = p;
			}
	// cumulative pointers	
	ivec pcni(np+1);
	pcni[0] = 0;
	std::partial_sum(pni.begin(),pni.end(),pcni.begin()+1);

	ivec icnp(nit+1);
	icnp[0] = 0;
	std::partial_sum(inp.begin(),inp.end(),icnp.begin()+1);

	imat icat = sum(icatg,2);
	ivec imax(nit,fill::zeros), isum(nit,fill::zeros), ncat(nit, fill::zeros);
	for(int i=0;i<nit;i++)
	{
		for(int k = max_score; k>=0; k--)
			if(icat.at(k,i)>0)
			{
				imax[i] = k;
				break;
			}
		for(int k=0; k<=max_score; k++) if(icat.at(k,i)>0)
		{
			isum[i] += k * icat.at(k,i);
			ncat[i]++;
		}
	}	
	
	return Rcpp::List::create(
		Named("pi") = pi, Named("px") = px, Named("ip") = ip, Named("ix") = ix,
		Named("inp") = inp, Named("icnp") = icnp, Named("pni") = pni, Named("pcni") = pcni,
		Named("icat") = icat, Named("icatg") = icatg, Named("ncat") = ncat, Named("imax") = imax, Named("isum") = isum, Named("psum") = psum);
}

// less memory grabbing and faster than going via dexter's get_resp_matrix
// and now just as sorted so in all ways equivalent

// input all 1-indexed, except pgroup & item_score, 0 indexed
// output all 0-indexed
// needs 4 loops over the data, O(n);

// [[Rcpp::export]]
Rcpp::List df_pre(const arma::ivec& person_id, const arma::ivec& pgroup, const int ng, const arma::ivec& item_id, const arma::ivec& item_score, 
					const int max_score, const int np, const int nit, const bool sorted=true)
{
	const int sz = person_id.n_elem;
	
	icube icatg(max_score+1,nit,ng, fill::zeros);	
	
	ivec inp(nit,fill::zeros);
	ivec psum(np,fill::zeros), pni(np,fill::zeros);
	ivec ip(sz), ix(sz);
	ivec pi(sz),px(sz);

	// need counts first
	for(int i=0; i<sz; i++)
	{
		const int p = person_id[i]-1;
		inp[item_id[i]-1]++;
		pni[p]++;	
		psum[p] += item_score[i];
		//icat.at(item_score[i], item_id[i]-1)++;
		icatg.at(item_score[i], item_id[i]-1, pgroup[p])++;
	}
	
	// cumulative pointers	
	ivec pcni(np+1);
	pcni[0] = 0;
	std::partial_sum(pni.begin(),pni.end(),pcni.begin()+1);

	ivec icnp(nit+1);
	icnp[0] = 0;
	std::partial_sum(inp.begin(),inp.end(),icnp.begin()+1);
	
	// effectively a variation of a radix sort
	
	// working copy
	ivec icnp2=icnp;
		
	// fill items
	for(int i=0; i<sz; i++)
	{
		int ii = icnp2[item_id[i]-1]++;
			
		ip[ii] = person_id[i]-1;
		ix[ii] = item_score[i];		
	}
	
	// fill persons, items become sorted within person
	ivec pcni2 = pcni;	
	for(int i=0; i<nit; i++)
	{
		for(int indx=icnp[i]; indx<icnp[i+1]; indx++)
		{
			int pp = pcni2[ip[indx]]++;
			pi[pp] = i;
			px[pp] = ix[indx];
		}
	}
	if(sorted)
	{
		// overwrite items, persons become sorted within items completing the radix for item, person
		icnp2 = icnp;
		for(int p=0; p<np; p++)
		{
			for(int indx=pcni[p]; indx<pcni[p+1]; indx++)
			{
				int ii = icnp2[pi[indx]]++;
				ip[ii] = p;
				ix[ii] = px[indx];
			}
		}
	}
	imat icat = sum(icatg,2);	

	ivec imax(nit,fill::zeros), isum(nit,fill::zeros), ncat(nit, fill::zeros);
	for(int i=0;i<nit;i++)
	{
		for(int k = max_score; k>=0; k--)
			if(icat.at(k,i)>0)
			{
				imax[i] = k;
				break;
			}
		for(int k=0; k<=max_score; k++) if(icat.at(k,i)>0)
		{
			isum[i] += k * icat.at(k,i);
			ncat[i]++;
		}
	}
	
	return Rcpp::List::create(
		Named("pi") = pi, Named("px") = px, Named("ip") = ip, Named("ix") = ix,
		Named("inp") = inp, Named("icnp") = icnp, Named("pni") = pni, Named("pcni") = pcni,
		Named("icat") = icat, Named("icatg") = icatg, Named("ncat") = ncat, Named("imax") = imax, Named("isum") = isum, Named("psum") = psum);
}

// [[Rcpp::export]]
bool duplicate_person_item(const arma::ivec& ip, const arma::ivec& icnp)
{
	const int nit = icnp.n_elem;
	int cnt = 0;
	
#pragma omp parallel for reduction(+: cnt)
	for(int i=0; i<nit; i++)
	{
		for(int j=icnp[i]; j<icnp[i+1]-1; j++)
		{
			cnt += (ip[j]==ip[j+1]);
		}	
	}
	return cnt > 0;
}



// px changed in place if necessary
// returns matrix a
// [[Rcpp::export]]
arma::imat categorize(const arma::ivec& pni,
						const arma::ivec& pcni,
						const arma::ivec& icnp,
						const arma::ivec& pi,				
						const arma::imat& icat, const arma::ivec& imax,
						const int max_cat,
						arma::ivec& px,
						arma::ivec& ix)
{
	
	const int nit = icat.n_cols, np = pni.n_elem;
	imat a(max_cat, nit, fill::zeros);
	imat ai(icat.n_rows, icat.n_cols,fill::zeros);
	
	arma::ivec recode_items(nit);
	int nir=0;
	for(int i=0; i<nit; i++)
	{
		int k=0;
		for(int j=0; j<=imax[i]; j++)
		{
			if(icat.at(j,i) == 0) 
			{
				if(nir==0 || recode_items[nir-1] != i)
					recode_items[nir++] = i;
			}
			else
			{
				a.at(k,i) = j;
				ai.at(j,i) = k;
				k++;
			}
		}			
	}
	
	if(nir > 0)
	{	
#pragma omp parallel for	
		for(int p=0; p<np; p++)
		{
			for(int j=pcni[p]; j<pcni[p+1]; j++)
				px[j] = ai.at(px[j],pi[j]);
		}
		
#pragma omp parallel for	
		for(int ii=0; ii<nir; ii++)
		{
			const int i = recode_items[ii];
			for(int j=icnp[i]; j<icnp[i+1]; j++)
			{
				ix[j] = ai.at(ix[j],i);
			}
		}
		
	}	
	return a;
}

// ai, this one assumes ordered person id's within items
// 2 col datamatrix for any 2 items and persons who did both
void persons_ii(const int item1, const int item2, const ivec& ix,
				const ivec& inp, const ivec& icnp, const ivec& ip,
				ivec& persons, ivec& x1, ivec& x2, int& np)
{
	np=0;
	int pp1=icnp[item1];
	int pp2=icnp[item2];
	while (pp1 < icnp[item1+1] && pp2 <icnp[item2+1])
	{
		if( ip[pp1] == ip[pp2])
		{
			//add output
			x1[np] = ix[pp1];
			x2[np] = ix[pp2];
			persons[np++] = ip[pp1];
			pp1++;
			pp2++;
		}
		else if(ip[pp1] < ip[pp2])
			pp1++;
		else
			pp2++;
	}	
}



/****************************
* Designs
*****************************/


// [[Rcpp::export]]
Rcpp::List design_matrices(const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& pg, const int nit, const int ng)
{
	const int np = pni.n_elem;
	imat item(nit,nit,fill::zeros), group(ng,nit,fill::zeros);
	
	//case complete data
	if(np * nit == pi.n_elem)
	{
		item.ones();
		group.ones();		
	}
	else
	{	
#pragma omp parallel for reduction(||: item, group)
		for(int p=0; p<np; p++)
		{
			const int g = pg[p];
			const int r2 = pni[p] % 2;
			if(r2>0)
			{
				const int itm=pi[pcni[p]];
				group.at(g,itm) = 1;
				for(int j=pcni[p]+1;j<pcni[p+1]; j++) item.at(itm,pi[j]) = 1;	
			}		
			
			for(int i=pcni[p]+r2; i<pcni[p+1]; i+=2)
			{
				const int itm1=pi[i], itm2=pi[i+1];
				
				group.at(g,itm1) = 1;
				group.at(g,itm2) = 1;
				item.at(itm1,itm2) = 1;
				
				for(int j=i+2; j<pcni[p+1]; j++)
				{
					item.at(itm1,pi[j]) = 1;	
					item.at(itm2,pi[j]) = 1;
				}
			}
		}

		item += item.t();
		item.diag().ones(); 
	}
	
	return Rcpp::List::create(Named("items")=item, Named("groups")=group);
}




// bitflag, 0=unidentified; 1=connected; 2,4,8 and combinations are tenuously identified (should give warnings)
// [[Rcpp::export]]
int check_connected_c(const arma::imat& item, const arma::imat& group, const arma::ivec& item_fixed)
{
	const int nit = item.n_cols, ng=group.n_rows;

	// items
	arma::ivec item_groups(nit);
	item_groups.fill(-1);
	int nig=-1, out=0;
	std::stack<int> st;

	for(int j=0;j<nit;j++)
	{
		if(item_groups[j]<0)
		{
			st.push(j);
			item_groups[j] = ++nig;
			while(!st.empty())
			{
				int s = st.top();
				st.pop();
				for(int i=0;i<nit;i++)
				{
					if(item.at(i,s)>0 && item_groups[i]<0)
					{
						item_groups[i]=nig;
						st.push(i);
					}
				}	
			}
		}
	}
	nig++;
	if(nig==1) return 1; // connected via common items
	
	ivec ig_fixed(nig, fill::zeros);
	for(int i=0; i<nit; i++) 
		if(item_fixed[i]==1)
			ig_fixed[item_groups[i]] = 1;
	
	if(all(ig_fixed==1))
		out += 2; // connected via fixed items
	
	if(ng==1)
		out += 4; //connected via groups
	else
	{		
		imat igg(nig,ng, fill::zeros);
		for(int g=0;g<ng;g++)
			for(int i=0; i<nit; i++)
				if(group.at(g,i) == 1)
					igg.at(item_groups[i],g) = 1;

		ivec group_connected(ng, fill::zeros);
		st.push(0);
		group_connected[0]=1;
		while(!st.empty())
		{
			int g = st.top();
			st.pop();		
			for(int ig=0; ig<nig; ig++) if(igg.at(ig,g) == 1)
			{
				for(int g2=0; g2<ng; g2++) if(group_connected[g2]==0 && igg.at(ig,g2)==1) 
				{
					group_connected[g2] = 1;
					st.push(g2);
				}		
			}		
		}

		if(all(group_connected==1))
			out += 4; //connected via groups
		else if(out!=2 && accu(item_fixed)>1)
		{
			group_connected.zeros();
			for(int i=0; i<nit; i++) if(item_fixed[i]==1)
			{
				for(int g=0;g<ng;g++) if(group_connected[g]==0)
				{
					group_connected[g]=1;
					st.push(g);
				}
			}
			while(!st.empty())
			{
				int g = st.top();
				st.pop();		
				for(int ig=0; ig<nig; ig++) if(igg.at(ig,g) == 1)
				{
					for(int g2=0; g2<ng; g2++) if(group_connected[g2]==0 && igg.at(ig,g2)==1) 
					{
						group_connected[g2] = 1;
						st.push(g2);
					}		
				}		
			}
			if(all(group_connected==1))
				out += 8; //all populations that have no item overlap have 1 or more fixed items
		}
	}
	return out;
}


// very unsure about the polytomous correctness of this after categorize a, test with weird a
// when understood, add comments in several places


// heuristic to estimate beta sequentially per population and add to one scale
// not sure if polytomous is done very correclty
// [[Rcpp::export]]
arma::mat start_beta(const arma::imat& a, const arma::ivec& ncat, const arma::icube& icatg, const int ref_group,
					 const arma::ivec& item_fixed, const arma::mat& fixed_beta)
{
	const int nit = a.n_cols, ng=icatg.n_slices;
	const double nc = -std::sqrt(1.702);
	
	imat ign = sum(icatg,0);
	
	if(ign.n_rows == 1) ign = ign.t(); // dimension drops with ng==1 for an unfathomable reason, maybe make a bugreport with armadillo
	
	mat beta(a.n_rows, nit, fill::zeros);
	mat g_beta = beta;
	ivec estimated(nit, fill::zeros), g_estimated(nit, fill::zeros), g_used(ng,fill::zeros);
	
	int max_iter = ng;
	
	if(ref_group < 0) // fixed items
	{		
		for(int i=0; i<nit; i++) if(item_fixed[i] == 1)
			beta.col(i) = fixed_beta.col(i);
		
		estimated = item_fixed;
	}
	else
	{
		int g = ref_group;
		
		for(int i=0; i<nit; i++) if(ign.at(i,g) > 0)
		{
			bool incl = true;
			for(int j=0;j<ncat[i]; j++)
				incl = incl && icatg.at(a.at(j,i),i,g) > 5;
			if(incl)
			{
				for(int j=1;j<ncat[i]; j++)
				{
					beta.at(j,i) = nc * (std::log(icatg.at(a.at(j,i),i,g)) - std::log(icatg.at(0,i,g)))/j;
				}
				estimated[i] = 1;		
			}
		}
		g_used[g] = 1;
		max_iter--;
	}
		
	for(int iter=0; iter<max_iter; iter++)
	{
		int max_ovl = -1, g=-1;

		for(int gg=0; gg<ng; gg++) if(g_used[gg] == 0)
		{
			const int m = accu(ign.col(gg) % estimated);
			if(m > max_ovl)
			{
				max_ovl = m;
				g = gg;
			}
		}
		g_used[g] = 1;
		g_estimated.zeros(); g_beta.zeros();
		for(int i=0; i<nit; i++) if(ign.at(i,g) > 0)
		{			
			bool incl = true;

			for(int j=0;j<ncat[i]; j++)
				incl = incl && icatg.at(a.at(j,i),i,g) > 5;
			
			if(incl)
			{
				for(int j=1;j<ncat[i]; j++)
				{
					g_beta.at(j,i) = nc * (std::log(icatg.at(a.at(j,i),i,g)) - std::log(icatg.at(0,i,g)))/j;
				}
				g_estimated[i] = 1;		
			}
		}
		ivec overlap = g_estimated % estimated;
		double df = 0;
		int nn=0;

		for(int i=0; i<nit; i++) if(overlap[i] == 1)
		{
			for(int j=1;j<ncat[i]; j++)
			{
				df += g_beta.at(j,i) - beta.at(j,i) ;
				nn++;
			}			
		}
		
		g_beta = g_beta - df/nn;
		g_beta.row(0).zeros();
		for(int i=0; i<nit; i++) if(overlap[i] == 0 && g_estimated[i] == 1)
			beta.col(i) = g_beta.col(i);

		estimated = estimated + g_estimated - overlap;			
	}

	return beta;
}
