#include <RcppArmadillo.h>
#include <stack>

using namespace arma;
using Rcpp::Named;



// to do: min should be >=0
// gives two data sets:
// person, item, score
// item, person, score
// to do: overkill ip is not used anywhere
// [[Rcpp::export]]
Rcpp::List mat_pre(arma::imat& dat, const int max_score)
{
	const int nit = dat.n_cols, np = dat.n_rows;
	
	imat icat(max_score+1,nit, fill::zeros);
	ivec inp(nit,fill::zeros);
	ivec psum(np,fill::zeros), pni(np,fill::zeros);

// margins
	for(int i=0; i<nit; i++)
	{
		const ivec rsp(dat.colptr(i), np, false, true);

		for(int p=0; p<np; p++)
		{
			if(rsp[p]>=0) // NA test
			{
				psum[p] += rsp[p];				
				inp[i]++;
				pni[p]++;
				icat.at(rsp[p],i)++;			
			}		
		}	
	}

	// cumulative pointers	
	ivec icnp(nit+1), pcni(np+1);
	icnp[0] = 0;
	pcni[0] = 0;
	std::partial_sum(inp.begin(),inp.end(),icnp.begin()+1);
	std::partial_sum(pni.begin(),pni.end(),pcni.begin()+1);

	// response vectors organized two ways
	const int sz = icnp[nit];
	
	ivec ip(sz), ix(sz), pi(sz),px(sz);
	
	ivec pindx(pcni.memptr(), np, true, true);

	for(int i=0; i<nit; i++)
	{
		const ivec rsp(dat.colptr(i), np, false, true);
		int indx = icnp[i];

		for(int p=0; p<np; p++)
		{
			if(rsp[p]>=0)
			{
				ip[indx] = p;
				ix[indx++] = rsp[p];
				pi[pindx[p]] = i;
				px[pindx[p]++] = rsp[p]; 
			}
		}
	}
	
	ivec imax(nit), isum(nit,fill::zeros), ncat(nit, fill::zeros);
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
		Named("icat") = icat, Named("ncat") = ncat, Named("imax") = imax, Named("isum") = isum, Named("psum") = psum);
}



// ix and px changed in place if necessary
// returns matrix a
// [[Rcpp::export]]
arma::imat categorize(const arma::ivec& inp, const arma::ivec& pni,
						const arma::ivec& icnp, const arma::ivec& pcni,
						const arma::ivec& ip, const arma::ivec& pi,				
						const arma::imat& icat, const arma::ivec& imax,
						const int max_cat,
						arma::ivec& ix, arma::ivec& px)
{
	
	const int nit = icat.n_cols, np = pni.n_elem;
	imat a(max_cat, nit, fill::zeros);
	imat ai(icat.n_rows, icat.n_cols,fill::zeros);
		
	bool recode = false;
	for(int i=0; i<nit; i++)
	{
		int k=0;
		for(int j=0; j<=imax[i]; j++)
		{
			if(icat.at(j,i) == 0) 
				recode = true;
			else
			{
				a.at(k,i) = j;
				ai.at(j,i) = k;
				k++;
			}
		}			
	}
	
	if(recode)
	{	
#pragma omp parallel for	
		for(int i=0; i<nit; i++)
		{
			for(int j=icnp[i]; j<icnp[i+1]; j++)
				ix[j] = ai.at(ix[j],i);
		}
#pragma omp parallel for	
		for(int p=0; p<np; p++)
		{
			for(int j=pcni[p]; j<pcni[p+1]; j++)
				px[j] = ai.at(px[j],pi[j]);
		}
	}	
	return a;
}

// at the moment only for NRM parametrisation !!
// also not correct for poly yet
// [[Rcpp::export]]
void scale_b(arma::mat& b, const arma::ivec& ncat, const arma::ivec& item_fixed)
{
	double mean_free=0, mean_fixed=0;
	const int nit = b.n_cols, n_fixed = accu(item_fixed);
	const int n_free =  nit - n_fixed;
	
	for(int i=0; i<nit; i++)
		for(int j=1; j<ncat[i]; j++)
		{
			mean_free += b.at(i,j) * (1-item_fixed[i]);
			mean_fixed += b.at(i,j) * item_fixed[i];	
		}
	mean_free /= n_free;
	mean_fixed /= n_fixed;
	double adj = mean_fixed - mean_free;
	
	for(int i=0; i<nit; i++)
		if(item_fixed[i]==0)
			for(int j=1; j<ncat[i]; j++)
				b.at(i,j) += adj;
}


/****************************
* Designs
*****************************/

// [[Rcpp::export]]
Rcpp::List design_matrices(const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& pg, const int nit, const int ng)
{
	const int np = pni.n_elem;
	imat item(nit,nit,fill::zeros), group(ng,nit,fill::zeros);
	
	std::vector<bool> bk(nit);
	std::unordered_map<std::vector<bool>, int> booklets;
	
	for(int p=0; p<np; p++)
	{
		std::fill(bk.begin(), bk.end(), false);
		int g = pg[p];
		for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
		{
			bk[pi[indx]] = true;
			group.at(g,pi[indx]) = 1;
		}
		booklets.insert(std::make_pair(bk,p));		
	}
	
	for(auto& iter: booklets )
	{
		int p = iter.second;
		for(int i= pcni[p]; i<pcni[p+1]; i++)
			for(int j=i+1; j<pcni[p+1]; j++)
				item.at(pi[i],pi[j]) = 1;	
	}
	item += item.t();
	item.diag().ones(); 
	return Rcpp::List::create(Named("items")=item, Named("groups")=group);
}

// bitflag, 0=unidentified; 1=connected; 2,4,8 and combinations are tenously identified (should give warnings)
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

